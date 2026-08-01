#include "ILP_gurobi_path_cover.h"

#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <stdexcept>
#include <limits>
#include <cmath>

static void applyWarmStartPath_Common(
    const std::vector<int>& route,
    int n,
    int start,
    std::vector<GRBVar>& x,
    std::vector<std::vector<GRBVar>>& z ) {
    if (route.empty()) return;

    // Route should start with `start` (if not, we still try but warn by behavior)
    // (You can enforce strictly if you want.)
    std::vector<char> inRoute(n, 0);
    for (int v : route) {
        if (v < 0 || v >= n) continue;
        inRoute[v] = 1;
    }

    // x warm start
    for (int i = 0; i < n; ++i) {
        x[i].set(GRB_DoubleAttr_Start, inRoute[i] ? 1.0 : 0.0);
    }
    // enforce start visited in warm start
    if (start >= 0 && start < n) x[start].set(GRB_DoubleAttr_Start, 1.0);

    // z warm start: set arcs along consecutive route pairs
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            if (z[i][j].get(GRB_CharAttr_VType) == GRB_BINARY) {
                z[i][j].set(GRB_DoubleAttr_Start, 0.0);
            }
        }
    }
    for (size_t k = 0; k + 1 < route.size(); ++k) {
        int i = route[k], j = route[k + 1];
        if (i < 0 || i >= n || j < 0 || j >= n || i == j) continue;
        z[i][j].set(GRB_DoubleAttr_Start, 1.0);
        x[i].set(GRB_DoubleAttr_Start, 1.0);
        x[j].set(GRB_DoubleAttr_Start, 1.0);
    }
}

static void applyWarmStartPath_ContinuousTime(
    const std::vector<int>& route,
    int n,
    int start,
    double D,
    const std::vector<double>& dwell,
    const std::vector<std::vector<double>>& slew,
    std::vector<GRBVar>& x,
    std::vector<std::vector<GRBVar>>& z,
    std::vector<GRBVar>& tau ) {
    if (route.empty()) return;

    applyWarmStartPath_Common(route, n, start, x, z);

    // tau warm start (cumulative time along the route)
    std::vector<double> t(n, 0.0);
    std::vector<char> hasTau(n, 0);

    if (start >= 0 && start < n) {
        t[start] = 0.0;
        hasTau[start] = 1;
        tau[start].set(GRB_DoubleAttr_Start, 0.0);
    }

    double cur = 0.0;
    for (size_t k = 0; k + 1 < route.size(); ++k) {
        int i = route[k], j = route[k + 1];
        if (i < 0 || i >= n || j < 0 || j >= n || i == j) continue;
        // start time at j = start time at i + dwell_i + slew_ij
        // (this is consistent with your time-flow constraint)
        double Tij = dwell[i] + slew[i][j];
        cur = cur + Tij;
        if (cur > D) cur = D; // keep within bounds for warm start feasibility-ish
        if (!hasTau[j]) {
            t[j] = cur;
            hasTau[j] = 1;
            tau[j].set(GRB_DoubleAttr_Start, t[j]);
        }
    }

    // For nodes not in route, give a harmless start value
    for (int i = 0; i < n; ++i) {
        if (!hasTau[i]) tau[i].set(GRB_DoubleAttr_Start, 0.0);
    }
}

static void applyWarmStartPath_MTZ(
    const std::vector<int>& route,
    int n,
    int start,
    std::vector<GRBVar>& x,
    std::vector<std::vector<GRBVar>>& z,
    std::vector<GRBVar>& u ) {
    if (route.empty()) return;

    applyWarmStartPath_Common(route, n, start, x, z);

    // u warm start: order along the route
    for (int i = 0; i < n; ++i) u[i].set(GRB_DoubleAttr_Start, 0.0);
    for (size_t k = 0; k < route.size(); ++k) {
        int v = route[k];
        if (v < 0 || v >= n) continue;
        u[v].set(GRB_DoubleAttr_Start, (double)k);
    }
    if (start >= 0 && start < n) u[start].set(GRB_DoubleAttr_Start, 0.0);
}

static std::vector<int> extractPathFromZ(
    int n, int start,
    const std::vector<std::vector<GRBVar>>& z,
    double eps = 0.5 ) {
    std::vector<int> succ(n, -1);

    for (int i = 0; i < n; ++i) {
        int bestj = -1;
        double bestv = 0.0;
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            // z[i][j] exists only if i!=j in your code
            double val = z[i][j].get(GRB_DoubleAttr_X);
            if (val > bestv) { bestv = val; bestj = j; }
        }
        if (bestv >= eps) succ[i] = bestj;
    }

    std::vector<int> path;
    std::vector<char> seen(n, 0);

    int cur = start;
    while (cur >= 0 && cur < n && !seen[cur]) {
        path.push_back(cur);
        seen[cur] = 1;
        int nxt = succ[cur];
        if (nxt == -1) break;
        cur = nxt;
    }
    return path;
}


std::vector<int> solveContinuousTimeModel(
    const ProblemInstance& inst,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute // can be empty
) {

    for (std::size_t i = 0; i < inst.dwell_times.size(); ++i) {
    if (!std::isfinite(inst.dwell_times[i])) {
        std::cerr << "Bad dwell_times["<<i<<"]="<<inst.dwell_times[i]<<"\n";
        break;
    }
    for (std::size_t j = 0; j < inst.slew_times.size(); ++j) {
        if (i==j) continue;
        if (!std::isfinite(inst.slew_times[i][j])) {
        std::cerr << "Bad slew_times["<<i<<"]["<<j<<"]="<<inst.slew_times[i][j]<<"\n";
        break;
        }
    }
    }

    
    GRBEnv env(true);
    env.start();
    GRBModel model(env);
    // model.set(GRB_DoubleParam_OptimalityTol, 1e-7);
    // model.set(GRB_DoubleParam_ObjScale, 0.01);    

    int n = inst.nTiles;
    int m = inst.nPixels;
    int start = inst.start;
    double D = inst.Deadline;

    // try {

    // --- Vars ---
    std::vector<GRBVar> x(n);
    for (int i = 0; i < n; ++i)
        x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(i));

    std::vector<GRBVar> y(m);
    for (int p = 0; p < m; ++p)
        y[p] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + std::to_string(p));

    std::vector<std::vector<GRBVar>> z(n, std::vector<GRBVar>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                z[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                       "z_" + std::to_string(i) + "_" + std::to_string(j));

    std::vector<GRBVar> tau(n);
    for (int i = 0; i < n; ++i)
        tau[i] = model.addVar(0.0, D, 0.0, GRB_CONTINUOUS, "tau_" + std::to_string(i));

    model.update();

    // --- Objective ---
    GRBLinExpr obj = 0.0;
    for (int p = 0; p < m; ++p) obj += inst.prizes[p] * y[p];
    model.setObjective(obj, GRB_MAXIMIZE);

    // --- Coverage y[p] <= sum_i x[i] ---
    std::vector<GRBLinExpr> coverExpr(m, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int p : inst.coverage[i]) {
            if (0 <= p && p < m) coverExpr[p] += x[i];
        }
    }
    for (int p = 0; p < m; ++p)
        model.addConstr(y[p] <= coverExpr[p], "cover_once_p" + std::to_string(p));

    // --- Start ---
    model.addConstr(x[start] == 1, "start_visited");
    for (int i = 0; i < n; ++i) {
        model.addConstr(tau[i] <= D * x[i], "tau_active_" + std::to_string(i));
    }
    model.addConstr(tau[start] == 0.0, "start_time");  // keep this

    // model.addConstr(tau[start] == 0.0, "start_time");

    // start out <= 1
    {
        GRBLinExpr out = 0.0;
        for (int j = 0; j < n; ++j) if (j != start) out += z[start][j];
        model.addConstr(out <= 1.0, "start_out");
    }

    // non-start nodes: exactly one predecessor iff visited
    for (int k = 0; k < n; ++k) {
        if (k == start) continue;
        GRBLinExpr in = 0.0;
        for (int i = 0; i < n; ++i) if (i != k) in += z[i][k];
        model.addConstr(in == x[k], "one_pred_" + std::to_string(k));
    }

    // depart only if visited (out <= x)
    for (int i = 0; i < n; ++i) {
        GRBLinExpr out = 0.0;
        for (int j = 0; j < n; ++j) if (i != j) out += z[i][j];
        model.addConstr(out <= x[i], "depart_visit_" + std::to_string(i));
    }

    // time flow
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            double Tij = inst.dwell_times[i] + inst.slew_times[i][j];
            double Mij  = D + Tij;
            model.addConstr(tau[j] >= tau[i] + Tij - Mij * (1.0 - z[i][j]),
                            "time_flow_" + std::to_string(i) + "_" + std::to_string(j));

        }
    }

    // deadline per tile
    for (int i = 0; i < n; ++i) {
        // tau[i] is zero when x[i] is zero, so this activates the dwell
        // deadline exactly when tile i is selected, without another big-M.
        model.addConstr(tau[i] + inst.dwell_times[i] * x[i] <= D,
                        "deadline_" + std::to_string(i));
    }

    // --- Params ---
    model.set(GRB_DoubleParam_MIPGap, mipGap);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    
    // --- Warm start ---
    if (!initialRoute.empty()) {
        applyWarmStartPath_ContinuousTime(initialRoute, n, start, D,
                                          inst.dwell_times, inst.slew_times, x, z, tau);
    }

    // --- Solve ---
    // model.optimize();
    try {
        model.optimize();
    } catch (GRBException& e) {
        std::cerr << "Gurobi error code = " << e.getErrorCode() << "\n";
        std::cerr << e.getMessage() << "\n";
        model.write("continuous_time_fail.lp");
        return {};
    }


    int status = model.get(GRB_IntAttr_Status);
    if (!(status == GRB_OPTIMAL || status == GRB_TIME_LIMIT || status == GRB_SUBOPTIMAL)) {
        // Return empty if infeasible / etc.
        return {};
    }
    if (model.get(GRB_IntAttr_SolCount) == 0) {
        // A time limit does not guarantee that an incumbent was found.
        return {};
    }

    return extractPathFromZ(n, start, z);
    // } catch (GRBException& e) {
    // std::cerr << "Gurobi error code = " << e.getErrorCode() << "\n";
    // std::cerr << e.getMessage() << "\n";
    // try { model.write("continuous_time_fail.lp"); } catch (...) {}
    // return {};
    // }
}


std::vector<int> solveMTZModel(
    const ProblemInstance& inst,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute
) {
    GRBEnv env(true);
    env.start();
    GRBModel model(env);
    model.set(GRB_DoubleParam_OptimalityTol, 1e-7);
    model.set(GRB_DoubleParam_ObjScale, 0.01);

    int n = inst.nTiles;
    int m = inst.nPixels;
    int start = inst.start;
    double D = inst.Deadline;

    // --- Vars ---
    std::vector<GRBVar> x(n);
    for (int i = 0; i < n; ++i)
        x[i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "x_" + std::to_string(i));

    std::vector<GRBVar> y(m);
    for (int p = 0; p < m; ++p)
        y[p] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, "y_" + std::to_string(p));

    std::vector<std::vector<GRBVar>> z(n, std::vector<GRBVar>(n));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            if (i != j)
                z[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                       "z_" + std::to_string(i) + "_" + std::to_string(j));

    std::vector<GRBVar> u(n);
    for (int i = 0; i < n; ++i)
        u[i] = model.addVar(0.0, n, 0.0, GRB_CONTINUOUS, "u_" + std::to_string(i));

    model.update();

    // --- Objective ---
    GRBLinExpr obj = 0.0;
    for (int p = 0; p < m; ++p) obj += inst.prizes[p] * y[p];
    model.setObjective(obj, GRB_MAXIMIZE);

    // --- Coverage ---
    std::vector<GRBLinExpr> coverExpr(m, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int p : inst.coverage[i]) {
            if (0 <= p && p < m) coverExpr[p] += x[i];
        }
    }
    for (int p = 0; p < m; ++p)
        model.addConstr(y[p] <= coverExpr[p], "cover_once_p" + std::to_string(p));

    // --- Budget ---
    {
        GRBLinExpr totalTime = 0.0;
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                if (i != j)
                    totalTime += inst.slew_times[i][j] * z[i][j];
        for (int i = 0; i < n; ++i)
            totalTime += inst.dwell_times[i] * x[i];
        model.addConstr(totalTime <= D, "budget");
    }

    // --- Start ---
    model.addConstr(x[start] == 1, "start_visited");

    // start out <= 1
    {
        GRBLinExpr out = 0.0;
        for (int j = 0; j < n; ++j) if (j != start) out += z[start][j];
        model.addConstr(out <= 1.0, "start_out");
    }

    // non-start: exactly one predecessor iff visited  (IMPORTANT: avoid isolated visited tiles)
    for (int k = 0; k < n; ++k) {
        if (k == start) continue;
        GRBLinExpr in = 0.0;
        for (int i = 0; i < n; ++i) if (i != k) in += z[i][k];
        model.addConstr(in == x[k], "one_pred_" + std::to_string(k));
    }

    // depart only if visited (out <= x)
    for (int i = 0; i < n; ++i) {
        GRBLinExpr out = 0.0;
        for (int j = 0; j < n; ++j) if (i != j) out += z[i][j];
        model.addConstr(out <= x[i], "depart_visit_" + std::to_string(i));
    }

    // --- MTZ ---
    model.addConstr(u[start] == 0.0, "u_start");
    for (int i = 0; i < n; ++i) {
        if (i == start) continue;
        model.addConstr(u[i] <= n * x[i], "u_ub_" + std::to_string(i));
        model.addConstr(u[i] >= x[i],     "u_lb_" + std::to_string(i));
    }

    double Mmtz = n + 1.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) continue;
            model.addConstr(u[j] >= u[i] + 1.0 - Mmtz * (1.0 - z[i][j]),
                            "mtz_" + std::to_string(i) + "_" + std::to_string(j));
        }
    }

    // --- Params ---
    model.set(GRB_DoubleParam_MIPGap, mipGap);
    model.set(GRB_DoubleParam_TimeLimit, timeLimit);

    // --- Warm start ---
    if (!initialRoute.empty()) {
        applyWarmStartPath_MTZ(initialRoute, n, start, x, z, u);
    }

    // --- Solve ---
    try {
        model.optimize();
    } catch (GRBException& e) {
        std::cerr << "Gurobi error code = " << e.getErrorCode() << "\n";
        std::cerr << e.getMessage() << "\n";
        return {};
    }

    int status = model.get(GRB_IntAttr_Status);
    if (!(status == GRB_OPTIMAL || status == GRB_TIME_LIMIT || status == GRB_SUBOPTIMAL)) {
        double bestObj = model.get(GRB_DoubleAttr_ObjBound);
        return {};
    }
    if (model.get(GRB_IntAttr_SolCount) == 0) {
        return {};
    }

    return extractPathFromZ(n, start, z);
}


std::vector<int> gurobiPathCover(
    const std::vector<std::vector<double>>& slew_times,
    const std::vector<double>& dwell_times,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start, double Budget,
    double mipGap, double timeLimit,
    const std::vector<int>& initialRoute) {
        
    ProblemInstance inst(
        start,
        Budget,
        pixel_probs,
        dwell_times,
        slew_times,
        member_pixels
    );

    // return solveContinuousTimeModel(inst, mipGap, timeLimit, initialRoute);
    return solveMTZModel(inst, mipGap, timeLimit, initialRoute);
}