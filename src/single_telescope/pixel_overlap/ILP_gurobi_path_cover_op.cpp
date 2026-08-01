#include "ILP_gurobi_path_cover_op.h"

#include "gurobi_c++.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

// Contract shared with buildGraph/buildGraphOrienteering:
//   * costs already contains slew, acceleration/deceleration, settle, dwell,
//     and the dummy-end padding;
//   * start and end are zero-prize dummy nodes;
//   * Budget already includes the padding;
//   * coverage uses compact HEALPixel ids.
struct STPixelProblem {
    const std::vector<std::vector<double>>& costs;
    const std::vector<double>& pixelProbs;
    const std::vector<std::vector<int>>& coverage;
    int start;
    int end;
    double budget;

    int nodeCount() const {
        return static_cast<int>(costs.size());
    }

    int pixelCount() const {
        return static_cast<int>(pixelProbs.size());
    }
};

void validateProblem(const STPixelProblem& problem) {
    const int n = problem.nodeCount();
    const int m = problem.pixelCount();

    if (n < 2) {
        throw std::invalid_argument(
            "gurobiPathCover: the s-t graph must contain at least two nodes.");
    }
    if (problem.start < 0 || problem.start >= n ||
        problem.end < 0 || problem.end >= n ||
        problem.start == problem.end) {
        throw std::invalid_argument(
            "gurobiPathCover: invalid or identical start/end nodes.");
    }
    if (!std::isfinite(problem.budget) || problem.budget < 0.0) {
        throw std::invalid_argument(
            "gurobiPathCover: budget must be finite and nonnegative.");
    }
    if (static_cast<int>(problem.coverage.size()) != n) {
        throw std::invalid_argument(
            "gurobiPathCover: coverage must contain one list per graph node.");
    }

    for (int i = 0; i < n; ++i) {
        if (static_cast<int>(problem.costs[i].size()) != n) {
            throw std::invalid_argument(
                "gurobiPathCover: cost matrix must be square.");
        }
        for (int j = 0; j < n; ++j) {
            if (!std::isfinite(problem.costs[i][j]) ||
                problem.costs[i][j] < 0.0) {
                throw std::invalid_argument(
                    "gurobiPathCover: costs must be finite and nonnegative.");
            }
        }

        for (int pixel : problem.coverage[i]) {
            if (pixel < 0 || pixel >= m) {
                throw std::out_of_range(
                    "gurobiPathCover: coverage contains an invalid pixel id.");
            }
        }
    }

    for (double probability : problem.pixelProbs) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "gurobiPathCover: pixel probabilities must be finite and "
                "nonnegative.");
        }
    }

    if (!problem.coverage[problem.start].empty() ||
        !problem.coverage[problem.end].empty()) {
        throw std::invalid_argument(
            "gurobiPathCover: dummy start/end coverage lists must be empty.");
    }
}

bool isUsableArc(
    const STPixelProblem& problem,
    int from,
    int to) {
    if (from == to || from == problem.end || to == problem.start) {
        return false;
    }
    return problem.costs[from][to] <= problem.budget + 1e-12;
}

bool isValidInitialRoute(
    const STPixelProblem& problem,
    const std::vector<int>& route,
    std::string& reason) {
    if (route.empty()) {
        return false;
    }
    if (route.front() != problem.start || route.back() != problem.end) {
        reason = "route must begin at start and end at the dummy end";
        return false;
    }

    const int n = problem.nodeCount();
    std::vector<char> seen(n, 0);
    double totalCost = 0.0;

    for (std::size_t k = 0; k < route.size(); ++k) {
        const int node = route[k];
        if (node < 0 || node >= n || seen[node]) {
            reason = "route contains an invalid or repeated node";
            return false;
        }
        seen[node] = 1;

        if (k + 1 < route.size()) {
            const int next = route[k + 1];
            if (next < 0 || next >= n ||
                !isUsableArc(problem, node, next)) {
                reason = "route contains an invalid or over-budget arc";
                return false;
            }
            totalCost += problem.costs[node][next];
        }
    }

    if (totalCost > problem.budget + 1e-9) {
        reason = "route exceeds the already padded budget";
        return false;
    }
    return true;
}

void applyWarmStartCommon(
    const std::vector<int>& route,
    int n,
    std::vector<GRBVar>& visit,
    std::vector<std::vector<GRBVar>>& arc ) {
    std::vector<char> inRoute(n, 0);
    for (int node : route) {
        inRoute[node] = 1;
    }

    for (int i = 0; i < n; ++i) {
        visit[i].set(GRB_DoubleAttr_Start, inRoute[i] ? 1.0 : 0.0);
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                arc[i][j].set(GRB_DoubleAttr_Start, 0.0);
            }
        }
    }

    for (std::size_t k = 0; k + 1 < route.size(); ++k) {
        arc[route[k]][route[k + 1]].set(GRB_DoubleAttr_Start, 1.0);
    }
}

void applyWarmStartMTZ(
    const std::vector<int>& route,
    int n,
    std::vector<GRBVar>& visit,
    std::vector<std::vector<GRBVar>>& arc,
    std::vector<GRBVar>& order ) {
    applyWarmStartCommon(route, n, visit, arc);

    for (int i = 0; i < n; ++i) {
        order[i].set(GRB_DoubleAttr_Start, 0.0);
    }
    for (std::size_t k = 0; k < route.size(); ++k) {
        order[route[k]].set(
            GRB_DoubleAttr_Start, static_cast<double>(k));
    }
}

void applyWarmStartContinuousTime(
    const STPixelProblem& problem,
    const std::vector<int>& route,
    std::vector<GRBVar>& visit,
    std::vector<std::vector<GRBVar>>& arc,
    std::vector<GRBVar>& order,
    std::vector<GRBVar>& completion ) {
    const int n = problem.nodeCount();
    applyWarmStartCommon(route, n, visit, arc);

    for (int i = 0; i < n; ++i) {
        order[i].set(GRB_DoubleAttr_Start, 0.0);
        completion[i].set(GRB_DoubleAttr_Start, 0.0);
    }

    double elapsed = 0.0;
    for (std::size_t k = 0; k < route.size(); ++k) {
        const int node = route[k];
        order[node].set(
            GRB_DoubleAttr_Start, static_cast<double>(k));
        completion[node].set(GRB_DoubleAttr_Start, elapsed);

        if (k + 1 < route.size()) {
            elapsed += problem.costs[node][route[k + 1]];
        }
    }
}

void addCoverageConstraints(
    GRBModel& model,
    const STPixelProblem& problem,
    const std::vector<GRBVar>& visit,
    const std::vector<GRBVar>& pixelCovered ) {
    const int m = problem.pixelCount();
    std::vector<GRBLinExpr> coverExpression(m, 0.0);

    for (int node = 0; node < problem.nodeCount(); ++node) {
        std::unordered_set<int> pixelsSeenInNode;
        for (int pixel : problem.coverage[node]) {
            if (pixelsSeenInNode.insert(pixel).second) {
                coverExpression[pixel] += visit[node];
            }
        }
    }

    // There is exactly one variable per compact HEALPixel. Multiple selected
    // tiles can raise the right-hand side, but the [0,1] pixel variable and
    // objective count its probability only once.
    for (int pixel = 0; pixel < m; ++pixel) {
        model.addConstr(
            pixelCovered[pixel] <= coverExpression[pixel],
            "pixel_covered_" + std::to_string(pixel));
    }
}

void addSTPathConstraints(
    GRBModel& model,
    const STPixelProblem& problem,
    const std::vector<GRBVar>& visit,
    const std::vector<std::vector<GRBVar>>& arc ) {
    const int n = problem.nodeCount();
    const int start = problem.start;
    const int end = problem.end;

    model.addConstr(visit[start] == 1.0, "start_active");
    model.addConstr(visit[end] == 1.0, "end_active");

    GRBLinExpr startIn = 0.0;
    GRBLinExpr startOut = 0.0;
    GRBLinExpr endIn = 0.0;
    GRBLinExpr endOut = 0.0;
    for (int node = 0; node < n; ++node) {
        if (node != start) {
            startIn += arc[node][start];
            startOut += arc[start][node];
        }
        if (node != end) {
            endIn += arc[node][end];
            endOut += arc[end][node];
        }
    }
    model.addConstr(startIn == 0.0, "start_in_zero");
    model.addConstr(startOut == 1.0, "start_out_one");
    model.addConstr(endIn == 1.0, "end_in_one");
    model.addConstr(endOut == 0.0, "end_out_zero");

    for (int node = 0; node < n; ++node) {
        if (node == start || node == end) {
            continue;
        }

        GRBLinExpr incoming = 0.0;
        GRBLinExpr outgoing = 0.0;
        for (int other = 0; other < n; ++other) {
            if (other != node) {
                incoming += arc[other][node];
                outgoing += arc[node][other];
            }
        }
        model.addConstr(
            incoming == visit[node],
            "node_in_" + std::to_string(node));
        model.addConstr(
            outgoing == visit[node],
            "node_out_" + std::to_string(node));
    }
}

void addOrderingConstraints(
    GRBModel& model,
    const STPixelProblem& problem,
    const std::vector<GRBVar>& visit,
    const std::vector<std::vector<GRBVar>>& arc,
    const std::vector<GRBVar>& order,
    const std::string& prefix ) {
    const int n = problem.nodeCount();
    model.addConstr(order[problem.start] == 0.0, prefix + "_start");

    for (int node = 0; node < n; ++node) {
        if (node == problem.start) {
            continue;
        }
        model.addConstr(
            order[node] <= (n - 1) * visit[node],
            prefix + "_upper_" + std::to_string(node));
        model.addConstr(
            order[node] >= visit[node],
            prefix + "_lower_" + std::to_string(node));
    }

    const double bigM = static_cast<double>(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!isUsableArc(problem, i, j)) {
                continue;
            }
            model.addConstr(
                order[j] >= order[i] + 1.0 - bigM * (1.0 - arc[i][j]),
                prefix + "_arc_" + std::to_string(i) + "_" +
                    std::to_string(j));
        }
    }
}

void addBudgetConstraint(
    GRBModel& model,
    const STPixelProblem& problem,
    const std::vector<std::vector<GRBVar>>& arc ) {
    GRBLinExpr totalCost = 0.0;
    for (int i = 0; i < problem.nodeCount(); ++i) {
        for (int j = 0; j < problem.nodeCount(); ++j) {
            if (isUsableArc(problem, i, j)) {
                totalCost += problem.costs[i][j] * arc[i][j];
            }
        }
    }
    model.addConstr(totalCost <= problem.budget, "budget");
}

std::vector<int> extractSTPath(
    const STPixelProblem& problem,
    const std::vector<std::vector<GRBVar>>& arc ) {
    const int n = problem.nodeCount();
    std::vector<int> path;
    std::vector<char> seen(n, 0);
    int current = problem.start;

    while (current >= 0 && current < n && !seen[current]) {
        path.push_back(current);
        seen[current] = 1;
        if (current == problem.end) {
            return path;
        }

        int successor = -1;
        double largestValue = 0.5;
        for (int next = 0; next < n; ++next) {
            if (next == current) {
                continue;
            }
            const double value = arc[current][next].get(GRB_DoubleAttr_X);
            if (value > largestValue) {
                largestValue = value;
                successor = next;
            }
        }
        if (successor < 0) {
            break;
        }
        current = successor;
    }

    std::cerr << "gurobiPathCover: failed to extract a complete s-t path.\n";
    return {};
}

std::vector<int> optimizeAndExtract(
    GRBModel& model,
    const STPixelProblem& problem,
    const std::vector<std::vector<GRBVar>>& arc,
    const std::string& modelName ) {
    try {
        model.optimize();
        if (model.get(GRB_IntAttr_SolCount) <= 0) {
            std::cerr << modelName << ": no incumbent solution.\n";
            return {};
        }

        const int status = model.get(GRB_IntAttr_Status);
        const double incumbent = model.get(GRB_DoubleAttr_ObjVal);
        const double bestBound = model.get(GRB_DoubleAttr_ObjBound);
        const double gap = model.get(GRB_DoubleAttr_MIPGap);
        std::cout << modelName
                  << ": incumbent=" << incumbent
                  << ", best_bound=" << bestBound
                  << ", gap=" << gap
                  << ", status=" << status << "\n";

        return extractSTPath(problem, arc);
    } catch (GRBException& error) {
        std::cerr << modelName << ": Gurobi error "
                  << error.getErrorCode() << ": "
                  << error.getMessage() << "\n";
        return {};
    }
}

std::vector<int> solveMTZST(
    const STPixelProblem& problem,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute ) {
    validateProblem(problem);

    GRBEnv environment(true);
    environment.start();
    GRBModel model(environment);
    model.set(GRB_StringAttr_ModelName, "ST_Pixel_Path_Cover_MTZ");

    const int n = problem.nodeCount();
    const int m = problem.pixelCount();

    // Only routing arcs are binary. visit is determined by incoming/outgoing
    // binary arcs, and pixelCovered is driven to 1 by its nonnegative prize.
    std::vector<GRBVar> visit(n);
    for (int node = 0; node < n; ++node) {
        visit[node] = model.addVar(
            0.0, 1.0, 0.0, GRB_CONTINUOUS,
            "visit_" + std::to_string(node));
    }

    std::vector<GRBVar> pixelCovered(m);
    for (int pixel = 0; pixel < m; ++pixel) {
        pixelCovered[pixel] = model.addVar(
            0.0, 1.0, 0.0, GRB_CONTINUOUS,
            "pixel_" + std::to_string(pixel));
    }

    std::vector<std::vector<GRBVar>> arc(
        n, std::vector<GRBVar>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            arc[i][j] = model.addVar(
                0.0, isUsableArc(problem, i, j) ? 1.0 : 0.0,
                0.0, GRB_BINARY,
                "arc_" + std::to_string(i) + "_" + std::to_string(j));
        }
    }

    std::vector<GRBVar> order(n);
    for (int node = 0; node < n; ++node) {
        order[node] = model.addVar(
            0.0, static_cast<double>(n - 1),
            0.0, GRB_CONTINUOUS,
            "order_" + std::to_string(node));
    }

    model.update();

    GRBLinExpr objective = 0.0;
    for (int pixel = 0; pixel < m; ++pixel) {
        objective += problem.pixelProbs[pixel] * pixelCovered[pixel];
    }
    model.setObjective(objective, GRB_MAXIMIZE);

    addCoverageConstraints(model, problem, visit, pixelCovered);
    addSTPathConstraints(model, problem, visit, arc);
    addOrderingConstraints(model, problem, visit, arc, order, "mtz");
    addBudgetConstraint(model, problem, arc);

    if (mipGap >= 0.0) {
        model.set(GRB_DoubleParam_MIPGap, mipGap);
    }
    if (timeLimit > 0.0) {
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    }

    if (!initialRoute.empty()) {
        std::string reason;
        if (isValidInitialRoute(problem, initialRoute, reason)) {
            applyWarmStartMTZ(initialRoute, n, visit, arc, order);
        } else {
            std::cerr << "MTZ warm start skipped: " << reason << ".\n";
        }
    }

    return optimizeAndExtract(model, problem, arc, "MTZ path-cover model");
}

std::vector<int> solveContinuousTimeST(
    const STPixelProblem& problem,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute ) {
    validateProblem(problem);

    GRBEnv environment(true);
    environment.start();
    GRBModel model(environment);
    model.set(
        GRB_StringAttr_ModelName,
        "ST_Pixel_Path_Cover_Continuous_Time");

    const int n = problem.nodeCount();
    const int m = problem.pixelCount();

    std::vector<GRBVar> visit(n);
    for (int node = 0; node < n; ++node) {
        visit[node] = model.addVar(
            0.0, 1.0, 0.0, GRB_CONTINUOUS,
            "visit_" + std::to_string(node));
    }

    std::vector<GRBVar> pixelCovered(m);
    for (int pixel = 0; pixel < m; ++pixel) {
        pixelCovered[pixel] = model.addVar(
            0.0, 1.0, 0.0, GRB_CONTINUOUS,
            "pixel_" + std::to_string(pixel));
    }

    std::vector<std::vector<GRBVar>> arc(
        n, std::vector<GRBVar>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            arc[i][j] = model.addVar(
                0.0, isUsableArc(problem, i, j) ? 1.0 : 0.0,
                0.0, GRB_BINARY,
                "arc_" + std::to_string(i) + "_" + std::to_string(j));
        }
    }

    // order guarantees connectivity even if a caller supplies a zero-cost
    // off-diagonal edge; completion carries the graph's cumulative cost.
    std::vector<GRBVar> order(n);
    std::vector<GRBVar> completion(n);
    for (int node = 0; node < n; ++node) {
        order[node] = model.addVar(
            0.0, static_cast<double>(n - 1),
            0.0, GRB_CONTINUOUS,
            "order_" + std::to_string(node));
        completion[node] = model.addVar(
            0.0, problem.budget,
            0.0, GRB_CONTINUOUS,
            "completion_" + std::to_string(node));
    }

    model.update();

    GRBLinExpr objective = 0.0;
    for (int pixel = 0; pixel < m; ++pixel) {
        objective += problem.pixelProbs[pixel] * pixelCovered[pixel];
    }
    model.setObjective(objective, GRB_MAXIMIZE);

    addCoverageConstraints(model, problem, visit, pixelCovered);
    addSTPathConstraints(model, problem, visit, arc);
    addOrderingConstraints(
        model, problem, visit, arc, order, "time_order");
    addBudgetConstraint(model, problem, arc);

    model.addConstr(
        completion[problem.start] == 0.0,
        "completion_start");
    for (int node = 0; node < n; ++node) {
        model.addConstr(
            completion[node] <= problem.budget * visit[node],
            "completion_active_" + std::to_string(node));
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (!isUsableArc(problem, i, j)) {
                continue;
            }
            const double cost = problem.costs[i][j];
            const double bigM = problem.budget + cost;
            model.addConstr(
                completion[j] >= completion[i] + cost -
                    bigM * (1.0 - arc[i][j]),
                "time_" + std::to_string(i) + "_" +
                    std::to_string(j));
        }
    }

    if (mipGap >= 0.0) {
        model.set(GRB_DoubleParam_MIPGap, mipGap);
    }
    if (timeLimit > 0.0) {
        model.set(GRB_DoubleParam_TimeLimit, timeLimit);
    }

    if (!initialRoute.empty()) {
        std::string reason;
        if (isValidInitialRoute(problem, initialRoute, reason)) {
            applyWarmStartContinuousTime(
                problem, initialRoute, visit, arc, order, completion);
        } else {
            std::cerr << "Continuous-time warm start skipped: "
                      << reason << ".\n";
        }
    }

    return optimizeAndExtract(
        model, problem, arc, "Continuous-time path-cover model");
}

}  // namespace

// MTZ entry point. costs is used directly; do not pass a pure slew matrix and
// do not add dwell or padding again inside the solver.
std::vector<int> gurobiPathCoverPixelOverlapMTZ(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double Budget,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute ) {
    try {
        const STPixelProblem problem{
            costs, pixel_probs, member_pixels, start, end, Budget};
        return solveMTZST(problem, mipGap, timeLimit, initialRoute);
    } catch (GRBException& error) {
        std::cerr << "MTZ path-cover model: Gurobi error "
                  << error.getErrorCode() << ": "
                  << error.getMessage() << "\n";
        return {};
    }
}


// Continuous-time entry point with the same graph contract.
std::vector<int> gurobiPathCoverPixelOverlapContinuousTime(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double Budget,
    double mipGap,
    double timeLimit,
    const std::vector<int>& initialRoute ) {
    try {
        const STPixelProblem problem{
            costs, pixel_probs, member_pixels, start, end, Budget};
        return solveContinuousTimeST(
            problem, mipGap, timeLimit, initialRoute);
    } catch (GRBException& error) {
        std::cerr << "Continuous-time path-cover model: Gurobi error "
                  << error.getErrorCode() << ": "
                  << error.getMessage() << "\n";
        return {};
    }
}