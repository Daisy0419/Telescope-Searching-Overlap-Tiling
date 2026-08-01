#include "ILP_gurobi_mink.h"
 
#include <gurobi_c++.h>
 
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <unordered_map>
#include <cmath>


//common start
std::pair<std::vector<std::vector<int>>, double> solveMinRoutesCover(const ProblemDataMinKCommonR& data,
                                double mipGap,
                                double timeLimit,
                                const std::vector<std::vector<int>>& initialRoutes) {
    const int N      = data.N;
    const int K_max  = data.K;
    const int s      = data.startIndex;
    const int t      = data.endIndex;
    const double Budget = data.BudgetPerTeam;
    const auto& C    = data.Cost;
 
    std::vector<std::vector<int>> routes;
    double bestObj = 0.0;
 
    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_min_routes.log");
        env.start();
 
        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Min_Routes_Cover");
 
        if (mipGap >= 0.0)
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);
 
        // MINIMIZE total number of active routes
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        // model.set(GRB_DoubleParam_ObjScale, 0.001);
 
        const int M = std::max(1, N - 1);  // flow capacity upper bound
 

        //  Variables
        // u[k] ∈ {0,1}: whether route k is active.  obj coeff = 1.
        std::vector<GRBVar> u(K_max);
        for (int k = 0; k < K_max; ++k)
            u[k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
                                "u_" + std::to_string(k));
 
        // x[k][i][j] ∈ {0,1}: edge (i,j) used by route k
        std::vector<std::vector<std::vector<GRBVar>>> x(
            K_max, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));
 
        // z[k][i] ∈ {0,1}: node i visited by route k
        std::vector<std::vector<GRBVar>> z(K_max, std::vector<GRBVar>(N));
 
        // f[k][i][j] ≥ 0: SCF flow on edge (i,j) in route k
        std::vector<std::vector<std::vector<GRBVar>>> f(
            K_max, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));
 
        for (int k = 0; k < K_max; ++k) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    x[k][i][j] = model.addVar(
                        0.0, 1.0, 0.0, GRB_BINARY,
                        "x_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                    f[k][i][j] = model.addVar(
                        0.0, M, 0.0, GRB_CONTINUOUS,
                        "f_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                }
            }
            for (int i = 0; i < N; ++i)
                z[k][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                    "z_" + std::to_string(k) + "_" + std::to_string(i));
        }
 
        //  Constraints
        for (int k = 0; k < K_max; ++k) {
            // --- (1a) Start node: out-degree = u[k], in-degree = 0 ---
            GRBLinExpr startOut = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) startOut += x[k][s][j];
            model.addConstr(startOut == u[k],
                "StartOut_k" + std::to_string(k));
 
            GRBLinExpr startIn = 0.0;
            for (int i = 0; i < N; ++i) if (i != s) startIn += x[k][i][s];
            model.addConstr(startIn == 0.0,
                "StartInZero_k" + std::to_string(k));
 
            // --- (1b) End node: in-degree = u[k], out-degree = 0 ---
            GRBLinExpr endIn = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) endIn += x[k][i][t];
            model.addConstr(endIn == u[k],
                "EndIn_k" + std::to_string(k));
 
            GRBLinExpr endOut = 0.0;
            for (int j = 0; j < N; ++j) if (j != t) endOut += x[k][t][j];
            model.addConstr(endOut == 0.0,
                "EndOutZero_k" + std::to_string(k));
 
            // --- (1c) Flow conservation at intermediate nodes ---
            for (int v = 0; v < N; ++v) {
                if (v == s || v == t) continue;
                GRBLinExpr out = 0.0, in = 0.0;
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.addConstr(out == in,
                    "DegBal_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
 
            // --- (1d) Link z to x ---
            for (int v = 0; v < N; ++v) {
                if (v == s || v == t) {
                    // s and t are visited iff route is active
                    model.addConstr(z[k][v] == u[k],
                        "zFixed_k" + std::to_string(k) + "_v" + std::to_string(v));
                } else {
                    GRBLinExpr out = 0.0;
                    for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                    model.addConstr(out == z[k][v],
                        "DefZ_k" + std::to_string(k) + "_v" + std::to_string(v));
                }
            }
        }
 
        // --- (2) Cover all customer nodes exactly once ---
        for (int v = 0; v < N; ++v) {
            if (v == s || v == t) continue;
            GRBLinExpr sumz = 0.0;
            for (int k = 0; k < K_max; ++k) sumz += z[k][v];
            model.addConstr(sumz == 1.0, "Cover_v" + std::to_string(v));
        }
 
        // --- (3) Per-route budget ---
        for (int k = 0; k < K_max; ++k) {
            GRBLinExpr cost = 0.0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    cost += C[i][j] * x[k][i][j];
            model.addConstr(cost <= Budget * u[k],
                "Budget_k" + std::to_string(k));
        }
 
        // --- (4) SCF subtour elimination ---
        for (int k = 0; k < K_max; ++k) {
            // Flow capacity
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    model.addConstr(f[k][i][j] <= M * x[k][i][j],
                        "FlowCap_k" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
 
            // Flow supply at s: net outflow = number of nodes visited (excl s)
            GRBLinExpr fOutS = 0.0, fInS = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) fOutS += f[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) fInS  += f[k][i][s];
 
            GRBLinExpr visitsExclS = 0.0;
            for (int v = 0; v < N; ++v) if (v != s) visitsExclS += z[k][v];
            model.addConstr(fOutS - fInS == visitsExclS,
                "FlowSupply_k" + std::to_string(k));
 
            // Flow demand at every non-s node
            for (int v = 0; v < N; ++v) {
                if (v == s) continue;
                GRBLinExpr fIn = 0.0, fOut = 0.0;
                for (int i = 0; i < N; ++i) if (i != v) fIn  += f[k][i][v];
                for (int j = 0; j < N; ++j) if (j != v) fOut += f[k][v][j];
                model.addConstr(fIn - fOut == z[k][v],
                    "FlowDemand_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
        }
 
        // --- (5) Symmetry breaking: u[0] ≥ u[1] ≥ … ≥ u[K_max-1] ---
        for (int k = 0; k + 1 < K_max; ++k)
            model.addConstr(u[k] >= u[k + 1],
                "SymBreak_" + std::to_string(k));

        //  MIP warm start from initialRoutes
        if (!initialRoutes.empty()) {
            // Zero out all variables first
            for (int k = 0; k < K_max; ++k) {
                u[k].set(GRB_DoubleAttr_Start, 0.0);
                for (int i = 0; i < N; ++i) {
                    z[k][i].set(GRB_DoubleAttr_Start, 0.0);
                    for (int j = 0; j < N; ++j) if (i != j) {
                        x[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                        f[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                    }
                }
            }
 
            int nInit = std::min(static_cast<int>(initialRoutes.size()), K_max);
            for (int k = 0; k < nInit; ++k) {
                const auto& route = initialRoutes[k];
                if (route.size() < 2) continue;
 
                u[k].set(GRB_DoubleAttr_Start, 1.0);
 
                // Set z for every node in the route
                for (int node : route)
                    z[k][node].set(GRB_DoubleAttr_Start, 1.0);
 
                // Set x for consecutive edges, and compute SCF flow values.
                // Flow value on edge (route[i], route[i+1]) = #nodes remaining
                // after route[i] (excluding s).
                int routeLen = static_cast<int>(route.size());
                for (int idx = 0; idx + 1 < routeLen; ++idx) {
                    int from = route[idx], to = route[idx + 1];
                    x[k][from][to].set(GRB_DoubleAttr_Start, 1.0);
 
                    // Number of nodes downstream of this edge (including 'to'),
                    // excluding s.
                    int nodesAhead = 0;
                    for (int r = idx + 1; r < routeLen; ++r)
                        if (route[r] != s) ++nodesAhead;
                    f[k][from][to].set(GRB_DoubleAttr_Start,
                                       static_cast<double>(nodesAhead));
                }
            }
        }
 

        //  Solve
        model.update();
        model.optimize();
 
        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            // Check that a feasible solution exists
            if (model.get(GRB_IntAttr_SolCount) == 0) {
                std::cerr << "No feasible solution found within time limit.\n";
                return {routes, bestObj};
            }
 
            bestObj = model.get(GRB_DoubleAttr_ObjBound);
            // double incumbentObj = model.get(GRB_DoubleAttr_ObjVal);
            // double bestBound    = model.get(GRB_DoubleAttr_ObjBound);
            // double mipGapValue  = model.get(GRB_DoubleAttr_MIPGap);
 
            for (int k = 0; k < K_max; ++k) {
                if (u[k].get(GRB_DoubleAttr_X) < 0.5) continue;
 
                std::vector<int> path;
                path.push_back(s);
                int cur = s;
                int safetyCounter = 0;
                while (cur != t && safetyCounter < N) {
                    bool found = false;
                    for (int j = 0; j < N; ++j) {
                        if (j == cur) continue;
                        if (x[k][cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            path.push_back(j);
                            cur = j;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        std::cerr << "Warning: route " << k
                                  << " broken at node " << cur << "\n";
                        break;
                    }
                    ++safetyCounter;
                }
                routes.push_back(std::move(path));
            }
        } else if (status == GRB_INFEASIBLE) {
            std::cerr << "Model is infeasible — cannot cover all nodes "
                         "within the given budget.\n";
            // Optionally compute IIS for debugging:
            // model.computeIIS();
            // model.write("infeasible.ilp");
        } else {
            std::cerr << "Solve failed. Gurobi status: " << status << "\n";
        }
 
    } catch (GRBException& e) {
        std::cerr << "Gurobi error " << e.getErrorCode()
                  << ": " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
 
    return {routes, bestObj};
}


//different starts -> common t
std::pair<std::vector<std::vector<int>>, double>
                    solveMinRoutesCoverDiffR(const ProblemDataMinKDiffR& data,
                                        double mipGap,
                                        double timeLimit,
                                        const std::vector<std::vector<int>>& initialRoutes) {
    const int N      = data.N;
    const int K_max  = data.K;
    const auto& S    = data.startIndices;   // S[k] = start node for route k
    const int t      = data.endIndex;
    const double Budget = data.BudgetPerTeam;
    const auto& C    = data.Cost;

    std::vector<std::vector<int>> routes;
    double bestObj = 0.0;

    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_min_routes.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Min_Routes_Cover");

        if (mipGap >= 0.0)
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);

        // MINIMIZE total number of active routes
        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        // model.set(GRB_DoubleParam_ObjScale, 0.001);

        const int M = std::max(1, N - 1);  // flow capacity upper bound


        //  Variables
        // u[k] ∈ {0,1}: whether route k is active.  obj coeff = 1.
        std::vector<GRBVar> u(K_max);
        for (int k = 0; k < K_max; ++k)
            u[k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
                                "u_" + std::to_string(k));

        // x[k][i][j] ∈ {0,1}: edge (i,j) used by route k
        std::vector<std::vector<std::vector<GRBVar>>> x(
            K_max, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));

        // z[k][i] ∈ {0,1}: node i visited by route k
        std::vector<std::vector<GRBVar>> z(K_max, std::vector<GRBVar>(N));

        // f[k][i][j] ≥ 0: SCF flow on edge (i,j) in route k
        std::vector<std::vector<std::vector<GRBVar>>> f(
            K_max, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));

        for (int k = 0; k < K_max; ++k) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    x[k][i][j] = model.addVar(
                        0.0, 1.0, 0.0, GRB_BINARY,
                        "x_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                    f[k][i][j] = model.addVar(
                        0.0, M, 0.0, GRB_CONTINUOUS,
                        "f_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                }
            }
            for (int i = 0; i < N; ++i)
                z[k][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                    "z_" + std::to_string(k) + "_" + std::to_string(i));
        }

        //  Constraints
        for (int k = 0; k < K_max; ++k) {
            const int s = S[k];   // <-- per-route start node

            // --- (1a) Start node s: out-degree = u[k], in-degree = 0 ---
            GRBLinExpr startOut = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) startOut += x[k][s][j];
            model.addConstr(startOut == u[k],
                "StartOut_k" + std::to_string(k));

            GRBLinExpr startIn = 0.0;
            for (int i = 0; i < N; ++i) if (i != s) startIn += x[k][i][s];
            model.addConstr(startIn == 0.0,
                "StartInZero_k" + std::to_string(k));

            // --- (1b) End node t: in-degree = u[k], out-degree = 0 ---
            GRBLinExpr endIn = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) endIn += x[k][i][t];
            model.addConstr(endIn == u[k],
                "EndIn_k" + std::to_string(k));

            GRBLinExpr endOut = 0.0;
            for (int j = 0; j < N; ++j) if (j != t) endOut += x[k][t][j];
            model.addConstr(endOut == 0.0,
                "EndOutZero_k" + std::to_string(k));

            // --- (1c) Flow conservation at intermediate nodes ---
            for (int v = 0; v < N; ++v) {
                if (v == s || v == t) continue;
                GRBLinExpr out = 0.0, in = 0.0;
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.addConstr(out == in,
                    "DegBal_k" + std::to_string(k) + "_v" + std::to_string(v));
            }

            // --- (1d) Link z to x ---
            for (int v = 0; v < N; ++v) {
                if (v == s || v == t) {
                    // s and t are visited iff route is active
                    model.addConstr(z[k][v] == u[k],
                        "zFixed_k" + std::to_string(k) + "_v" + std::to_string(v));
                } else {
                    GRBLinExpr out = 0.0;
                    for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                    model.addConstr(out == z[k][v],
                        "DefZ_k" + std::to_string(k) + "_v" + std::to_string(v));
                }
            }

            // --- Route k must not visit other routes' start nodes ---
            for (int k2 = 0; k2 < K_max; ++k2) {
                if (k2 == k) continue;
                int s2 = S[k2];
                if (s2 == t) continue;          // degenerate case guard
                if (s2 == s) continue;          // same start node, already handled
                model.addConstr(z[k][s2] == 0,
                    "NoVisitStart_k" + std::to_string(k) + "_s" + std::to_string(s2));
            }
        }

        // --- (2) Cover all customer nodes exactly once ---
        //     Customer = any node that is NOT the end node and NOT a start node.
        std::vector<bool> isStartOrEnd(N, false);
        isStartOrEnd[t] = true;
        for (int k = 0; k < K_max; ++k)
            isStartOrEnd[S[k]] = true;

        for (int v = 0; v < N; ++v) {
            if (isStartOrEnd[v]) continue;
            GRBLinExpr sumz = 0.0;
            for (int k = 0; k < K_max; ++k) sumz += z[k][v];
            model.addConstr(sumz == 1.0, "Cover_v" + std::to_string(v));
        }

        // --- (3) Per-route budget ---
        for (int k = 0; k < K_max; ++k) {
            GRBLinExpr cost = 0.0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    cost += C[i][j] * x[k][i][j];
            model.addConstr(cost <= Budget * u[k],
                "Budget_k" + std::to_string(k));
        }

        // --- (4) SCF subtour elimination ---
        for (int k = 0; k < K_max; ++k) {
            const int s = S[k];

            // Flow capacity
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    model.addConstr(f[k][i][j] <= M * x[k][i][j],
                        "FlowCap_k" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));

            // Flow supply at s: net outflow = number of nodes visited (excl s)
            GRBLinExpr fOutS = 0.0, fInS = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) fOutS += f[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) fInS  += f[k][i][s];

            GRBLinExpr visitsExclS = 0.0;
            for (int v = 0; v < N; ++v) if (v != s) visitsExclS += z[k][v];
            model.addConstr(fOutS - fInS == visitsExclS,
                "FlowSupply_k" + std::to_string(k));

            // Flow demand at every non-s node
            for (int v = 0; v < N; ++v) {
                if (v == s) continue;
                GRBLinExpr fIn = 0.0, fOut = 0.0;
                for (int i = 0; i < N; ++i) if (i != v) fIn  += f[k][i][v];
                for (int j = 0; j < N; ++j) if (j != v) fOut += f[k][v][j];
                model.addConstr(fIn - fOut == z[k][v],
                    "FlowDemand_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
        }

        // --- (5) Symmetry breaking: u[0] ≥ u[1] ≥ … ≥ u[K_max-1] ---
        for (int k = 0; k + 1 < K_max; ++k)
            model.addConstr(u[k] >= u[k + 1],
                "SymBreak_" + std::to_string(k));

        //  MIP warm start from initialRoutes
        if (!initialRoutes.empty()) {
            // Zero out all variables first
            for (int k = 0; k < K_max; ++k) {
                u[k].set(GRB_DoubleAttr_Start, 0.0);
                for (int i = 0; i < N; ++i) {
                    z[k][i].set(GRB_DoubleAttr_Start, 0.0);
                    for (int j = 0; j < N; ++j) if (i != j) {
                        x[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                        f[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                    }
                }
            }

            int nInit = std::min(static_cast<int>(initialRoutes.size()), K_max);
            for (int k = 0; k < nInit; ++k) {
                const auto& route = initialRoutes[k];
                if (route.size() < 2) continue;

                u[k].set(GRB_DoubleAttr_Start, 1.0);

                // Set z for every node in the route
                for (int node : route)
                    z[k][node].set(GRB_DoubleAttr_Start, 1.0);

                // Set x for consecutive edges, and compute SCF flow values.
                const int s = S[k];
                int routeLen = static_cast<int>(route.size());
                for (int idx = 0; idx + 1 < routeLen; ++idx) {
                    int from = route[idx], to = route[idx + 1];
                    x[k][from][to].set(GRB_DoubleAttr_Start, 1.0);

                    int nodesAhead = 0;
                    for (int r = idx + 1; r < routeLen; ++r)
                        if (route[r] != s) ++nodesAhead;
                    f[k][from][to].set(GRB_DoubleAttr_Start,
                                       static_cast<double>(nodesAhead));
                }
            }
        }


        //  Solve
        model.update();
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            if (model.get(GRB_IntAttr_SolCount) == 0) {
                std::cerr << "No feasible solution found within time limit.\n";
                return {routes, bestObj};
            }

            // bestObj = model.get(GRB_DoubleAttr_ObjVal);
            bestObj = model.get(GRB_DoubleAttr_ObjBound);

            for (int k = 0; k < K_max; ++k) {
                if (u[k].get(GRB_DoubleAttr_X) < 0.5) continue;

                const int s = S[k];
                std::vector<int> path;
                path.push_back(s);
                int cur = s;
                int safetyCounter = 0;
                while (cur != t && safetyCounter < N) {
                    bool found = false;
                    for (int j = 0; j < N; ++j) {
                        if (j == cur) continue;
                        if (x[k][cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            path.push_back(j);
                            cur = j;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        std::cerr << "Warning: route " << k
                                  << " broken at node " << cur << "\n";
                        break;
                    }
                    ++safetyCounter;
                }
                routes.push_back(std::move(path));
            }
        } else if (status == GRB_INFEASIBLE) {
            std::cerr << "Model is infeasible — cannot cover all nodes "
                         "within the given budget.\n";
        } else {
            std::cerr << "Solve failed. Gurobi status: " << status << "\n";
        }

    } catch (GRBException& e) {
        std::cerr << "Gurobi error " << e.getErrorCode()
                  << ": " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return {routes, bestObj};
}


//unrooted + open
std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverUnrooted(const ProblemDataMinKUnrooted& data,
                            double mipGap,
                            double timeLimit,
                            const std::vector<std::vector<int>>& initialRoutes) {
    const int N       = data.N;
    const int K_max   = data.K;
    const double B    = data.BudgetPerTeam;
    const auto& SC    = data.SlewCost;
    const auto& DC    = data.DwellCost;

    std::vector<std::vector<int>> routes;
    double bestObj = 0.0;

    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_min_routes_unrooted.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Min_Routes_Cover_Unrooted");

        if (mipGap >= 0.0)
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0.0)
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);

        model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        // model.set(GRB_DoubleParam_ObjScale, 0.001);

        const int D = N;          // virtual depot index
        const int NE = N + 1;     // extended node count
        const int M = N;          // flow capacity upper bound

        // ---- Variables ----

        // u[k] in {0,1}: route k is active.  obj coeff = 1.
        std::vector<GRBVar> u(K_max);
        for (int k = 0; k < K_max; ++k)
            u[k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
                                "u_" + std::to_string(k));

        // x[k][i][j] in {0,1}: directed edge (i,j) used by route k
        //   i,j in 0..NE-1 (includes depot D)
        std::vector<std::vector<std::vector<GRBVar>>> x(
            K_max, std::vector<std::vector<GRBVar>>(NE, std::vector<GRBVar>(NE)));

        // z[k][i] in {0,1}: customer i visited by route k (i in 0..N-1 only)
        std::vector<std::vector<GRBVar>> z(K_max, std::vector<GRBVar>(N));

        // f[k][i][j] >= 0: SCF flow on edge (i,j) in route k
        std::vector<std::vector<std::vector<GRBVar>>> f(
            K_max, std::vector<std::vector<GRBVar>>(NE, std::vector<GRBVar>(NE)));

        for (int k = 0; k < K_max; ++k) {
            for (int i = 0; i < NE; ++i) {
                for (int j = 0; j < NE; ++j) {
                    if (i == j) continue;
                    x[k][i][j] = model.addVar(
                        0.0, 1.0, 0.0, GRB_BINARY,
                        "x_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                    f[k][i][j] = model.addVar(
                        0.0, M, 0.0, GRB_CONTINUOUS,
                        "f_" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                }
            }
            for (int i = 0; i < N; ++i)
                z[k][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                    "z_" + std::to_string(k) + "_" + std::to_string(i));
        }

        // ---- Constraints ----

        for (int k = 0; k < K_max; ++k) {

            // (1a) Depot out-degree = u[k]
            GRBLinExpr depotOut = 0.0;
            for (int j = 0; j < N; ++j) depotOut += x[k][D][j];
            model.addConstr(depotOut == u[k],
                "DepotOut_k" + std::to_string(k));

            // (1b) Depot in-degree = u[k]
            GRBLinExpr depotIn = 0.0;
            for (int i = 0; i < N; ++i) depotIn += x[k][i][D];
            model.addConstr(depotIn == u[k],
                "DepotIn_k" + std::to_string(k));

            // (1c) Flow conservation at customer nodes:
            //      in-degree = out-degree, and both <= 1
            //      A path visits v iff exactly one edge enters and one leaves.
            for (int v = 0; v < N; ++v) {
                GRBLinExpr out = 0.0, in = 0.0;
                for (int j = 0; j < NE; ++j) {
                    if (j == v) continue;
                    out += x[k][v][j];
                    in  += x[k][j][v];
                }
                model.addConstr(out == in,
                    "DegBal_k" + std::to_string(k) + "_v" + std::to_string(v));

                // out-degree <= 1 (simple path, at most one successor)
                model.addConstr(out <= 1,
                    "MaxOut_k" + std::to_string(k) + "_v" + std::to_string(v));
            }

            // (1d) Link z to x: z[k][v] = out-degree of v (0 or 1)
            for (int v = 0; v < N; ++v) {
                GRBLinExpr out = 0.0;
                for (int j = 0; j < NE; ++j) if (j != v) out += x[k][v][j];
                model.addConstr(out == z[k][v],
                    "DefZ_k" + std::to_string(k) + "_v" + std::to_string(v));
            }

            // No edge between depot and itself (already excluded by i==j check)
            // No customer-to-customer edge where i==j (already excluded)
        }

        // (2) Cover every customer exactly once
        for (int v = 0; v < N; ++v) {
            GRBLinExpr sumz = 0.0;
            for (int k = 0; k < K_max; ++k) sumz += z[k][v];
            model.addConstr(sumz == 1.0, "Cover_v" + std::to_string(v));
        }

        // (3) Per-route budget: slew costs on edges + dwell costs on nodes <= B
        for (int k = 0; k < K_max; ++k) {
            GRBLinExpr cost = 0.0;

            // Edge slew costs (customer-to-customer edges only)
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) {
                    if (i == j) continue;
                    cost += SC[i][j] * x[k][i][j];
                }
            // Depot edges have zero slew cost (no physical slew to/from depot)

            // Node dwell costs
            for (int i = 0; i < N; ++i)
                cost += DC[i] * z[k][i];

            model.addConstr(cost <= B * u[k],
                "Budget_k" + std::to_string(k));
        }

        // (4) SCF subtour elimination (source = depot D)
        for (int k = 0; k < K_max; ++k) {
            // Flow capacity
            for (int i = 0; i < NE; ++i)
                for (int j = 0; j < NE; ++j) {
                    if (i == j) continue;
                    model.addConstr(f[k][i][j] <= M * x[k][i][j],
                        "FlowCap_k" + std::to_string(k) + "_" +
                        std::to_string(i) + "_" + std::to_string(j));
                }

            // Flow supply at depot: net outflow = total nodes visited
            GRBLinExpr fOutD = 0.0, fInD = 0.0;
            for (int j = 0; j < N; ++j) fOutD += f[k][D][j];
            for (int i = 0; i < N; ++i) fInD  += f[k][i][D];

            GRBLinExpr totalVisited = 0.0;
            for (int v = 0; v < N; ++v) totalVisited += z[k][v];
            model.addConstr(fOutD - fInD == totalVisited,
                "FlowSupply_k" + std::to_string(k));

            // Flow demand at every customer node: consumes 1 unit if visited
            for (int v = 0; v < N; ++v) {
                GRBLinExpr fIn = 0.0, fOut = 0.0;
                for (int i = 0; i < NE; ++i) {
                    if (i == v) continue;
                    fIn  += f[k][i][v];
                    fOut += f[k][v][i];
                }
                model.addConstr(fIn - fOut == z[k][v],
                    "FlowDemand_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
        }

        // (5) Symmetry breaking (routes are interchangeable since unrooted)
        for (int k = 0; k + 1 < K_max; ++k)
            model.addConstr(u[k] >= u[k + 1],
                "SymBreak_" + std::to_string(k));

        // ---- MIP warm start ----
        if (!initialRoutes.empty()) {
            for (int k = 0; k < K_max; ++k) {
                u[k].set(GRB_DoubleAttr_Start, 0.0);
                for (int i = 0; i < N; ++i) {
                    z[k][i].set(GRB_DoubleAttr_Start, 0.0);
                }
                for (int i = 0; i < NE; ++i)
                    for (int j = 0; j < NE; ++j) {
                        if (i == j) continue;
                        x[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                        f[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
                    }
            }

            int nInit = std::min(static_cast<int>(initialRoutes.size()), K_max);
            for (int k = 0; k < nInit; ++k) {
                const auto& route = initialRoutes[k];
                if (route.empty()) continue;

                u[k].set(GRB_DoubleAttr_Start, 1.0);

                for (int node : route)
                    z[k][node].set(GRB_DoubleAttr_Start, 1.0);

                int routeLen = static_cast<int>(route.size());

                // Depot -> first node
                x[k][D][route[0]].set(GRB_DoubleAttr_Start, 1.0);
                // Customer edges
                for (int idx = 0; idx + 1 < routeLen; ++idx) {
                    x[k][route[idx]][route[idx + 1]].set(GRB_DoubleAttr_Start, 1.0);
                }
                // Last node -> depot
                x[k][route[routeLen - 1]][D].set(GRB_DoubleAttr_Start, 1.0);

                // SCF flow values: depot->first carries routeLen,
                // each subsequent edge carries one less
                f[k][D][route[0]].set(GRB_DoubleAttr_Start,
                                      static_cast<double>(routeLen));
                for (int idx = 0; idx + 1 < routeLen; ++idx) {
                    f[k][route[idx]][route[idx + 1]].set(
                        GRB_DoubleAttr_Start,
                        static_cast<double>(routeLen - 1 - idx));
                }
                // last->depot carries 0 flow (already zeroed)
            }
        }

        // ---- Solve ----
        model.update();
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            if (model.get(GRB_IntAttr_SolCount) == 0) {
                std::cerr << "No feasible solution found within time limit.\n";
                return {routes, bestObj};
            }

            // bestObj = model.get(GRB_DoubleAttr_ObjVal);
            bestObj = model.get(GRB_DoubleAttr_ObjBound);

            for (int k = 0; k < K_max; ++k) {
                if (u[k].get(GRB_DoubleAttr_X) < 0.5) continue;

                // Find the first node: depot -> first
                int first = -1;
                for (int j = 0; j < N; ++j) {
                    if (x[k][D][j].get(GRB_DoubleAttr_X) > 0.5) {
                        first = j;
                        break;
                    }
                }
                if (first < 0) {
                    std::cerr << "Warning: route " << k
                              << " active but no depot-out edge found.\n";
                    continue;
                }

                // Trace the path: first -> ... -> last (stop when we reach depot)
                std::vector<int> path;
                path.push_back(first);
                int cur = first;
                int safetyCounter = 0;
                while (safetyCounter < N) {
                    bool found = false;
                    for (int j = 0; j < NE; ++j) {
                        if (j == cur) continue;
                        if (x[k][cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            if (j == D) {
                                found = true;  // reached depot, path ends
                                break;
                            }
                            path.push_back(j);
                            cur = j;
                            found = true;
                            break;
                        }
                    }
                    if (!found || cur == first) break;  // back to depot or stuck
                    ++safetyCounter;
                    // Check if next hop is depot
                    if (x[k][cur][D].get(GRB_DoubleAttr_X) > 0.5) break;
                }
                routes.push_back(std::move(path));
            }
        } else if (status == GRB_INFEASIBLE) {
            std::cerr << "Model is infeasible — possible causes:\n"
                         "  - per-route budget too tight\n"
                         "  - not enough routes (K_max="
                      << K_max << ") to cover all " << N << " nodes\n";
        } else {
            std::cerr << "Solve failed. Gurobi status: " << status << "\n";
        }

    } catch (GRBException& e) {
        std::cerr << "Gurobi error " << e.getErrorCode()
                  << ": " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return {routes, bestObj};
}


// //unrooted + open
// std::pair<std::vector<std::vector<int>>, double>
// solveMinRoutesCoverRootedOpen(const ProblemDataMinKRootedOpen& data,
//                               double mipGap,
//                               double timeLimit,
//                               const std::vector<std::vector<int>>& initialRoutes) {
//     const int N       = data.N;
//     const int K_max   = data.K;
//     const auto& S     = data.startIndices;
//     const double B    = data.BudgetPerTeam;
//     const auto& SC    = data.SlewCost;
//     const auto& DC    = data.DwellCost;

//     // Extended graph: nodes 0..N-1 are real, node D=N is virtual depot
//     const int D  = N;
//     const int NE = N + 1;
//     const int M  = N;          // SCF flow capacity upper bound

//     std::vector<std::vector<int>> routes;
//     double bestObj = 0.0;

//     try {
//         GRBEnv env(true);
//         env.set("LogFile", "gurobi_min_routes_rooted_open.log");
//         env.start();

//         GRBModel model(env);
//         model.set(GRB_StringAttr_ModelName, "Min_Routes_Cover_Rooted_Open");

//         if (mipGap >= 0.0)
//             model.set(GRB_DoubleParam_MIPGap, mipGap);
//         if (timeLimit > 0.0)
//             model.set(GRB_DoubleParam_TimeLimit, timeLimit);

//         model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
//         model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
//         model.set(GRB_DoubleParam_ObjScale, 0.001);
        
//         // ---- Variables ----

//         // u[k] in {0,1}: route k is active
//         std::vector<GRBVar> u(K_max);
//         for (int k = 0; k < K_max; ++k)
//             u[k] = model.addVar(0.0, 1.0, 1.0, GRB_BINARY,
//                                 "u_" + std::to_string(k));

//         // x[k][i][j] in {0,1}: directed edge (i,j) used by route k
//         std::vector<std::vector<std::vector<GRBVar>>> x(
//             K_max, std::vector<std::vector<GRBVar>>(NE, std::vector<GRBVar>(NE)));

//         // z[k][i] in {0,1}: node i visited by route k (i in 0..N-1)
//         std::vector<std::vector<GRBVar>> z(K_max, std::vector<GRBVar>(N));

//         // f[k][i][j] >= 0: SCF flow
//         std::vector<std::vector<std::vector<GRBVar>>> f(
//             K_max, std::vector<std::vector<GRBVar>>(NE, std::vector<GRBVar>(NE)));

//         for (int k = 0; k < K_max; ++k) {
//             for (int i = 0; i < NE; ++i) {
//                 for (int j = 0; j < NE; ++j) {
//                     if (i == j) continue;
//                     x[k][i][j] = model.addVar(
//                         0.0, 1.0, 0.0, GRB_BINARY,
//                         "x_" + std::to_string(k) + "_" +
//                         std::to_string(i) + "_" + std::to_string(j));
//                     f[k][i][j] = model.addVar(
//                         0.0, M, 0.0, GRB_CONTINUOUS,
//                         "f_" + std::to_string(k) + "_" +
//                         std::to_string(i) + "_" + std::to_string(j));
//                 }
//             }
//             for (int i = 0; i < N; ++i)
//                 z[k][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
//                     "z_" + std::to_string(k) + "_" + std::to_string(i));
//         }

//         // ---- Identify start and non-customer nodes ----
//         std::vector<bool> isStart(N, false);
//         for (int k = 0; k < K_max; ++k)
//             if (S[k] < N) isStart[S[k]] = true;

//         // ---- Constraints ----

//         for (int k = 0; k < K_max; ++k) {
//             const int s = S[k];

//             // (1a) Start node s_k: out-degree = u[k], in-degree = 0
//             //      (no edge from depot to s; s is the real origin)
//             GRBLinExpr startOut = 0.0;
//             for (int j = 0; j < NE; ++j) {
//                 if (j == s) continue;
//                 startOut += x[k][s][j];
//             }
//             model.addConstr(startOut == u[k],
//                 "StartOut_k" + std::to_string(k));

//             GRBLinExpr startIn = 0.0;
//             for (int i = 0; i < NE; ++i) {
//                 if (i == s) continue;
//                 startIn += x[k][i][s];
//             }
//             model.addConstr(startIn == 0,
//                 "StartInZero_k" + std::to_string(k));

//             // (1b) Virtual depot D: in-degree = u[k], out-degree = 0
//             //      (path ends by sending last customer -> D)
//             GRBLinExpr depotIn = 0.0;
//             for (int i = 0; i < N; ++i) depotIn += x[k][i][D];
//             model.addConstr(depotIn == u[k],
//                 "DepotIn_k" + std::to_string(k));

//             GRBLinExpr depotOut = 0.0;
//             for (int j = 0; j < N; ++j) depotOut += x[k][D][j];
//             model.addConstr(depotOut == 0,
//                 "DepotOutZero_k" + std::to_string(k));

//             // (1c) Flow conservation at customer nodes (excluding s_k)
//             for (int v = 0; v < N; ++v) {
//                 if (v == s) continue;
//                 GRBLinExpr out = 0.0, in = 0.0;
//                 for (int j = 0; j < NE; ++j) {
//                     if (j == v) continue;
//                     out += x[k][v][j];
//                     in  += x[k][j][v];
//                 }
//                 model.addConstr(out == in,
//                     "DegBal_k" + std::to_string(k) + "_v" + std::to_string(v));

//                 // At most one successor (simple path)
//                 model.addConstr(out <= 1,
//                     "MaxOut_k" + std::to_string(k) + "_v" + std::to_string(v));
//             }

//             // (1d) Link z to x
//             for (int v = 0; v < N; ++v) {
//                 if (v == s) {
//                     // Start node is visited iff route is active
//                     model.addConstr(z[k][v] == u[k],
//                         "zFixed_k" + std::to_string(k) + "_v" + std::to_string(v));
//                 } else {
//                     // For non-start customers, z = in-degree from real nodes
//                     // (equivalent to out-degree since flow is conserved)
//                     GRBLinExpr out = 0.0;
//                     for (int j = 0; j < NE; ++j) if (j != v) out += x[k][v][j];
//                     model.addConstr(out == z[k][v],
//                         "DefZ_k" + std::to_string(k) + "_v" + std::to_string(v));
//                 }
//             }

//             // (1e) Route k must not visit other routes' start nodes
//             for (int k2 = 0; k2 < K_max; ++k2) {
//                 if (k2 == k) continue;
//                 int s2 = S[k2];
//                 if (s2 == s) continue;
//                 model.addConstr(z[k][s2] == 0,
//                     "NoVisitStart_k" + std::to_string(k) + "_s" + std::to_string(s2));
//             }
//         }

//         // (2) Cover every customer node exactly once
//         //     Start nodes are excluded (they are depots, not targets)
//         for (int v = 0; v < N; ++v) {
//             if (isStart[v]) continue;
//             GRBLinExpr sumz = 0.0;
//             for (int k = 0; k < K_max; ++k) sumz += z[k][v];
//             model.addConstr(sumz == 1.0, "Cover_v" + std::to_string(v));
//         }

//         // (3) Per-route budget: slew costs + dwell costs <= B
//         //     Only customer-to-customer edges carry slew cost.
//         //     Only non-start visited nodes carry dwell cost.
//         //     Edges to/from depot D have zero cost.
//         for (int k = 0; k < K_max; ++k) {
//             const int s = S[k];
//             GRBLinExpr cost = 0.0;

//             // Edge slew costs (customer-to-customer only)
//             for (int i = 0; i < N; ++i)
//                 for (int j = 0; j < N; ++j) {
//                     if (i == j) continue;
//                     cost += SC[i][j] * x[k][i][j];
//                 }

//             // Node dwell costs (skip start node — it's a depot)
//             for (int i = 0; i < N; ++i) {
//                 if (i == s) continue;
//                 cost += DC[i] * z[k][i];
//             }

//             model.addConstr(cost <= B * u[k],
//                 "Budget_k" + std::to_string(k));
//         }

//         // (4) SCF subtour elimination (source = start node s_k)
//         for (int k = 0; k < K_max; ++k) {
//             const int s = S[k];

//             // Flow capacity
//             for (int i = 0; i < NE; ++i)
//                 for (int j = 0; j < NE; ++j) {
//                     if (i == j) continue;
//                     model.addConstr(f[k][i][j] <= M * x[k][i][j],
//                         "FlowCap_k" + std::to_string(k) + "_" +
//                         std::to_string(i) + "_" + std::to_string(j));
//                 }

//             // Flow supply at s: net outflow = nodes visited excluding s
//             GRBLinExpr fOutS = 0.0, fInS = 0.0;
//             for (int j = 0; j < NE; ++j) {
//                 if (j == s) continue;
//                 fOutS += f[k][s][j];
//                 fInS  += f[k][j][s];  // should be zero (no in-edges to s)
//             }

//             GRBLinExpr visitsExclS = 0.0;
//             for (int v = 0; v < N; ++v) if (v != s) visitsExclS += z[k][v];
//             // +1 for the depot (depot is always "visited" when route is active)
//             visitsExclS += u[k];

//             model.addConstr(fOutS - fInS == visitsExclS,
//                 "FlowSupply_k" + std::to_string(k));

//             // Flow demand at every non-s node (customers + depot)
//             for (int v = 0; v < NE; ++v) {
//                 if (v == s) continue;
//                 GRBLinExpr fIn = 0.0, fOut = 0.0;
//                 for (int i = 0; i < NE; ++i) {
//                     if (i == v) continue;
//                     fIn  += f[k][i][v];
//                     fOut += f[k][v][i];
//                 }
//                 if (v == D) {
//                     // Depot consumes 1 unit of flow when route is active
//                     model.addConstr(fIn - fOut == u[k],
//                         "FlowDemand_k" + std::to_string(k) + "_depot");
//                 } else {
//                     model.addConstr(fIn - fOut == z[k][v],
//                         "FlowDemand_k" + std::to_string(k) + "_v" + std::to_string(v));
//                 }
//             }
//         }


//         if (!initialRoutes.empty()) {
//             // Zero out all variables first
//             for (int k = 0; k < K_max; ++k) {
//                 u[k].set(GRB_DoubleAttr_Start, 0.0);
//                 for (int i = 0; i < N; ++i)
//                     z[k][i].set(GRB_DoubleAttr_Start, 0.0);
//                 for (int i = 0; i < NE; ++i)
//                     for (int j = 0; j < NE; ++j) {
//                         if (i == j) continue;
//                         x[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
//                         f[k][i][j].set(GRB_DoubleAttr_Start, 0.0);
//                     }
//             }

//             // Build a map from start node -> ILP route index
//             std::unordered_map<int, int> startToK;
//             for (int k = 0; k < K_max; ++k)
//                 startToK[S[k]] = k;

//             for (const auto& route : initialRoutes) {
//                 if (route.empty()) continue;

//                 int routeStart = route.front();
//                 auto it = startToK.find(routeStart);
//                 if (it == startToK.end()) continue;  // no matching ILP route

//                 int k = it->second;
//                 u[k].set(GRB_DoubleAttr_Start, 1.0);

//                 for (int node : route)
//                     z[k][node].set(GRB_DoubleAttr_Start, 1.0);

//                 int routeLen = static_cast<int>(route.size());

//                 // Edges along the path
//                 for (int idx = 0; idx + 1 < routeLen; ++idx)
//                     x[k][route[idx]][route[idx + 1]].set(GRB_DoubleAttr_Start, 1.0);

//                 // Last node -> depot
//                 x[k][route[routeLen - 1]][D].set(GRB_DoubleAttr_Start, 1.0);

//                 // SCF flow values
//                 int totalFlow = routeLen;  // non-s customers + depot
//                 for (int idx = 0; idx + 1 < routeLen; ++idx) {
//                     f[k][route[idx]][route[idx + 1]].set(
//                         GRB_DoubleAttr_Start,
//                         static_cast<double>(totalFlow - idx - 1));
//                 }
//                 f[k][route[routeLen - 1]][D].set(GRB_DoubleAttr_Start, 1.0);
//             }
//         }

//         // ---- Solve ----
//         model.update();
//         model.optimize();

//         int status = model.get(GRB_IntAttr_Status);
//         if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
//             if (model.get(GRB_IntAttr_SolCount) == 0) {
//                 std::cerr << "No feasible solution found within time limit.\n";
//                 return {routes, bestObj};
//             }

//             // bestObj = model.get(GRB_DoubleAttr_ObjVal);
//             bestObj = model.get(GRB_DoubleAttr_ObjBound);

//             for (int k = 0; k < K_max; ++k) {
//                 if (u[k].get(GRB_DoubleAttr_X) < 0.5) continue;

//                 const int s = S[k];
//                 std::vector<int> path;
//                 path.push_back(s);
//                 int cur = s;
//                 int safetyCounter = 0;
//                 while (safetyCounter < N) {
//                     bool found = false;
//                     for (int j = 0; j < NE; ++j) {
//                         if (j == cur) continue;
//                         if (x[k][cur][j].get(GRB_DoubleAttr_X) > 0.5) {
//                             if (j == D) break;  // reached virtual depot, done
//                             path.push_back(j);
//                             cur = j;
//                             found = true;
//                             break;
//                         }
//                     }
//                     if (!found) break;
//                     ++safetyCounter;
//                 }
//                 routes.push_back(std::move(path));
//             }
//         } else if (status == GRB_INFEASIBLE) {
//             std::cerr << "Model is infeasible — possible causes:\n"
//                          "  - per-route budget too tight\n"
//                          "  - not enough routes (K_max="
//                       << K_max << ") to cover all " << N << " nodes\n";
//         } else {
//             std::cerr << "Solve failed. Gurobi status: " << status << "\n";
//         }

//     } catch (GRBException& e) {
//         std::cerr << "Gurobi error " << e.getErrorCode()
//                   << ": " << e.getMessage() << "\n";
//     } catch (std::exception& e) {
//         std::cerr << "Error: " << e.what() << "\n";
//     }

//     return {routes, bestObj};
// }


std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverRootedOpen(const ProblemDataMinKRootedOpen& data,
                              double mipGap,
                              double timeLimit,
                              const std::vector<std::vector<int>>& initialRoutes) {
    const int N       = data.N;
    const int K_max   = data.K;
    const auto& S     = data.startIndices;
    const double B    = data.BudgetPerTeam;
    const auto& SC    = data.SlewCost;
    const auto& DC    = data.DwellCost;

    constexpr double EPS = 1e-7;
    const double INF = std::numeric_limits<double>::infinity();

    std::vector<std::vector<int>> routes;
    double bestObjBound = std::numeric_limits<double>::quiet_NaN();

    std::vector<std::vector<int>> fallbackRoutes;
    bool warmStartValid = false;

    try {
        // ------------------------------------------------------------------
        // Basic validation
        // ------------------------------------------------------------------
        if (N <= 0 || K_max <= 0) {
            std::cerr << "[compact-flow ILP] Invalid problem size: N="
                      << N << ", K=" << K_max << "\n";
            return {routes, bestObjBound};
        }

        if ((int)S.size() != K_max) {
            throw std::runtime_error(
                "data.K does not match data.startIndices.size().");
        }

        if ((int)SC.size() != N) {
            throw std::runtime_error("SlewCost has wrong number of rows.");
        }

        for (int i = 0; i < N; ++i) {
            if ((int)SC[i].size() != N) {
                throw std::runtime_error("SlewCost is not an N x N matrix.");
            }
        }

        if ((int)DC.size() != N) {
            throw std::runtime_error("DwellCost has wrong size.");
        }

        if (!(B > 0.0) || !std::isfinite(B)) {
            throw std::runtime_error("Budget must be positive and finite.");
        }

        for (int k = 0; k < K_max; ++k) {
            if (S[k] < 0 || S[k] >= N) {
                throw std::runtime_error("Invalid start index in startIndices.");
            }
        }

        for (int i = 0; i < N; ++i) {
            if (std::isfinite(DC[i]) && DC[i] < -EPS) {
                throw std::runtime_error("Negative dwell cost detected.");
            }
        }

        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                if (std::isfinite(SC[i][j]) && SC[i][j] < -EPS) {
                    throw std::runtime_error("Negative slew cost detected.");
                }
            }
        }

        // ------------------------------------------------------------------
        // Customer set = all nodes except physical start nodes.
        // Duplicate starts are allowed as separate route slots, but the
        // physical start node is still not a customer.
        // ------------------------------------------------------------------
        std::vector<char> isStart(N, 0);
        for (int k = 0; k < K_max; ++k) {
            isStart[S[k]] = 1;
        }

        std::vector<int> customers;
        std::vector<int> custPos(N, -1);

        for (int v = 0; v < N; ++v) {
            if (!isStart[v]) {
                if (!std::isfinite(DC[v])) {
                    throw std::runtime_error(
                        "A customer node has non-finite dwell cost.");
                }

                custPos[v] = static_cast<int>(customers.size());
                customers.push_back(v);
            }
        }

        const int C = static_cast<int>(customers.size());

        if (C == 0) {
            std::cout << "[compact-flow ILP] No customers to cover.\n";
            bestObjBound = 0.0;
            return {routes, bestObjBound};
        }

        // ------------------------------------------------------------------
        // Helper: validate route collection in greedy/output format.
        // A route is [start, customer, customer, ...].
        // ------------------------------------------------------------------
        auto validateRoutes = [&](const std::vector<std::vector<int>>& candRoutes,
                                  const std::string& label) -> bool {
            std::unordered_map<int, int> availableStartSlots;
            for (int k = 0; k < K_max; ++k) {
                availableStartSlots[S[k]] += 1;
            }

            std::unordered_map<int, int> usedStartSlots;
            std::vector<int> cover(C, 0);

            bool ok = true;

            for (int r = 0; r < (int)candRoutes.size(); ++r) {
                const auto& path = candRoutes[r];

                if (path.size() < 2) {
                    std::cerr << "[" << label << "] route " << r
                              << " has fewer than 2 nodes.\n";
                    ok = false;
                    continue;
                }

                const int s = path.front();

                if (!availableStartSlots.count(s)) {
                    std::cerr << "[" << label << "] route " << r
                              << " starts at non-start node " << s << ".\n";
                    ok = false;
                    continue;
                }

                usedStartSlots[s] += 1;
                if (usedStartSlots[s] > availableStartSlots[s]) {
                    std::cerr << "[" << label << "] too many routes use start "
                              << s << ".\n";
                    ok = false;
                }

                std::vector<char> seenInRoute(C, 0);
                double cost = 0.0;

                for (int p = 1; p < (int)path.size(); ++p) {
                    const int v = path[p];

                    if (v < 0 || v >= N || custPos[v] < 0) {
                        std::cerr << "[" << label << "] route " << r
                                  << " contains invalid customer node "
                                  << v << ".\n";
                        ok = false;
                        continue;
                    }

                    const int c = custPos[v];

                    if (seenInRoute[c]) {
                        std::cerr << "[" << label << "] route " << r
                                  << " repeats customer " << v << ".\n";
                        ok = false;
                    }

                    seenInRoute[c] = 1;
                    cover[c] += 1;

                    const int prev = path[p - 1];

                    if (prev < 0 || prev >= N) {
                        std::cerr << "[" << label << "] route " << r
                                  << " has invalid predecessor node "
                                  << prev << ".\n";
                        ok = false;
                        continue;
                    }

                    if (!std::isfinite(SC[prev][v])) {
                        std::cerr << "[" << label << "] route " << r
                                  << " uses non-finite slew edge "
                                  << prev << " -> " << v << ".\n";
                        ok = false;
                        continue;
                    }

                    cost += SC[prev][v] + DC[v];
                }

                if (cost > B + 1e-6) {
                    std::cerr << "[" << label << "] route " << r
                              << " exceeds budget: cost=" << cost
                              << ", B=" << B << ".\n";
                    ok = false;
                }
            }

            int missing = 0;
            int duplicate = 0;

            for (int c = 0; c < C; ++c) {
                if (cover[c] == 0) ++missing;
                if (cover[c] > 1) ++duplicate;
            }

            if (missing > 0 || duplicate > 0) {
                std::cerr << "[" << label << "] missing customers=" << missing
                          << ", duplicate customers=" << duplicate << ".\n";
                ok = false;
            }

            return ok;
        };

        if (!initialRoutes.empty()) {
            warmStartValid = validateRoutes(initialRoutes, "initialRoutes");
            if (warmStartValid) {
                fallbackRoutes = initialRoutes;
            }
        }

        // ------------------------------------------------------------------
        // Compute an optimistic lower bound to reach each customer.
        // This is used only for safe pruning of impossible customer arcs.
        // ------------------------------------------------------------------
        std::vector<double> reachLB(C, INF);

        for (int b = 0; b < C; ++b) {
            const int j = customers[b];

            for (int k = 0; k < K_max; ++k) {
                const int s = S[k];

                if (!std::isfinite(SC[s][j])) continue;

                const double c = SC[s][j] + DC[j];

                if (std::isfinite(c) && c < reachLB[b]) {
                    reachLB[b] = c;
                }
            }
        }

        // Dense Dijkstra over customers. Costs are nonnegative by validation.
        {
            std::vector<char> done(C, 0);

            for (int iter = 0; iter < C; ++iter) {
                int a = -1;
                double best = INF;

                for (int i = 0; i < C; ++i) {
                    if (!done[i] && reachLB[i] < best) {
                        best = reachLB[i];
                        a = i;
                    }
                }

                if (a < 0 || !std::isfinite(best)) break;

                done[a] = 1;

                const int u = customers[a];

                for (int b = 0; b < C; ++b) {
                    if (done[b] || b == a) continue;

                    const int v = customers[b];

                    if (!std::isfinite(SC[u][v])) continue;

                    const double step = SC[u][v] + DC[v];

                    if (!std::isfinite(step)) continue;

                    if (reachLB[a] + step < reachLB[b]) {
                        reachLB[b] = reachLB[a] + step;
                    }
                }
            }
        }

        // ------------------------------------------------------------------
        // Gurobi model
        // ------------------------------------------------------------------
        GRBEnv env(true);
        env.set("LogFile", "gurobi_min_routes_rooted_open_compact_flow.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName,
                  "Min_Routes_Cover_Rooted_Open_Compact_Flow");

        if (mipGap >= 0.0) {
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        }

        if (timeLimit > 0.0) {
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        }

        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);

        // Memory-safer defaults. Tune these for your machine.
        model.set(GRB_IntParam_Threads, 8);
        model.set(GRB_DoubleParam_NodefileStart, 1.0);
        // Let Gurobi choose the root algorithm. Set Method=1 if dual simplex
        // is better on your instances.
        model.set(GRB_IntParam_Method, -1);

        // ------------------------------------------------------------------
        // Arc structures
        // ------------------------------------------------------------------
        struct StartArc {
            int k;       // route/start slot
            int to;      // customer position
            int sNode;   // original start node
            int jNode;   // original customer node
            double cost; // SC[s][j] + DC[j]
            GRBVar x;
        };

        struct CustArc {
            int from;    // customer position
            int to;      // customer position
            int iNode;   // original customer node
            int jNode;   // original customer node
            double cost; // SC[i][j] + DC[j]
            GRBVar x;
        };

        std::vector<StartArc> startArcs;
        std::vector<CustArc> custArcs;

        std::vector<std::vector<int>> startOut(K_max);
        std::vector<std::vector<int>> inFromStart(C);

        std::vector<std::vector<int>> outCust(C);
        std::vector<std::vector<int>> inCust(C);

        std::vector<std::vector<int>> startArcId(
            K_max, std::vector<int>(C, -1));

        std::vector<std::vector<int>> custArcId(
            C, std::vector<int>(C, -1));

        // ------------------------------------------------------------------
        // Main variables
        // ------------------------------------------------------------------

        // y[k] = 1 iff start/team slot k is used.
        std::vector<GRBVar> y(K_max);

        for (int k = 0; k < K_max; ++k) {
            y[k] = model.addVar(
                0.0, 1.0, 0.0, GRB_BINARY,
                "y_" + std::to_string(k)
            );
        }

        // t[c] = cumulative route cost when reaching customer c.
        // ord[c] = path order, used to cut subtours.
        // endVar[c] = 1 iff customer c is terminal.
        std::vector<GRBVar> t(C);
        std::vector<GRBVar> ord(C);
        std::vector<GRBVar> endVar(C);

        for (int c = 0; c < C; ++c) {
            const int node = customers[c];

            t[c] = model.addVar(
                0.0, B, 0.0, GRB_CONTINUOUS,
                "t_" + std::to_string(node)
            );

            ord[c] = model.addVar(
                1.0, static_cast<double>(C), 0.0, GRB_CONTINUOUS,
                "ord_" + std::to_string(node)
            );

            endVar[c] = model.addVar(
                0.0, 1.0, 0.0, GRB_BINARY,
                "end_" + std::to_string(node)
            );
        }

        // Start -> customer arcs.
        for (int k = 0; k < K_max; ++k) {
            const int s = S[k];

            for (int b = 0; b < C; ++b) {
                const int j = customers[b];

                if (!std::isfinite(SC[s][j])) continue;

                const double c = SC[s][j] + DC[j];

                // A one-customer route must fit.
                if (!std::isfinite(c) || c > B + EPS) continue;

                StartArc a;
                a.k = k;
                a.to = b;
                a.sNode = s;
                a.jNode = j;
                a.cost = c;
                a.x = model.addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "xs_" + std::to_string(k) + "_" + std::to_string(j)
                );

                const int id = static_cast<int>(startArcs.size());
                startArcs.push_back(a);
                startArcId[k][b] = id;
                startOut[k].push_back(id);
                inFromStart[b].push_back(id);
            }
        }

        // Customer -> customer arcs.
        for (int a = 0; a < C; ++a) {
            const int i = customers[a];

            // If no finite way to reach i even in the optimistic relaxation,
            // outgoing arcs from i cannot be used in any rooted path.
            if (!std::isfinite(reachLB[a])) continue;

            for (int b = 0; b < C; ++b) {
                if (a == b) continue;

                const int j = customers[b];

                if (!std::isfinite(SC[i][j])) continue;

                const double step = SC[i][j] + DC[j];

                if (!std::isfinite(step)) continue;

                // If one step alone exceeds budget, impossible because t[i] >= 0.
                if (step > B + EPS) continue;

                // Safe pruning: every route reaching i costs at least reachLB[i].
                if (reachLB[a] + step > B + EPS) continue;

                CustArc e;
                e.from = a;
                e.to = b;
                e.iNode = i;
                e.jNode = j;
                e.cost = step;
                e.x = model.addVar(
                    0.0, 1.0, 0.0, GRB_BINARY,
                    "xc_" + std::to_string(i) + "_" + std::to_string(j)
                );

                const int id = static_cast<int>(custArcs.size());
                custArcs.push_back(e);
                custArcId[a][b] = id;
                outCust[a].push_back(id);
                inCust[b].push_back(id);
            }
        }

        std::cout << "[compact-flow ILP] N=" << N
                  << ", customers=" << C
                  << ", start slots=" << K_max
                  << ", startArcs=" << startArcs.size()
                  << ", custArcs=" << custArcs.size()
                  << "\n";

        // ------------------------------------------------------------------
        // Objective
        // ------------------------------------------------------------------
        GRBLinExpr obj = 0.0;
        GRBLinExpr sumY = 0.0;

        for (int k = 0; k < K_max; ++k) {
            obj += y[k];
            sumY += y[k];
        }

        model.setObjective(obj, GRB_MINIMIZE);

        // ------------------------------------------------------------------
        // Degree/path constraints
        // ------------------------------------------------------------------

        // Each customer has exactly one predecessor:
        // either a start arc or a customer arc.
        for (int b = 0; b < C; ++b) {
            GRBLinExpr pred = 0.0;

            for (int id : inFromStart[b]) {
                pred += startArcs[id].x;
            }

            for (int id : inCust[b]) {
                pred += custArcs[id].x;
            }

            model.addConstr(
                pred == 1.0,
                "pred_" + std::to_string(customers[b])
            );
        }

        // Each customer has exactly one successor, or it is terminal.
        for (int a = 0; a < C; ++a) {
            GRBLinExpr succ = endVar[a];

            for (int id : outCust[a]) {
                succ += custArcs[id].x;
            }

            model.addConstr(
                succ == 1.0,
                "succ_" + std::to_string(customers[a])
            );
        }

        // Each used start launches exactly one path.
        for (int k = 0; k < K_max; ++k) {
            GRBLinExpr out = 0.0;

            for (int id : startOut[k]) {
                out += startArcs[id].x;
            }

            model.addConstr(
                out == y[k],
                "start_out_" + std::to_string(k)
            );
        }

        // Since C > 0, at least one route is needed. This is redundant with
        // the unit-flow constraints, but it helps the root relaxation.
        model.addConstr(sumY >= 1.0, "at_least_one_route");

        // ------------------------------------------------------------------
        // Resource/order constraints.
        // These are kept from the compact model. The new compact flows below
        // are the main strengthening, but these constraints still help prune.
        // ------------------------------------------------------------------

        for (const auto& a : startArcs) {
            // If x=1: t[j] >= SC[start][j] + DC[j].
            model.addConstr(
                t[a.to] >= a.cost - B * (1.0 - a.x),
                "res_start_" + std::to_string(a.k) + "_" +
                std::to_string(a.jNode)
            );

            model.addConstr(
                ord[a.to] >= 1.0 - static_cast<double>(C) * (1.0 - a.x),
                "ord_start_" + std::to_string(a.k) + "_" +
                std::to_string(a.jNode)
            );
        }

        for (const auto& e : custArcs) {
            // If x=1:
            //   t[j] >= t[i] + SC[i][j] + DC[j]
            // If x=0, this is relaxed for all t[i], t[j] in [0,B].
            const double Mres = B + e.cost;

            model.addConstr(
                t[e.to] >= t[e.from] + e.cost - Mres * (1.0 - e.x),
                "res_cust_" + std::to_string(e.iNode) + "_" +
                std::to_string(e.jNode)
            );

            model.addConstr(
                ord[e.to] >= ord[e.from] + 1.0
                             - static_cast<double>(C) * (1.0 - e.x),
                "ord_cust_" + std::to_string(e.iNode) + "_" +
                std::to_string(e.jNode)
            );
        }

        // ------------------------------------------------------------------
        // Compact strengthening flows
        // ------------------------------------------------------------------
        // h-flow: unit connectivity flow. Each customer consumes 1 unit.
        //         This kills disconnected customer-only fractional cycles.
        //
        // g-flow: downstream route-cost flow. On a path
        //         s -> a -> b -> c,
        //             g[s,a] = cost(s,a)+cost(a,b)+cost(b,c),
        //             g[a,b] = cost(a,b)+cost(b,c),
        //             g[b,c] = cost(b,c).
        //         Therefore every start arc carries the total cost of its
        //         route and is capped by B*x. This gives a useful route-count
        //         lower bound without route-indexing all customer arcs.
        // ------------------------------------------------------------------

        std::vector<GRBVar> hs(startArcs.size());
        std::vector<GRBVar> hc(custArcs.size());
        std::vector<GRBVar> gs(startArcs.size());
        std::vector<GRBVar> gc(custArcs.size());

        const double HCAP = static_cast<double>(C);

        for (int id = 0; id < (int)startArcs.size(); ++id) {
            const auto& a = startArcs[id];

            hs[id] = model.addVar(
                0.0, HCAP, 0.0, GRB_CONTINUOUS,
                "h_start_" + std::to_string(a.k) + "_" +
                std::to_string(a.jNode)
            );

            gs[id] = model.addVar(
                0.0, B, 0.0, GRB_CONTINUOUS,
                "g_start_" + std::to_string(a.k) + "_" +
                std::to_string(a.jNode)
            );

            model.addConstr(
                hs[id] <= HCAP * a.x,
                "hcap_start_" + std::to_string(id)
            );

            // If an arc is selected integrally, it carries at least the
            // customer reached by that arc.
            model.addConstr(
                hs[id] >= a.x,
                "hlb_start_" + std::to_string(id)
            );

            model.addConstr(
                gs[id] <= B * a.x,
                "gcap_start_" + std::to_string(id)
            );

            model.addConstr(
                gs[id] >= a.cost * a.x,
                "glb_start_" + std::to_string(id)
            );
        }

        for (int id = 0; id < (int)custArcs.size(); ++id) {
            const auto& e = custArcs[id];

            hc[id] = model.addVar(
                0.0, HCAP, 0.0, GRB_CONTINUOUS,
                "h_cust_" + std::to_string(e.iNode) + "_" +
                std::to_string(e.jNode)
            );

            gc[id] = model.addVar(
                0.0, B, 0.0, GRB_CONTINUOUS,
                "g_cust_" + std::to_string(e.iNode) + "_" +
                std::to_string(e.jNode)
            );

            model.addConstr(
                hc[id] <= HCAP * e.x,
                "hcap_cust_" + std::to_string(id)
            );

            model.addConstr(
                hc[id] >= e.x,
                "hlb_cust_" + std::to_string(id)
            );

            model.addConstr(
                gc[id] <= B * e.x,
                "gcap_cust_" + std::to_string(id)
            );

            model.addConstr(
                gc[id] >= e.cost * e.x,
                "glb_cust_" + std::to_string(id)
            );
        }

        GRBLinExpr totalSelectedArcCost = 0.0;

        for (const auto& a : startArcs) {
            totalSelectedArcCost += a.cost * a.x;
        }

        for (const auto& e : custArcs) {
            totalSelectedArcCost += e.cost * e.x;
        }

        // Very cheap valid inequality. The g-flow constraints imply a similar
        // relation, but this explicit cut often helps root presolve.
        model.addConstr(
            totalSelectedArcCost <= B * sumY,
            "global_budget_cut"
        );

        // Flow balances at every customer.
        for (int b = 0; b < C; ++b) {
            GRBLinExpr hIn = 0.0;
            GRBLinExpr hOut = 0.0;

            GRBLinExpr gIn = 0.0;
            GRBLinExpr gOut = 0.0;
            GRBLinExpr incomingArcCost = 0.0;

            for (int id : inFromStart[b]) {
                hIn += hs[id];
                gIn += gs[id];
                incomingArcCost += startArcs[id].cost * startArcs[id].x;
            }

            for (int id : inCust[b]) {
                hIn += hc[id];
                gIn += gc[id];
                incomingArcCost += custArcs[id].cost * custArcs[id].x;
            }

            for (int id : outCust[b]) {
                hOut += hc[id];
                gOut += gc[id];
            }

            // Every customer consumes one unit of connectivity flow.
            model.addConstr(
                hIn - hOut == 1.0,
                "unit_flow_balance_" + std::to_string(customers[b])
            );

            // Every customer consumes the cost of the selected incoming arc.
            model.addConstr(
                gIn - gOut == incomingArcCost,
                "cost_flow_balance_" + std::to_string(customers[b])
            );
        }

        // ------------------------------------------------------------------
        // MIP start from initial routes
        // ------------------------------------------------------------------
        if (warmStartValid) {
            std::vector<int> warmY(K_max, 0);
            std::vector<int> warmEnd(C, 0);
            std::vector<double> warmT(C, 0.0);
            std::vector<double> warmOrd(C, 1.0);

            std::vector<double> warmHs(startArcs.size(), 0.0);
            std::vector<double> warmHc(custArcs.size(), 0.0);
            std::vector<double> warmGs(startArcs.size(), 0.0);
            std::vector<double> warmGc(custArcs.size(), 0.0);

            std::vector<int> chosenStartArcIds;
            std::vector<int> chosenCustArcIds;

            std::unordered_map<int, std::vector<int>> startToSlots;
            for (int k = 0; k < K_max; ++k) {
                startToSlots[S[k]].push_back(k);
            }

            std::vector<char> usedSlot(K_max, 0);
            bool canInstallStart = true;
            int installedRoutes = 0;

            for (const auto& route : initialRoutes) {
                const int routeStart = route.front();

                auto it = startToSlots.find(routeStart);
                if (it == startToSlots.end()) {
                    canInstallStart = false;
                    break;
                }

                int k = -1;
                for (int cand : it->second) {
                    if (!usedSlot[cand]) {
                        k = cand;
                        break;
                    }
                }

                if (k < 0) {
                    canInstallStart = false;
                    break;
                }

                usedSlot[k] = 1;
                warmY[k] = 1;
                ++installedRoutes;

                std::vector<int> routeCust;
                routeCust.reserve(route.size() - 1);

                for (int p = 1; p < (int)route.size(); ++p) {
                    routeCust.push_back(custPos[route[p]]);
                }

                const int m = static_cast<int>(routeCust.size());
                if (m <= 0) {
                    canInstallStart = false;
                    break;
                }

                const int firstC = routeCust.front();
                const int sid = startArcId[k][firstC];

                if (sid < 0) {
                    std::cerr << "[compact-flow ILP MIP start] missing start arc "
                              << S[k] << " -> " << customers[firstC] << ".\n";
                    canInstallStart = false;
                    break;
                }

                std::vector<std::pair<int, double>> custArcAndCumBefore;
                custArcAndCumBefore.reserve(m > 0 ? static_cast<size_t>(m - 1) : 0);

                double cumulative = startArcs[sid].cost;

                if (cumulative > B + EPS) {
                    canInstallStart = false;
                    break;
                }

                chosenStartArcIds.push_back(sid);
                warmT[firstC] = cumulative;
                warmOrd[firstC] = 1.0;

                for (int p = 1; p < m; ++p) {
                    const int prevC = routeCust[p - 1];
                    const int curC  = routeCust[p];

                    const int cid = custArcId[prevC][curC];

                    if (cid < 0) {
                        std::cerr << "[compact-flow ILP MIP start] missing customer arc "
                                  << customers[prevC] << " -> "
                                  << customers[curC] << ".\n";
                        canInstallStart = false;
                        break;
                    }

                    const double cumulativeBefore = cumulative;
                    cumulative += custArcs[cid].cost;

                    if (cumulative > B + EPS) {
                        canInstallStart = false;
                        break;
                    }

                    chosenCustArcIds.push_back(cid);
                    custArcAndCumBefore.push_back({cid, cumulativeBefore});

                    warmT[curC] = cumulative;
                    warmOrd[curC] = static_cast<double>(p + 1);
                }

                if (!canInstallStart) break;

                const double totalRouteCost = cumulative;

                warmHs[sid] = static_cast<double>(m);
                warmGs[sid] = totalRouteCost;

                for (int p = 1; p < m; ++p) {
                    const int cid = custArcAndCumBefore[p - 1].first;
                    const double cumulativeBefore = custArcAndCumBefore[p - 1].second;

                    // Arc from customer p-1 to p carries customers p..m-1.
                    warmHc[cid] = static_cast<double>(m - p);

                    // It carries the remaining route cost from that arc onward.
                    warmGc[cid] = totalRouteCost - cumulativeBefore;
                }

                const int lastC = routeCust.back();
                warmEnd[lastC] = 1;
            }

            if (canInstallStart) {
                std::cout << "[compact-flow ILP MIP start] installing "
                          << installedRoutes << " greedy routes.\n";

                // Add a valid greedy upper bound.
                model.addConstr(
                    obj <= static_cast<double>(installedRoutes),
                    "GreedyUB"
                );

                for (int k = 0; k < K_max; ++k) {
                    y[k].set(GRB_DoubleAttr_Start,
                             static_cast<double>(warmY[k]));
                }

                for (int c = 0; c < C; ++c) {
                    t[c].set(GRB_DoubleAttr_Start, warmT[c]);
                    ord[c].set(GRB_DoubleAttr_Start, warmOrd[c]);
                    endVar[c].set(GRB_DoubleAttr_Start,
                                  static_cast<double>(warmEnd[c]));
                }

                for (auto& a : startArcs) {
                    a.x.set(GRB_DoubleAttr_Start, 0.0);
                }

                for (auto& e : custArcs) {
                    e.x.set(GRB_DoubleAttr_Start, 0.0);
                }

                for (int id : chosenStartArcIds) {
                    startArcs[id].x.set(GRB_DoubleAttr_Start, 1.0);
                }

                for (int id : chosenCustArcIds) {
                    custArcs[id].x.set(GRB_DoubleAttr_Start, 1.0);
                }

                for (int id = 0; id < (int)startArcs.size(); ++id) {
                    hs[id].set(GRB_DoubleAttr_Start, warmHs[id]);
                    gs[id].set(GRB_DoubleAttr_Start, warmGs[id]);
                }

                for (int id = 0; id < (int)custArcs.size(); ++id) {
                    hc[id].set(GRB_DoubleAttr_Start, warmHc[id]);
                    gc[id].set(GRB_DoubleAttr_Start, warmGc[id]);
                }
            } else {
                std::cerr << "[compact-flow ILP MIP start] greedy routes validated "
                          << "as paths but could not be installed in the pruned "
                          << "compact arc set. Solving without full MIP start.\n";
                warmStartValid = false;
                fallbackRoutes.clear();
            }
        }

        // ------------------------------------------------------------------
        // Solve
        // ------------------------------------------------------------------
        model.update();

        std::cout << "[compact-flow ILP stats] NumVars="
                  << model.get(GRB_IntAttr_NumVars)
                  << ", NumConstrs="
                  << model.get(GRB_IntAttr_NumConstrs)
                  << ", NumBinVars="
                  << model.get(GRB_IntAttr_NumBinVars)
                  << ", DNumNZs="
                  << model.get(GRB_DoubleAttr_DNumNZs)
                  << "\n";

        bool optimizeReturnedNormally = false;

        try {
            model.optimize();
            optimizeReturnedNormally = true;
        } catch (GRBException& e) {
            std::cerr << "[compact-flow ILP] optimize() raised Gurobi error "
                      << e.getErrorCode()
                      << ": "
                      << e.getMessage()
                      << "\n";
            std::cerr << "[compact-flow ILP] Trying to read available incumbent/bound attributes.\n";
        }

        int status = -1;
        int solCount = 0;
        double incumbentObj = std::numeric_limits<double>::quiet_NaN();

        try {
            status = model.get(GRB_IntAttr_Status);
        } catch (GRBException& e) {
            std::cerr << "[compact-flow ILP] Status unavailable: "
                      << e.getMessage() << "\n";
        }

        try {
            solCount = model.get(GRB_IntAttr_SolCount);
        } catch (GRBException& e) {
            std::cerr << "[compact-flow ILP] SolCount unavailable: "
                      << e.getMessage() << "\n";
            solCount = 0;
        }

        try {
            bestObjBound = model.get(GRB_DoubleAttr_ObjBound);
        } catch (GRBException& e) {
            std::cerr << "[compact-flow ILP] ObjBound unavailable: "
                      << e.getMessage() << "\n";
            bestObjBound = std::numeric_limits<double>::quiet_NaN();
        }

        if (solCount > 0) {
            try {
                incumbentObj = model.get(GRB_DoubleAttr_ObjVal);
            } catch (GRBException& e) {
                std::cerr << "[compact-flow ILP] ObjVal unavailable despite SolCount > 0: "
                          << e.getMessage() << "\n";
            }
        }

        std::cout << "[compact-flow ILP] optimizeReturnedNormally="
                  << (optimizeReturnedNormally ? "true" : "false")
                  << ", status=" << status
                  << ", SolCount=" << solCount
                  << ", ObjBound=" << bestObjBound;

        if (solCount > 0) {
            std::cout << ", ObjVal=" << incumbentObj;
        }

        std::cout << "\n";

        if (status == GRB_INFEASIBLE) {
            std::cerr << "[compact-flow ILP] Model reported infeasible. Computing IIS...\n";

            try {
                model.computeIIS();
                model.write("rooted_open_compact_flow_iis.ilp");
                std::cerr << "[compact-flow ILP] Wrote IIS to rooted_open_compact_flow_iis.ilp\n";
            } catch (GRBException& e) {
                std::cerr << "[compact-flow ILP] IIS failed: "
                          << e.getMessage() << "\n";
            }

            return {routes, bestObjBound};
        }

        // If Gurobi has no incumbent but the validated greedy solution exists,
        // return the validated feasible fallback. The bound is still the
        // Gurobi ObjBound if available.
        if (solCount == 0) {
            std::cerr << "[compact-flow ILP] No incumbent available from Gurobi. "
                      << "status=" << status
                      << ", ObjBound=" << bestObjBound << "\n";

            if (warmStartValid) {
                std::cerr << "[compact-flow ILP] Returning validated greedy fallback routes.\n";
                return {fallbackRoutes, bestObjBound};
            }

            return {routes, bestObjBound};
        }

        // ------------------------------------------------------------------
        // Extract incumbent routes
        // ------------------------------------------------------------------
        routes.clear();

        for (int k = 0; k < K_max; ++k) {
            if (y[k].get(GRB_DoubleAttr_X) < 0.5) continue;

            int first = -1;

            for (int id : startOut[k]) {
                if (startArcs[id].x.get(GRB_DoubleAttr_X) > 0.5) {
                    first = startArcs[id].to;
                    break;
                }
            }

            if (first < 0) {
                std::cerr << "[compact-flow ILP extract] active start slot "
                          << k
                          << " has no selected first customer.\n";
                continue;
            }

            std::vector<int> path;
            path.push_back(S[k]);

            std::vector<char> seen(C, 0);
            int cur = first;
            bool reachedTerminal = false;

            for (int safety = 0; safety <= C; ++safety) {
                if (cur < 0 || cur >= C) {
                    std::cerr << "[compact-flow ILP extract] invalid customer index.\n";
                    break;
                }

                if (seen[cur]) {
                    std::cerr << "[compact-flow ILP extract] repeated customer "
                              << customers[cur] << ".\n";
                    break;
                }

                seen[cur] = 1;
                path.push_back(customers[cur]);

                if (endVar[cur].get(GRB_DoubleAttr_X) > 0.5) {
                    reachedTerminal = true;
                    break;
                }

                int next = -1;

                for (int id : outCust[cur]) {
                    if (custArcs[id].x.get(GRB_DoubleAttr_X) > 0.5) {
                        next = custArcs[id].to;
                        break;
                    }
                }

                if (next < 0) {
                    std::cerr << "[compact-flow ILP extract] customer "
                              << customers[cur]
                              << " is not terminal but has no selected successor.\n";
                    break;
                }

                cur = next;
            }

            if (!reachedTerminal) {
                std::cerr << "[compact-flow ILP extract] path from start slot "
                          << k << " did not reach a terminal customer.\n";
            }

            routes.push_back(std::move(path));
        }

        bool extractedOK = validateRoutes(routes, "compact-flow ILP extracted");

        if (!extractedOK) {
            std::cerr << "[compact-flow ILP] Extracted incumbent failed validation.\n";

            if (warmStartValid) {
                std::cerr << "[compact-flow ILP] Returning validated greedy fallback instead.\n";
                return {fallbackRoutes, bestObjBound};
            }

            routes.clear();
            return {routes, bestObjBound};
        }

        return {routes, bestObjBound};

    } catch (GRBException& e) {
        std::cerr << "[compact-flow ILP] Gurobi error "
                  << e.getErrorCode()
                  << ": "
                  << e.getMessage()
                  << "\n";

        if (warmStartValid) {
            std::cerr << "[compact-flow ILP] Returning validated greedy fallback after exception.\n";
            return {fallbackRoutes, bestObjBound};
        }

        return {routes, bestObjBound};

    } catch (const std::exception& e) {
        std::cerr << "[compact-flow ILP] Error: "
                  << e.what()
                  << "\n";

        if (warmStartValid) {
            std::cerr << "[compact-flow ILP] Returning validated greedy fallback after exception.\n";
            return {fallbackRoutes, bestObjBound};
        }

        return {routes, bestObjBound};
    }
}


// wrappers
std::pair<std::vector<std::vector<int>>, double>
MinKCoverUnrooted_ILP(double budget,
                      const std::vector<std::vector<double>>& slew_costs,
                      const std::vector<double>& dwell_costs,
                      double mipGap, double timeLimit,
                      const std::vector<std::vector<int>>& initialRoutes) {

    int N = (int)slew_costs.size();

    for (int i = 0; i < N; ++i) {
        if (dwell_costs[i] > budget + 1e-9) {
            std::cerr << "Immediately infeasible: node " << i
                      << " has dwell " << dwell_costs[i]
                      << " > budget " << budget << "\n";
            return {{}, 0.0};
        }
    }

    // int N = (int)slew_costs.size();
    ProblemDataMinKUnrooted data(N, N, budget, slew_costs, dwell_costs);

    return solveMinRoutesCoverUnrooted(data, mipGap, timeLimit, initialRoutes);
}


std::pair<std::vector<std::vector<int>>, double>
MinKCoverRootedOpen_ILP(const std::vector<int>& starts, double budget,
                        const std::vector<std::vector<double>>& slew_costs,
                        const std::vector<double>& dwell_costs,
                        double mipGap, double timeLimit,
                        const std::vector<std::vector<int>>& initialRoutes) {

    int N = (int)slew_costs.size();
    int K = (int)starts.size();
    ProblemDataMinKRootedOpen data(N, K, starts, budget, slew_costs, dwell_costs);

    return solveMinRoutesCoverRootedOpen(data, mipGap, timeLimit, initialRoutes);
}
