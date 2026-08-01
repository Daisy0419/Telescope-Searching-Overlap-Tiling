#include "gurobi_c++.h"
#include "ILP_gurobi_top.h"
#include <vector>
#include <string>
#include <iostream>
#include <set>
#include <map>
#include <stdexcept>
#include <unordered_set>
#include <algorithm>

inline bool validateRouteST(const std::vector<int>& route, int N, int s, int t,
                            const std::vector<std::vector<double>>& Cost, double Budget) {
    if (route.empty() || route.front() != s || route.back() != t) return false;
    double cost = 0.0;
    std::set<int> seen;
    seen.insert(s);
    for (size_t i = 0; i + 1 < route.size(); ++i) {
        int a = route[i], b = route[i+1];
        if (a < 0 || a >= N || b < 0 || b >= N || a == b) return false;
        if (seen.count(b) && b != t) return false;
        cost += Cost[a][b];
        seen.insert(b);
    }
    return cost <= Budget + 1e-9;
}

// Solve team s-t orienteering with exactly K teams.
// Optional MIP starts: initialRoutes.size()==K, each is an s–t path within budget (others ignored if invalid).
std::pair<std::vector<std::vector<int>>, double>
solveTeamST_TOP(const ProblemDataTeamST& data,
                double mipGap, double timeLimit,
                const std::vector<std::vector<int>>& initialRoutes /* size K or empty */) {
    std::vector<std::vector<int>> routes(data.K);
    double BestBound = 0;

    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_top_st.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "Team_Orienteering_ST");
        if (mipGap >= 0)   model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0) model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        model.set(GRB_DoubleParam_ObjScale, 0.001);

        const int N = data.N, K = data.K;
        const int s = data.startIndex, t = data.endIndex;
        const double Budget = data.BudgetPerTeam;
        const auto& C = data.Cost;
        const auto& P = data.Prize;

        const int M = std::max(1, N - 1); // SCF capacity big-M

        // ---------- Variables ----------
        // x[k][i][j] in {0,1}
        std::vector<std::vector<std::vector<GRBVar>>> x(
            K, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));

        // z[k][i] in {0,1} (team k visits node i)
        std::vector<std::vector<GRBVar>> z(K, std::vector<GRBVar>(N));

        // y[i] in {0,1} (node i visited by any team) – objective counts prize once
        std::vector<GRBVar> y(N);

        // f[k][i][j] >= 0 (single-commodity flow for subtour elimination & connectivity)
        std::vector<std::vector<std::vector<GRBVar>>> f(
            K, std::vector<std::vector<GRBVar>>(N, std::vector<GRBVar>(N)));

        // Create x and f
        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < N; ++j) if (i != j) {
                    // x: objective coefficient 0 (objective uses y, not x)
                    x[k][i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                        "x_" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j));
                    // f
                    f[k][i][j] = model.addVar(0.0, M, 0.0, GRB_CONTINUOUS,
                        "f_" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j));
                }
            }
        }
        // y and z
        for (int i = 0; i < N; ++i) {
            const double obj = (i == s || i == t) ? 0.0 : P[i];
            y[i] = model.addVar(0.0, 1.0, obj, GRB_BINARY, "y_" + std::to_string(i));
            for (int k = 0; k < K; ++k) {
                z[k][i] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                    "z_" + std::to_string(k) + "_" + std::to_string(i));
            }
        }

        // ---------- Constraints ----------
        // One s->t path per team (degree constraints)
        for (int k = 0; k < K; ++k) {
            // start: exactly one outgoing, no incoming
            GRBLinExpr startOut = 0.0, startIn = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) startOut += x[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) startIn  += x[k][i][s];
            model.addConstr(startOut == 1.0, "StartOut_k" + std::to_string(k));
            model.addConstr(startIn  == 0.0, "StartInZero_k" + std::to_string(k));

            // end: exactly one incoming, no outgoing
            GRBLinExpr endIn = 0.0, endOut = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) endIn  += x[k][i][t];
            for (int j = 0; j < N; ++j) if (j != t) endOut += x[k][t][j];
            model.addConstr(endIn  == 1.0, "EndIn_k" + std::to_string(k));
            model.addConstr(endOut == 0.0, "EndOutZero_k" + std::to_string(k));

            // intermediate nodes: in-degree == out-degree (0 if not used)
            for (int v = 0; v < N; ++v) if (v != s && v != t) {
                GRBLinExpr out = 0.0, in = 0.0;
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.addConstr(out - in == 0.0, "DegBal_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
        }

        // Define z from x degrees (for v != s,t): sum_out == z and sum_in == z
        for (int k = 0; k < K; ++k) {
            for (int v = 0; v < N; ++v) if (v != s && v != t) {
                GRBLinExpr out = 0.0, in = 0.0;
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.addConstr(out == z[k][v], "DefZout_k" + std::to_string(k) + "_v" + std::to_string(v));
                model.addConstr(in  == z[k][v], "DefZin_k"  + std::to_string(k) + "_v" + std::to_string(v));
            }
            // every team uses s,t
            model.addConstr(z[k][s] == 1.0, "zStart_k" + std::to_string(k));
            model.addConstr(z[k][t] == 1.0, "zEnd_k"   + std::to_string(k));
        }

        // No duplicate visits across teams for intermediate nodes
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            GRBLinExpr sumz = 0.0;
            for (int k = 0; k < K; ++k) sumz += z[k][v];
            model.addConstr(sumz <= 1.0, "NoDup_v" + std::to_string(v));
        }

        // Link y to team z: y[v] = OR_k z[k][v]  (linearized by bounds)
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            GRBLinExpr sumz = 0.0;
            for (int k = 0; k < K; ++k) sumz += z[k][v];
            model.addConstr(y[v] <= sumz, "yUpper_v" + std::to_string(v));
            for (int k = 0; k < K; ++k)
                model.addConstr(y[v] >= z[k][v], "yLower_k" + std::to_string(k) + "_v" + std::to_string(v));
        }

        // Budget per team
        for (int k = 0; k < K; ++k) {
            GRBLinExpr cost = 0.0;
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    cost += C[i][j] * x[k][i][j];
            model.addConstr(cost <= Budget, "Budget_k" + std::to_string(k));
        }

        // SCF connectivity: flow only on selected arcs & supplies/demands
        for (int k = 0; k < K; ++k) {
            // capacity f <= M * x
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    model.addConstr(f[k][i][j] <= M * x[k][i][j],
                        "FlowCap_k" + std::to_string(k) + "_" + std::to_string(i) + "_" + std::to_string(j));

            // supply at start = number of visited nodes by team k excluding start
            GRBLinExpr fOutS = 0.0, fInS = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) fOutS += f[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) fInS  += f[k][i][s];
            GRBLinExpr visitsExclStart = 0.0;
            for (int v = 0; v < N; ++v) if (v != s) visitsExclStart += z[k][v];
            model.addConstr(fOutS - fInS == visitsExclStart, "FlowSupplyStart_k" + std::to_string(k));

            // demand at all other nodes = z[k][v]
            for (int v = 0; v < N; ++v) if (v != s) {
                GRBLinExpr fIn = 0.0, fOut = 0.0;
                for (int i = 0; i < N; ++i) if (i != v) fIn  += f[k][i][v];
                for (int j = 0; j < N; ++j) if (j != v) fOut += f[k][v][j];
                model.addConstr(fIn - fOut == z[k][v], "FlowDemand_k" + std::to_string(k) + "_v" + std::to_string(v));
            }
        }

        model.update();

        // ---------- Optional MIP starts ----------
        if (!initialRoutes.empty()) {
            for (int k = 0; k < K && k < (int)initialRoutes.size(); ++k) {
                const auto& route = initialRoutes[k];
                if (!validateRouteST(route, N, s, t, C, Budget)) {
                    std::cerr << "Team " << k << " MIP start ignored (invalid route or budget).\n";
                    continue;
                }
                // Set x starts
                for (size_t r = 0; r + 1 < route.size(); ++r) {
                    int i = route[r], j = route[r+1];
                    x[k][i][j].set(GRB_DoubleAttr_Start, 1.0);
                }
                // Set z starts
                std::set<int> used(route.begin(), route.end());
                for (int v : used) z[k][v].set(GRB_DoubleAttr_Start, 1.0);

                // Rough flow start along the path
                for (size_t r = 0; r + 1 < route.size(); ++r) {
                    int i = route[r], j = route[r+1];
                    f[k][i][j].set(GRB_DoubleAttr_Start, 1.0);
                }
                std::cout << "MIP start provided for team " << k << ".\n";
            }
        }

        // ---------- Optimize ----------
        model.optimize();

        const int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            double obj = model.get(GRB_DoubleAttr_ObjVal);
            BestBound = model.get(GRB_DoubleAttr_ObjBound);
            double gap = model.get(GRB_DoubleAttr_MIPGap);
            std::cout << "Best incumbent: " << obj
                      << ", Best bound: " << BestBound
                      << ", GAP: " << gap
                      << ", Status: " << status << "\n";

            // ---------- Extract routes ----------
            for (int k = 0; k < K; ++k) {
                std::vector<int> path;
                path.reserve(N);
                path.push_back(s);
                int cur = s;
                // Follow x[k][cur][j] = 1
                while (cur != t) {
                    bool found = false;
                    for (int j = 0; j < N; ++j) if (j != cur) {
                        if (x[k][cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            path.push_back(j);
                            cur = j;
                            found = true;
                            break;
                        }
                    }
                    if (!found) {
                        std::cerr << "Warning: could not complete path reconstruction for team " << k << "\n";
                        break;
                    }
                }
                routes[k] = std::move(path);
            }
        } else {
            std::cerr << "Solve failed. Status: " << status << "\n";
        }
    } catch (GRBException& e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }
    return {routes, BestBound};
}


std::pair<std::vector<std::vector<int>>, double>
solveTeamST_TOP_2index_GRB(const ProblemDataTeamST& data,
                           double mipGap, double timeLimit,
                           const std::vector<std::vector<int>>& initialRoutes /* size K or empty */) {
    const int N = data.N;
    const int s = data.startIndex, t = data.endIndex;
    const int m = data.K;                 // number of routes
    const double B = data.BudgetPerTeam;  // per-route budget
    const auto& C = data.Cost;            // NxN matrix
    const auto& P = data.Prize;           // size N

    std::vector<std::vector<int>> routes; // will hold m disjoint s->t paths
    routes.reserve(m);
    double bestBound = 0.0;

    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_top_st_2index.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TOP_2index_time");
        if (mipGap >= 0)   model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0) model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        model.set(GRB_DoubleParam_ObjScale, 0.001);

        // Precompute T_from0[i] = C[s][i], T_rem[j] = B - C[j][t]
        std::vector<double> T_from0(N, 0.0), T_rem(N, 0.0);
        for (int i = 0; i < N; ++i) {
            T_from0[i] = C[s][i];
            T_rem[i]   = B - C[i][t];
        }

        // ---------- Variables ----------
        // x[i][j] in {0,1}, y[i] in {0,1}, z[i][j] >= 0
        std::vector<std::vector<GRBVar>> x(N, std::vector<GRBVar>(N));
        std::vector<GRBVar> y(N);
        std::vector<std::vector<GRBVar>> z(N, std::vector<GRBVar>(N));

        for (int i = 0; i < N; ++i) {
            // objective coefficient on y[i] (exclude s,t)
            const double obj = (i == s || i == t) ? 0.0 : P[i];
            y[i] = model.addVar(0.0, 1.0, obj, GRB_BINARY, "y_" + std::to_string(i));
            for (int j = 0; j < N; ++j) if (i != j) {
                x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                        "x_" + std::to_string(i) + "_" + std::to_string(j));
                z[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                                        "z_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        model.addConstr(y[s] == 0.0, "y_start_zero");
        model.addConstr(y[t] == 0.0, "y_end_zero");

        // ---------- Degree constraints ----------
        // Start: exactly m outgoing, no incoming
        {
            GRBLinExpr outS = 0.0, inS = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) outS += x[s][j];
            for (int i = 0; i < N; ++i) if (i != s) inS  += x[i][s];
            model.addConstr(outS == m, "Start_Out");
            model.addConstr(inS  == 0.0, "Start_In_Zero");
        }
        // End: exactly m incoming, no outgoing
        {
            GRBLinExpr inT = 0.0, outT = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) inT  += x[i][t];
            for (int j = 0; j < N; ++j) if (j != t) outT += x[t][j];
            model.addConstr(inT  == m, "End_In");
            model.addConstr(outT == 0.0, "End_Out_Zero");
        }
        // Customers: in-degree == out-degree == y[i]
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            GRBLinExpr inV = 0.0, outV = 0.0;
            for (int i = 0; i < N; ++i) if (i != v) inV  += x[i][v];
            for (int j = 0; j < N; ++j) if (j != v) outV += x[v][j];
            model.addConstr(inV  == y[v], "InDeg_"  + std::to_string(v));
            model.addConstr(outV == y[v], "OutDeg_" + std::to_string(v));
        }
        {
            GRBLinExpr tot = 0.0;
            for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j)
                tot += C[i][j] * x[i][j];
            model.addConstr(tot <= m * B, "GlobalTimeBound");  // OP: <= B
        }


        // // (Optional) forbid direct s->t so every route has >=1 customer
        // if (forbidEmptyRoutes) {
        //     model.addConstr(x[s][t] == 0.0, "No_Direct_s_t");
        // }

        // ---------- Time / arrival propagation (arc-time z) ----------
        // (A) z[s][j] = C[s][j] * x[s][j]
        for (int j = 0; j < N; ++j) if (j != s) {
            model.addConstr(z[s][j] == C[s][j] * x[s][j],
                            "Z_Start_" + std::to_string(j));
        }
        // (B) For each customer i: sum_j z[i][j] - sum_h z[h][i] = sum_j C[i][j] * x[i][j]
        for (int i = 0; i < N; ++i) if (i != s && i != t) {
            GRBLinExpr sumOutZ = 0.0, sumInZ = 0.0, rhs = 0.0;
            for (int j = 0; j < N; ++j) if (j != i) {
                sumOutZ += z[i][j];
                rhs     += C[i][j] * x[i][j];
            }
            for (int h = 0; h < N; ++h) if (h != i) {
                sumInZ += z[h][i];
            }
            model.addConstr(sumOutZ - sumInZ == rhs, "Z_Flow_" + std::to_string(i));
        }
        // (C) Upper bound: z[i][j] <= (B - C[j][t]) * x[i][j] = T_rem[j] * x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.addConstr(z[i][j] <= T_rem[j] * x[i][j],
                            "Z_UB_" + std::to_string(i) + "_" + std::to_string(j));
        }
        // (D) Lower bound: z[i][j] >= (C[s][i] + C[i][j]) * x[i][j] = (T_from0[i] + C[i][j]) * x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.addConstr(z[i][j] >= (T_from0[i] + C[i][j]) * x[i][j],
                            "Z_LB_" + std::to_string(i) + "_" + std::to_string(j));
        }

        model.update();

        // ---------- Optional MIP starts (aggregate) ----------
        // We set x on all arcs used by the provided routes; y on visited customers; z consistent along each path.
        if (!initialRoutes.empty()) {
            // NOTE: If forbidEmptyRoutes==false and you supply empty routes, we'd need x[s][t] as integer. Here we forbid empties.
            for (const auto& r : initialRoutes) {
                // Basic validity check (you may reuse your validateRouteST with m=1 budget B)
                if (r.empty() || r.front() != s || r.back() != t) continue;

                double acc = 0.0; // arrival time at current node
                for (size_t p = 0; p + 1 < r.size(); ++p) {
                    int i = r[p], j = r[p+1];
                    if (i == j) continue;
                    x[i][j].set(GRB_DoubleAttr_Start, 1.0);
                    if (i != s && i != t) y[i].set(GRB_DoubleAttr_Start, 1.0);
                    // set z along path:
                    double zij = (i == s) ? C[s][j] : (acc + C[i][j]);
                    z[i][j].set(GRB_DoubleAttr_Start, zij);
                    acc = zij; // arrival at j
                }
                // mark internal customers
                for (size_t p = 1; p + 1 < r.size(); ++p) {
                    int v = r[p];
                    if (v != s && v != t) y[v].set(GRB_DoubleAttr_Start, 1.0);
                }
            }
        }

        // (Required for addLazy)
        model.set(GRB_IntParam_LazyConstraints, 1);
        BCCallback cb(N, s, t, x, y);
        model.setCallback(&cb);


        // ---------- Optimize ----------
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            bestBound = model.get(GRB_DoubleAttr_ObjBound);

            // ---------- Extract m disjoint s->t routes from x ----------
            // Build a simple successor map; customers have out-degree 0 or 1; s has out-degree m.
            std::vector<std::vector<char>> used(N, std::vector<char>(N, 0));
            for (int rIdx = 0; rIdx < m; ++rIdx) {
                std::vector<int> path; path.reserve(N);
                int cur = s;
                path.push_back(cur);
                while (cur != t) {
                    bool found = false;
                    for (int j = 0; j < N; ++j) if (j != cur && !used[cur][j]) {
                        if (x[cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            used[cur][j] = 1;
                            cur = j;
                            path.push_back(cur);
                            found = true;
                            break;
                        }
                    }
                    if (!found) break; // safety: infeasible extraction
                }
                if (path.back() == t) routes.push_back(std::move(path));
            }
        } else {
            std::cerr << "Solve failed. Status: " << status << "\n";
        }
    } catch (GRBException& e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return {routes, bestBound};
}

std::pair<std::vector<std::vector<int>>, double>
solveTeamST_TOP_2index_multiS(const ProblemDataTeamMultiS& data,
                           double mipGap, double timeLimit,
                           const std::vector<std::vector<int>>& initialRoutes) {
    const int N = data.N;
    const int t = data.endIndex;
    const int m = data.K;                 // number of routes
    const double B = data.BudgetPerTeam;  // per-route budget
    const auto& C = data.Cost;            // NxN matrix
    const auto& P = data.Prize;           // size N
    const auto& starts = data.startIndices;  // size m; duplicates allowed

    std::vector<std::vector<int>> routes; // will hold m disjoint start->t paths
    routes.reserve(m);
    double bestBound = 0.0;

    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_top_st_2index.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TOP_2index_time_multiStart");
        if (mipGap >= 0)   model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0) model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);
        // model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        // model.set(GRB_DoubleParam_ObjScale, 0.001);

        // --- validate starts and build helpers ---
        if ((int)starts.size() != m) {
            std::cerr << "[solve] starts.size() != m\n";
        }
        std::vector<int> startCount(N, 0);
        std::unordered_set<int> startSet;
        for (int s : starts) {
            if (s < 0 || s >= N) {
                std::cerr << "[solve] start node out of range: " << s << "\n";
            }
            if (s == t) {
                std::cerr << "[solve] start equals terminal t=" << t << " (invalid)\n";
            }
            startCount[s] += 1;
            startSet.insert(s);
        }

        // Precompute remaining-time bound wrt t
        std::vector<double> T_rem(N, 0.0);
        for (int i = 0; i < N; ++i) T_rem[i] = B - C[i][t];

        // ---------- Variables ----------
        // x[i][j] in {0,1}, y[i] in {0,1}, z[i][j] >= 0
        std::vector<std::vector<GRBVar>> x(N, std::vector<GRBVar>(N));
        std::vector<GRBVar> y(N);
        std::vector<std::vector<GRBVar>> z(N, std::vector<GRBVar>(N));

        for (int i = 0; i < N; ++i) {
            // objective coefficient on y[i] (exclude all starts & t)
            const double obj = (startSet.count(i) || i == t) ? 0.0 : P[i];
            y[i] = model.addVar(0.0, 1.0, obj, GRB_BINARY, "y_" + std::to_string(i));
            for (int j = 0; j < N; ++j) if (i != j) {
                x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                        "x_" + std::to_string(i) + "_" + std::to_string(j));
                z[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                                        "z_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // starts & terminal cannot collect prize
        for (int s : starts) model.addConstr(y[s] == 0.0, "y_start_zero_" + std::to_string(s));
        model.addConstr(y[t] == 0.0, "y_end_zero");

        // ---------- Degree constraints ----------
        // Starts: for each start node v, exactly startCount[v] outgoing, zero incoming
        for (int v = 0; v < N; ++v) if (startCount[v] > 0) {
            GRBLinExpr outV = 0.0, inV = 0.0;
            for (int j = 0; j < N; ++j) if (j != v) outV += x[v][j];
            for (int i = 0; i < N; ++i) if (i != v) inV  += x[i][v];
            model.addConstr(outV == startCount[v], "Start_Out_" + std::to_string(v));
            model.addConstr(inV  == 0.0,           "Start_In_Zero_" + std::to_string(v));
        }

        // End: exactly m incoming, no outgoing
        {
            GRBLinExpr inT = 0.0, outT = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) inT  += x[i][t];
            for (int j = 0; j < N; ++j) if (j != t) outT += x[t][j];
            model.addConstr(inT  == m, "End_In");
            model.addConstr(outT == 0.0, "End_Out_Zero");
        }

        // Customers (exclude all starts & t): in-degree == out-degree == y[i]
        for (int v = 0; v < N; ++v) if (!startSet.count(v) && v != t) {
            GRBLinExpr inV = 0.0, outV = 0.0;
            for (int i = 0; i < N; ++i) if (i != v) inV  += x[i][v];
            for (int j = 0; j < N; ++j) if (j != v) outV += x[v][j];
            model.addConstr(inV  == y[v], "InDeg_"  + std::to_string(v));
            model.addConstr(outV == y[v], "OutDeg_" + std::to_string(v));
        }

        // Global budget (sum of selected arc costs across all m routes)
        {
            GRBLinExpr tot = 0.0;
            for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j)
                tot += C[i][j] * x[i][j];
            model.addConstr(tot <= m * B, "GlobalTimeBound");
        }

        // // (Optional) forbid direct start->t (to avoid empty routes)
        // for (int s : starts) {
        //     model.addConstr(x[s][t] == 0.0, "No_Direct_start_t_" + std::to_string(s));
        // }

        // ---------- Time / arrival propagation (arc-time z) ----------
        // (A) For each start v: z[v][j] = C[v][j] * x[v][j]
        for (int v = 0; v < N; ++v) if (startCount[v] > 0) {
            for (int j = 0; j < N; ++j) if (j != v) {
                model.addConstr(z[v][j] == C[v][j] * x[v][j],
                                "Z_Start_" + std::to_string(v) + "_" + std::to_string(j));
            }
        }

        // (B) For each non-start, non-terminal i:
        //     sum_j z[i][j] - sum_h z[h][i] = sum_j C[i][j] * x[i][j]
        for (int i = 0; i < N; ++i) if (!startSet.count(i) && i != t) {
            GRBLinExpr sumOutZ = 0.0, sumInZ = 0.0, rhs = 0.0;
            for (int j = 0; j < N; ++j) if (j != i) {
                sumOutZ += z[i][j];
                rhs     += C[i][j] * x[i][j];
            }
            for (int h = 0; h < N; ++h) if (h != i) {
                sumInZ += z[h][i];
            }
            model.addConstr(sumOutZ - sumInZ == rhs, "Z_Flow_" + std::to_string(i));
        }

        // (C) Upper bound: z[i][j] <= (B - C[j][t]) * x[i][j] = T_rem[j] * x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.addConstr(z[i][j] <= T_rem[j] * x[i][j],
                            "Z_UB_" + std::to_string(i) + "_" + std::to_string(j));
        }

        // (D) Generic lower bound (multi-start safe): z[i][j] >= C[i][j] * x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.addConstr(z[i][j] >= C[i][j] * x[i][j],
                            "Z_LB_" + std::to_string(i) + "_" + std::to_string(j));
        }

        model.update();

        // ---------- Optional MIP starts (aggregate) ----------
        if (!initialRoutes.empty()) {
            for (const auto& r : initialRoutes) {
                if (r.size() < 2 || r.back() != t || !startSet.count(r.front())) continue;

                double acc = 0.0; // arrival time at current node
                for (size_t p = 0; p + 1 < r.size(); ++p) {
                    int i = r[p], j = r[p+1];
                    if (i == j) continue;
                    x[i][j].set(GRB_DoubleAttr_Start, 1.0);
                    if (!startSet.count(i) && i != t) y[i].set(GRB_DoubleAttr_Start, 1.0);
                    double zij = (startSet.count(i) ? C[i][j] : (acc + C[i][j]));
                    z[i][j].set(GRB_DoubleAttr_Start, zij);
                    acc = zij; // arrival at j
                }
                for (size_t p = 1; p + 1 < r.size(); ++p) {
                    int v = r[p];
                    if (!startSet.count(v) && v != t) y[v].set(GRB_DoubleAttr_Start, 1.0);
                }
            }
        }

        // (Required for addLazy) — update your callback to handle multiple starts!
        model.set(GRB_IntParam_LazyConstraints, 1);
        // OLD (single-start): BCCallback cb(N, s, t, x, y);
        // NEW (suggested):    BCCallback cb(N, starts, t, x, y);
        // For now, keep your existing line if your callback still compiles:
        // BCCallback cb(N, /*single s*/ (starts.empty()?0:starts[0]), t, x, y);
        // model.setCallback(&cb);

        // ---------- Optimize ----------
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            bestBound = model.get(GRB_DoubleAttr_ObjBound);

            // ---------- Extract m routes, one for each entry of startIndices ----------
            std::vector<std::vector<char>> used(N, std::vector<char>(N, 0));
            for (int rIdx = 0; rIdx < m; ++rIdx) {
                int s = starts[rIdx];
                std::vector<int> path; path.reserve(N);
                int cur = s;
                path.push_back(cur);
                while (cur != t) {
                    bool found = false;
                    for (int j = 0; j < N; ++j) if (j != cur && !used[cur][j]) {
                        if (x[cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                            used[cur][j] = 1;
                            cur = j;
                            path.push_back(cur);
                            found = true;
                            break;
                        }
                    }
                    if (!found) break; // safety: infeasible extraction
                }
                if (!path.empty() && path.back() == t) routes.push_back(std::move(path));
            }
        } else {
            std::cerr << "Solve failed. Status: " << status << "\n";
        }
    } catch (GRBException& e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }

    return {routes, bestBound};
}



//wrapper 
std::pair<std::vector<std::vector<int>>, double>
gurobiSolveTeamST(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  int start, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes) {
    ProblemDataTeamST data((int)Prize.size(), mTeams, start, end, BudgetPerTeam, Cost, Prize);
    auto [paths, bound] = solveTeamST_TOP(data, mipGap, timeLimit, initialRoutes);
    double totalPrize = 0.0;
    std::vector<char> seen(Prize.size(), 0);
    for (const auto& path : paths) {
        for (int v : path) {
            if (v == start || v == end) continue;
            if (!seen[v]) {
                totalPrize += Prize[v];
                seen[v] = 1;
            }
        }
    }
    return {paths, totalPrize};    
}


std::pair<std::vector<std::vector<int>>, double>
gurobiSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  int start, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes /*size mTeams or {}*/) {
    ProblemDataTeamST data((int)Prize.size(), mTeams, start, end, BudgetPerTeam, Cost, Prize);
    // return solveTeamST_TOP(data, mipGap, timeLimit, initialRoutes);
     return solveTeamST_TOP_2index_GRB(data, mipGap, timeLimit, initialRoutes);
}


std::tuple<std::vector<std::vector<int>>, double, double>
gurobiSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  std::vector<int> starts, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes) {
    ProblemDataTeamMultiS data((int)Prize.size(), mTeams, starts, end, BudgetPerTeam, Cost, Prize);
       
    auto [paths, bound] = solveTeamST_TOP_2index_multiS(data, mipGap, timeLimit, initialRoutes);
    double totalPrize = 0.0;
    std::vector<char> seen(Prize.size(), 0);
    for (const auto& path : paths) {
        for (int v : path) {
            if (std::find(starts.begin(), starts.end(), v) != starts.end() || v == end) continue;
            if (!seen[v]) {
                totalPrize += Prize[v];
                seen[v] = 1;
            }
        }
    }
    return {paths, totalPrize, bound};    
    // return solveTeamST_TOP(data, mipGap, timeLimit, initialRoutes);
    //  return solveTeamST_TOP_2index_multiS(data, mipGap, timeLimit, initialRoutes);
}





