#include "ILP_gurobi_op.h"
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <set>  
#include <map> 


std::vector<int> solveWithGurobiST_time(const ProblemDataST& data,
                                        double mipGap, double timeLimit,
                                        const std::vector<int>& initialRoute) {
    std::vector<int> route;
    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_tsp_time.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TSP_ST_Path_TimeContinuous");
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        model.set(GRB_DoubleParam_ObjScale, 0.001);

        if (mipGap >= 0)   model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit > 0) model.set(GRB_DoubleParam_TimeLimit, timeLimit);
        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        const int N      = data.N;
        const int s      = data.startIndex;
        const int t      = data.endIndex;
        const double B   = data.Budget;
        const auto& C    = data.Cost;   // NxN
        const auto& P    = data.Prize;  // N

        // ---------- Variables ----------
        std::vector<std::vector<GRBVar>> x(N, std::vector<GRBVar>(N));
        std::vector<GRBVar>              y(N);
        std::vector<std::vector<GRBVar>> z(N, std::vector<GRBVar>(N));

        // y: visit indicator (no prize at s or t)
        for (int i = 0; i < N; ++i) {
            const double obj = (i == s || i == t) ? 0.0 : P[i];
            y[i] = model.addVar(0.0, 1.0, obj, GRB_BINARY, "y_" + std::to_string(i));
        }

        // x,z on all i!=j
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) if (i != j) {
                x[i][j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY,
                                       "x_" + std::to_string(i) + "_" + std::to_string(j));
                z[i][j] = model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS,
                                       "z_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        // constraints
        // Start: exactly one outgoing, zero incoming
        {
            GRBLinExpr outS = 0.0, inS = 0.0;
            for (int j = 0; j < N; ++j) if (j != s) outS += x[s][j];
            for (int i = 0; i < N; ++i) if (i != s) inS  += x[i][s];
            model.addConstr(outS == 1.0, "Start_Out");
            model.addConstr(inS  == 0.0, "Start_In_Zero");
        }

        // End: exactly one incoming, zero outgoing
        {
            GRBLinExpr inT = 0.0, outT = 0.0;
            for (int i = 0; i < N; ++i) if (i != t) inT  += x[i][t];
            for (int j = 0; j < N; ++j) if (j != t) outT += x[t][j];
            model.addConstr(inT  == 1.0, "End_In");
            model.addConstr(outT == 0.0, "End_Out_Zero");
        }

        // Customers (exclude s,t): in-degree == out-degree == y[i]
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            GRBLinExpr inV = 0.0, outV = 0.0;
            for (int i = 0; i < N; ++i) if (i != v) inV  += x[i][v];
            for (int j = 0; j < N; ++j) if (j != v) outV += x[v][j];
            model.addConstr(inV  == y[v], "InDeg_"  + std::to_string(v));
            model.addConstr(outV == y[v], "OutDeg_" + std::to_string(v));
        }

        // Starts & end cannot collect prize
        model.addConstr(y[s] == 0.0, "y_start_zero");
        model.addConstr(y[t] == 0.0, "y_end_zero");

        // Global budget
        GRBLinExpr tot = 0.0;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) if (i != j)
                tot += C[i][j] * x[i][j];
        model.addConstr(tot <= B, "Budget");

        // Time / arrival propagation 
        // Precompute T_rem[j] = B - C[j][t] (remaining time needed after reaching j)
        std::vector<double> T_rem(N, 0.0);
        for (int j = 0; j < N; ++j) T_rem[j] = B - C[j][t];

        // (A) Start arcs: z[s][j] == C[s][j] * x[s][j]
        for (int j = 0; j < N; ++j) if (j != s)
            model.addConstr(z[s][j] == C[s][j] * x[s][j],
                            "Z_Start_" + std::to_string(s) + "_" + std::to_string(j));

        // (B) Flow at intermediate nodes (v != s,t):
        //     sum_j z[v][j] - sum_h z[h][v] = sum_j C[v][j] x[v][j]
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            GRBLinExpr sumOutZ = 0.0, sumInZ = 0.0, rhs = 0.0;
            for (int j = 0; j < N; ++j) if (j != v) {
                sumOutZ += z[v][j];
                rhs     += C[v][j] * x[v][j];
            }
            for (int h = 0; h < N; ++h) if (h != v) sumInZ += z[h][v];
            model.addConstr(sumOutZ - sumInZ == rhs, "Z_Flow_" + std::to_string(v));
        }

        // (C) Upper bound wrt terminal: z[i][j] ≤ (B - C[j][t]) x[i][j]
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) if (i != j)
                model.addConstr(z[i][j] <= T_rem[j] * x[i][j],
                                "Z_UB_" + std::to_string(i) + "_" + std::to_string(j));

        // (D) Lower bound: z[i][j] ≥ C[i][j] x[i][j]
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) if (i != j)
                model.addConstr(z[i][j] >= C[i][j] * x[i][j],
                                "Z_LB_" + std::to_string(i) + "_" + std::to_string(j));

        model.update();

        // Optional MIP start from initialRoute 
        if (!initialRoute.empty()) {
            bool ok = (initialRoute.front() == s && initialRoute.back() == t);
            if (!ok) {
                std::cerr << "Warning: initialRoute must start at s and end at t; skipping MIP start.\n";
            } else {
                // Check budget and indices; set starts for x,y,z
                double acc = 0.0; // arrival at current node
                double routeCost = 0.0;
                std::vector<char> seen(N, 0);
                seen[s] = 1;

                for (size_t k = 0; k + 1 < initialRoute.size(); ++k) {
                    int i = initialRoute[k], j = initialRoute[k + 1];
                    if (i < 0 || i >= N || j < 0 || j >= N || i == j) { ok = false; break; }
                    routeCost += C[i][j];
                }
                if (routeCost > B) ok = false;

                if (!ok) {
                    std::cerr << "MIP start invalid or over budget; skipping.\n";
                } else {
                    for (size_t k = 0; k + 1 < initialRoute.size(); ++k) {
                        int i = initialRoute[k], j = initialRoute[k + 1];
                        x[i][j].set(GRB_DoubleAttr_Start, 1.0);
                        if (i != s && i != t) y[i].set(GRB_DoubleAttr_Start, 1.0);
                        double zij = (i == s) ? C[i][j] : (acc + C[i][j]);
                        z[i][j].set(GRB_DoubleAttr_Start, zij);
                        acc = zij; // arrival at j
                    }
                    // also mark interior nodes’ y
                    for (size_t k = 1; k + 1 < initialRoute.size(); ++k) {
                        int v = initialRoute[k];
                        if (v != s && v != t) y[v].set(GRB_DoubleAttr_Start, 1.0);
                    }
                    std::cout << "MIP start provided from initialRoute (time-continuous).\n";
                }
            }
        }

        // Optimize
        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            double best_bound = model.get(GRB_DoubleAttr_ObjBound);
            double gap        = model.get(GRB_DoubleAttr_MIPGap);
            std::cout << "Best bound: " << best_bound << "  Gap: " << gap << "\n";

            // Extract s->t path
            route.clear();
            route.push_back(s);
            int cur = s;
            while (cur != t) {
                bool found = false;
                for (int j = 0; j < N; ++j) if (j != cur) {
                    if (x[cur][j].get(GRB_DoubleAttr_X) > 0.5) {
                        route.push_back(j);
                        cur = j;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Warning: Failed to reconstruct complete path.\n";
                    break;
                }
            }
        } else {
            std::cerr << "Solver failed. Status: " << status << "\n";
        }

    } catch (GRBException& e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (...) {
        std::cerr << "Unknown error occurred.\n";
    }
    return route;
}


std::vector<int> gurobiSolveST_time(const std::vector<std::vector<double>>& Cost,
                                    const std::vector<double>& Prize,
                                    int start, int end, double Budget,
                                    double mipGap, double timeLimit,
                                    const std::vector<int>& initialRoute) {
    ProblemDataST data((int)Prize.size(), start, end, Budget, Cost, Prize);
    return solveWithGurobiST_time(data, mipGap, timeLimit, initialRoute);
}


