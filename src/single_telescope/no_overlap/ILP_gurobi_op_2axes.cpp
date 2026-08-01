#include "ILP_gurobi_op.h"
#include <string>
#include <iostream>
#include <cmath>
#include <limits>
#include <set>  
#include <map> 


std::vector<int> solveWithGurobiST(const ProblemDataST& data, double mipGap, double timeLimit,
                                   const std::vector<int>& initialRoute) {
    std::vector<int> route;
    try {
        GRBEnv env(true);
        env.set("LogFile", "gurobi_tsp.log");
        env.start();

        GRBModel model(env);
        model.set(GRB_StringAttr_ModelName, "TSP_ST_Path");
        model.set(GRB_DoubleParam_OptimalityTol, 1e-8);
        model.set(GRB_DoubleParam_ObjScale, 0.001);

        if (mipGap != -1)
            model.set(GRB_DoubleParam_MIPGap, mipGap);
        if (timeLimit != -1)
            model.set(GRB_DoubleParam_TimeLimit, timeLimit);

        int N = data.N;
        int start = data.startIndex;
        int end = data.endIndex;
        double Budget = data.Budget;
        const auto& Cost = data.Cost;
        const auto& Prize = data.Prize;

        std::vector<std::vector<GRBVar>> x(N, std::vector<GRBVar>(N));
        std::vector<GRBVar> u(N);

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                double objCoeff = (i == j) ? 0.0 : Prize[j];
                x[i][j] = model.addVar(0.0, 1.0, objCoeff, GRB_BINARY,
                                       "x_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }

        for (int i = 0; i < N; i++) {
            if (i == start) continue;
            u[i] = model.addVar(2.0, (double)N, 0.0, GRB_INTEGER,
                                "u_" + std::to_string(i));
        }

        model.set(GRB_IntAttr_ModelSense, GRB_MAXIMIZE);

        // Constraints
        // A) At_Most_One_Incoming and B) At_Most_One_Outgoing
        for (int i = 0; i < N; i++) {
            GRBLinExpr lhs_in = 0.0;
            GRBLinExpr lhs_out = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i) {
                    lhs_in += x[j][i];
                    lhs_out += x[i][j];
                }
            }
            model.addConstr(lhs_in <= 1.0, "AtMostOneIn_" + std::to_string(i));
            model.addConstr(lhs_out <= 1.0, "AtMostOneOut_" + std::to_string(i));
        }

        // C)Fixed start node
        GRBLinExpr startOut = 0.0;
        for (int j = 0; j < N; j++) {
            if (j != start) startOut += x[start][j];
        }
        model.addConstr(startOut == 1.0, "StartOut");

        for (int i = 0; i < N; i++) {
            if (i != start) model.addConstr(x[i][start] == 0.0);
        }

        // D) Fixed end node
        GRBLinExpr endIn = 0.0;
        for (int i = 0; i < N; i++) { 
            if (i != end) endIn += x[i][end];
        }
        model.addConstr(endIn == 1.0, "EndIn");

        for (int j = 0; j < N; j++) {
            if (j != end) model.addConstr(x[end][j] == 0.0);
        }

        // E) Flow balance
        for (int i = 0; i < N; i++) {
            if (i == start || i == end) continue;
            GRBLinExpr flow = 0.0;
            for (int j = 0; j < N; j++) {
                if (j != i) flow += x[j][i] - x[i][j];
            }
            model.addConstr(flow == 0.0, "FlowBalance_" + std::to_string(i));
        }

        // F) Budget Constraint
        GRBLinExpr totalCost = 0.0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i != j) totalCost += Cost[i][j] * x[i][j];
            }
        }
        model.addConstr(totalCost <= Budget, "Budget");

        // G) Subtour_Elimination (MTZ)
        for (int i = 0; i < N; i++) {
            if (i == start) continue;
            for (int j = 0; j < N; j++) {
                if (j == start || i == j) continue;
                model.addConstr(u[i] - u[j] + N * x[i][j] <= N - 1,
                                "MTZ_" + std::to_string(i) + "_" + std::to_string(j));
            }
        }


        model.update(); 
        // Set MIP start from initialRoute
        if (!initialRoute.empty()) {
            // Validate initialRoute
            bool isValid = true;
            if (initialRoute.front() != start || initialRoute.back() != end) {
                isValid = false;
                std::cerr << "Warning: initialRoute does not start with startIndex or end with endIndex.\n";
            } else {
                // Check if route forms a valid path and respects budget
                double routeCost = 0.0;
                std::set<int> visited;
                visited.insert(start);
                for (size_t k = 0; k < initialRoute.size() - 1; ++k) {
                    int i = initialRoute[k];
                    int j = initialRoute[k + 1];
                    if (i < 0 || i >= N || j < 0 || j >= N || i == j) {
                        isValid = false;
                        std::cerr << "Warning: Invalid node indices in initialRoute.\n";
                        break;
                    }
                    if (visited.count(j) && j != end) {
                        isValid = false;
                        std::cerr << "Warning: initialRoute contains a cycle before reaching end.\n";
                        break;
                    }
                    routeCost += Cost[i][j];
                    visited.insert(j);
                }
                if (routeCost > Budget) {
                    isValid = false;
                    std::cerr << "Warning: initialRoute exceeds budget (" << routeCost << " > " << Budget << ").\n";
                }
            }

            if (isValid) {

                for (size_t k = 0; k < initialRoute.size() - 1; ++k) {
                    int i = initialRoute[k];
                    int j = initialRoute[k + 1];
                    x[i][j].set(GRB_DoubleAttr_Start, 1.0);
                }

                // Set u[i].Start for MTZ variables
                std::map<int, int> nodeOrder;
                for (size_t k = 0; k < initialRoute.size(); ++k) {
                    nodeOrder[initialRoute[k]] = k + 2; // Start u at 2
                }
                for (int i = 0; i < N; i++) {
                    if (i == start) continue;
                    if (nodeOrder.count(i)) {
                        u[i].set(GRB_DoubleAttr_Start, (double)nodeOrder[i]);
                    } else {
                        u[i].set(GRB_DoubleAttr_Start, GRB_UNDEFINED); 
                    }
                }

                std::cout << "MIP start provided from initialRoute.\n";
            } else {
                std::cerr << "MIP start skipped due to invalid initialRoute.\n";
            }
        }

        model.optimize();

        int status = model.get(GRB_IntAttr_Status);
        if (status == GRB_OPTIMAL || status == GRB_TIME_LIMIT) {
            //print best bound
            double best_bound = model.get(GRB_DoubleAttr_ObjBoundC); 
            double result_gap = model.get(GRB_DoubleAttr_MIPGap);
            std::cout << "Objective value: " << best_bound << " (Status: " << status << ")\n";
            std::cout << "Optimality gap: " << result_gap << std::endl;

            // Reconstruct path
            route.push_back(start);
            int current = start;
            while (current != end) {
                bool foundNext = false;
                for (int j = 0; j < N; j++) {
                    if (j == current) continue;
                    if (x[current][j].get(GRB_DoubleAttr_X) > 0.5) {
                        route.push_back(j);
                        current = j;
                        foundNext = true;
                        break;
                    }
                }
                if (!foundNext) {
                    std::cerr << "Warning: Failed to reconstruct complete path.\n";
                    break;
                }
            }
        } else {
            std::cerr << "Solver failed. Status: " << status << "\n";
        }

    } catch (GRBException &e) {
        std::cerr << "Gurobi error: " << e.getMessage() << "\n";
    } catch (...) {
        std::cerr << "Unknown error occurred.\n";
    }
    return route;
}


std::vector<int> gurobiSolveST(const std::vector<std::vector<double>>& Cost,
                               const std::vector<double>& Prize,
                               int start, int end, double Budget,
                               double mipGap, double timeLimit,
                               const std::vector<int>& initialRoute) {
    int N = Prize.size();
    ProblemDataST data(N, start, end, Budget, Cost, Prize);
    return solveWithGurobiST(data, mipGap, timeLimit, initialRoute);
}

