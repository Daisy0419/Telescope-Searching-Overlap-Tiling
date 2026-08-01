#include "ILP_cplex_top.h"

#include <set>

inline bool validateRouteST(const std::vector<int>& route, int N, int s, int t,
                            const std::vector<std::vector<double>>& Cost,
                            double Budget) {
    if (route.empty() || route.front() != s || route.back() != t) return false;
    double cost = 0.0;
    std::set<int> seen; seen.insert(s);
    for (size_t i = 0; i + 1 < route.size(); ++i) {
        int a = route[i], b = route[i+1];
        if (a < 0 || a >= N || b < 0 || b >= N || a == b) return false;
        if (seen.count(b) && b != t) return false;
        cost += Cost[a][b];
        seen.insert(b);
    }
    return cost <= Budget + 1e-9;
}


// Solve Team s–t Orienteering with exactly K teams.
// Optional: initialRoutes.size()==K provides MIP starts (ignored if invalid).
std::pair<std::vector<std::vector<int>>, double>
solveTeamST_TOP_CPLEX(const ProblemDataTeam& data,
                      double mipGap, double timeLimitSec,
                      const std::vector<std::vector<int>>& initialRoutes /* or {} */) {
    std::vector<std::vector<int>> routes(data.K);
    double BestBound = 0;

    IloEnv env;
    try {
        const int N = data.N, K = data.K;
        const int s = data.startIndex, t = data.endIndex;
        const double Budget = data.BudgetPerTeam;
        const auto& C = data.Cost;
        const auto& P = data.Prize;
        const int M = std::max(1, N - 1); // SCF capacity

        IloModel model(env);

        // ---------- Variables ----------
        // x[k][i][j] in {0,1} for i!=j
        IloArray<IloArray<IloNumVarArray>> x(env, K);
        // f[k][i][j] in [0,M] for i!=j
        IloArray<IloArray<IloNumVarArray>> f(env, K);
        for (int k = 0; k < K; ++k) {
            x[k] = IloArray<IloNumVarArray>(env, N);
            f[k] = IloArray<IloNumVarArray>(env, N);
            for (int i = 0; i < N; ++i) {
                x[k][i] = IloNumVarArray(env, N);
                f[k][i] = IloNumVarArray(env, N);
                for (int j = 0; j < N; ++j) {
                    if (i == j) {
                        x[k][i][j] = IloNumVar(env, 0.0, 0.0, ILOINT); // fixed 0
                        f[k][i][j] = IloNumVar(env, 0.0, 0.0, ILOFLOAT);
                    } else {
                        x[k][i][j] = IloNumVar(env, 0.0, 1.0, ILOINT);
                        f[k][i][j] = IloNumVar(env, 0.0, M,   ILOFLOAT);
                    }
                }
            }
        }

        // z[k][i] in {0,1}
        IloArray<IloNumVarArray> z(env, K);
        for (int k = 0; k < K; ++k) {
            z[k] = IloNumVarArray(env, N);
            for (int i = 0; i < N; ++i)
                z[k][i] = IloNumVar(env, 0.0, 1.0, ILOINT);
        }

        // y[i] in {0,1} with prize in objective
        IloNumVarArray y(env, N);
        for (int i = 0; i < N; ++i) {
            y[i] = IloNumVar(env, 0.0, 1.0, ILOINT);
        }

        // ---------- Objective: maximize sum Prize[i]*y[i] (exclude s,t if desired) ----------
        IloExpr obj(env);
        for (int i = 0; i < N; ++i)
            if (i != s && i != t) obj += P[i] * y[i];
        model.add(IloMaximize(env, obj));
        obj.end();

        // ---------- Constraints ----------

        // One s->t path per team (degree constraints)
        for (int k = 0; k < K; ++k) {
            // start: sum_out == 1, sum_in == 0
            IloExpr startOut(env), startIn(env);
            for (int j = 0; j < N; ++j) if (j != s) startOut += x[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) startIn  += x[k][i][s];
            model.add(startOut == 1.0);
            model.add(startIn  == 0.0);
            startOut.end(); startIn.end();

            // end: sum_in == 1, sum_out == 0
            IloExpr endIn(env), endOut(env);
            for (int i = 0; i < N; ++i) if (i != t) endIn  += x[k][i][t];
            for (int j = 0; j < N; ++j) if (j != t) endOut += x[k][t][j];
            model.add(endIn  == 1.0);
            model.add(endOut == 0.0);
            endIn.end(); endOut.end();

            // intermediate nodes: in == out
            for (int v = 0; v < N; ++v) if (v != s && v != t) {
                IloExpr out(env), in(env);
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.add(out - in == 0.0);
                out.end(); in.end();
            }
        }

        // Define z from x (for v != s,t): sum_out == z and sum_in == z
        for (int k = 0; k < K; ++k) {
            for (int v = 0; v < N; ++v) if (v != s && v != t) {
                IloExpr out(env), in(env);
                for (int j = 0; j < N; ++j) if (j != v) out += x[k][v][j];
                for (int i = 0; i < N; ++i) if (i != v) in  += x[k][i][v];
                model.add(out == z[k][v]);
                model.add(in  == z[k][v]);
                out.end(); in.end();
            }
            // every team uses s,t
            model.add(z[k][s] == 1.0);
            model.add(z[k][t] == 1.0);
        }

        // No duplicate visits across teams (intermediate nodes)
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            IloExpr sumz(env);
            for (int k = 0; k < K; ++k) sumz += z[k][v];
            model.add(sumz <= 1.0);
            sumz.end();
        }

        // Link y to z: y[v] <= sum_k z[k][v], and y[v] >= z[k][v]
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            IloExpr sumz(env);
            for (int k = 0; k < K; ++k) sumz += z[k][v];
            model.add(y[v] <= sumz);
            sumz.end();
            for (int k = 0; k < K; ++k) model.add(y[v] >= z[k][v]);
        }

        // Budget per team
        for (int k = 0; k < K; ++k) {
            IloExpr cost(env);
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    cost += C[i][j] * x[k][i][j];
            model.add(cost <= Budget);
            cost.end();
        }

        // SCF connectivity: capacity + supply/demand
        for (int k = 0; k < K; ++k) {
            // capacity f <= M * x
            for (int i = 0; i < N; ++i)
                for (int j = 0; j < N; ++j) if (i != j)
                    model.add(f[k][i][j] <= M * x[k][i][j]);

            // supply at start = sum z[k][v] (v != s)
            IloExpr fOutS(env), fInS(env), visitsExclStart(env);
            for (int j = 0; j < N; ++j) if (j != s) fOutS += f[k][s][j];
            for (int i = 0; i < N; ++i) if (i != s) fInS  += f[k][i][s];
            for (int v = 0; v < N; ++v) if (v != s) visitsExclStart += z[k][v];
            model.add(fOutS - fInS == visitsExclStart);
            fOutS.end(); fInS.end(); visitsExclStart.end();

            // demand at v != s equals z[k][v]
            for (int v = 0; v < N; ++v) if (v != s) {
                IloExpr fIn(env), fOut(env);
                for (int i = 0; i < N; ++i) if (i != v) fIn  += f[k][i][v];
                for (int j = 0; j < N; ++j) if (j != v) fOut += f[k][v][j];
                model.add(fIn - fOut == z[k][v]);
                fIn.end(); fOut.end();
            }
        }

        // ---------- CPLEX solve ----------
        IloCplex cplex(model);
        if (mipGap >= 0)   cplex.setParam(IloCplex::EpGap, mipGap);
        if (timeLimitSec > 0) cplex.setParam(IloCplex::TiLim, timeLimitSec);
        // cplex.setParam(IloCplex::Threads, 1); // if you want determinism

        // ---------- (Optional) MIP starts from routes ----------
        if (!initialRoutes.empty()) {
            IloNumVarArray startVars(env);
            IloNumArray    startVals(env);

            auto pushStart = [&](IloNumVar v, double val){
                startVars.add(v); startVals.add(val);
            };

            for (int k = 0; k < K && k < (int)initialRoutes.size(); ++k) {
                const auto& r = initialRoutes[k];
                if (!validateRouteST(r, N, s, t, C, Budget)) {
                    std::cerr << "CPLEX: MIP start for team " << k << " ignored (invalid).\n";
                    continue;
                }
                // x and z along the path
                std::set<int> used(r.begin(), r.end());
                for (size_t p = 0; p + 1 < r.size(); ++p) {
                    int i = r[p], j = r[p+1];
                    pushStart(x[k][i][j], 1.0);
                }
                for (int v : used) pushStart(z[k][v], 1.0);
                // rough flow guess
                for (size_t p = 0; p + 1 < r.size(); ++p) {
                    int i = r[p], j = r[p+1];
                    pushStart(f[k][i][j], 1.0);
                }
            }
            if (startVars.getSize() > 0) cplex.addMIPStart(startVars, startVals);
            startVars.end(); startVals.end();
        }

        bool ok = cplex.solve();
        if (!ok) {
            std::cerr << "CPLEX failed to optimize.\n";
            env.end();
            return {routes, 0.0};
        }

        BestBound = cplex.getBestObjValue();
        std::cout << "CPLEX status: " << cplex.getStatus() << "\n";
        std::cout << "Obj: " << cplex.getObjValue()
                  << "  BestBound: " << BestBound
                  << "  Gap: " << cplex.getMIPRelativeGap() << "\n";

        // ---------- Extract routes ----------
        for (int k = 0; k < K; ++k) {
            std::vector<int> path; path.reserve(N);
            int cur = s;
            path.push_back(cur);
            while (cur != t) {
                bool found = false;
                for (int j = 0; j < N; ++j) if (j != cur) {
                    if (cplex.getValue(x[k][cur][j]) > 0.5) {
                        path.push_back(j);
                        cur = j;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    std::cerr << "Warning: route reconstruction failed for team " << k << "\n";
                    break;
                }
            }
            routes[k] = std::move(path);
        }
    } catch (const IloException& e) {
        std::cerr << "CPLEX exception: " << e << "\n";
    } catch (const std::exception& e) {
        std::cerr << "std::exception: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }
    env.end();
    return {routes, BestBound};
}

std::pair<std::vector<std::vector<int>>, double>
solveTeamST_TOP_2index_CPX(const ProblemDataTeam& data,
                           double mipGap, double timeLimitSec,
                           const std::vector<std::vector<int>>& initialRoutes) {
    const int N = data.N;
    const int s = data.startIndex, t = data.endIndex;
    const int m = data.K;                 // number of routes
    const double B = data.BudgetPerTeam;  // per-route budget
    const auto& C = data.Cost;
    const auto& P = data.Prize;

    std::vector<std::vector<int>> routes; routes.reserve(m);
    double bestBound = 0.0;

    IloEnv env;
    try {
        IloModel model(env);

        // Precompute T_from0, T_rem
        std::vector<double> T_from0(N, 0.0), T_rem(N, 0.0);
        for (int i = 0; i < N; ++i) {
            T_from0[i] = C[s][i];
            T_rem[i]   = B - C[i][t];
        }

        // ---------- Variables ----------
        // x[i][j], y[i], z[i][j]
        IloArray<IloNumVarArray> x(env, N);
        IloArray<IloNumVarArray> z(env, N);
        IloNumVarArray y(env, N);

        for (int i = 0; i < N; ++i) {
            x[i] = IloNumVarArray(env, N);
            z[i] = IloNumVarArray(env, N);
            const double obj = (i == s || i == t) ? 0.0 : P[i];
            y[i] = IloNumVar(env, 0.0, 1.0, ILOINT, ("y_" + std::to_string(i)).c_str());
            // y[i].setObjectCoef(obj);
            for (int j = 0; j < N; ++j) {
                if (i == j) {
                    x[i][j] = IloNumVar(env, 0.0, 0.0, ILOINT, ("x_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
                    z[i][j] = IloNumVar(env, 0.0, 0.0, ILOFLOAT, ("z_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
                } else {
                    x[i][j] = IloNumVar(env, 0.0, 1.0, ILOINT, ("x_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
                    z[i][j] = IloNumVar(env, 0.0, IloInfinity, ILOFLOAT, ("z_" + std::to_string(i) + "_" + std::to_string(j)).c_str());
                }
            }
        }

        // // ---------- Objective ----------
        // model.add(IloMaximize(env, IloScalProd(y, IloNumArray(env, N, [&](IloEnv e, int idx){ return (idx==s||idx==t)?0.0:P[idx]; })))); // helper lambda not valid in old compilers
        // // If the above lambda annoys your compiler, the simple loop version used in your code is fine.

        // ---------- Objective ----------
        IloExpr obj(env);
        for (int i = 0; i < N; ++i) {
            if (i != s && i != t) obj += P[i] * y[i];
        }
        model.add(IloMaximize(env, obj));
        obj.end();

        model.add(y[s] == 0.0);
        model.add(y[t] == 0.0);

        // ---------- Degree constraints ----------
        // Start
        {
            IloExpr outS(env), inS(env);
            for (int j = 0; j < N; ++j) if (j != s) outS += x[s][j];
            for (int i = 0; i < N; ++i) if (i != s) inS  += x[i][s];
            model.add(outS == m); outS.end();
            model.add(inS  == 0.0); inS.end();
        }
        // End
        {
            IloExpr inT(env), outT(env);
            for (int i = 0; i < N; ++i) if (i != t) inT  += x[i][t];
            for (int j = 0; j < N; ++j) if (j != t) outT += x[t][j];
            model.add(inT  == m); inT.end();
            model.add(outT == 0.0); outT.end();
        }
        // Customers
        for (int v = 0; v < N; ++v) if (v != s && v != t) {
            IloExpr inV(env), outV(env);
            for (int i = 0; i < N; ++i) if (i != v) inV  += x[i][v];
            for (int j = 0; j < N; ++j) if (j != v) outV += x[v][j];
            model.add(inV  == y[v]); inV.end();
            model.add(outV == y[v]); outV.end();
        }

        // global time bound
        {
            IloExpr tot(env);
            for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j)
                tot += C[i][j] * x[i][j];
            model.add(tot <= m * B);
            tot.end();
        }


        // // Optional: forbid empty routes
        // if (/*forbidEmptyRoutes*/ true) {
        //     model.add(x[s][t] == 0.0);
        // }

        // ---------- Time / arrival propagation ----------
        // (A) z[s][j] = C[s][j] * x[s][j]
        for (int j = 0; j < N; ++j) if (j != s) {
            model.add(z[s][j] == C[s][j] * x[s][j]);
        }
        // (B) time-flow at customers
        for (int i = 0; i < N; ++i) if (i != s && i != t) {
            IloExpr sumOutZ(env), sumInZ(env), rhs(env);
            for (int j = 0; j < N; ++j) if (j != i) { sumOutZ += z[i][j]; rhs += C[i][j] * x[i][j]; }
            for (int h = 0; h < N; ++h) if (h != i) { sumInZ  += z[h][i]; }
            model.add(sumOutZ - sumInZ == rhs);
            sumOutZ.end(); sumInZ.end(); rhs.end();
        }
        // (C) UB: z[i][j] <= (B - C[j][t]) x[i][j] = T_rem[j] x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.add(z[i][j] <= T_rem[j] * x[i][j]);
        }
        // (D) LB: z[i][j] >= (C[s][i] + C[i][j]) x[i][j] = (T_from0[i] + C[i][j]) x[i][j]
        for (int i = 0; i < N; ++i) for (int j = 0; j < N; ++j) if (i != j) {
            model.add(z[i][j] >= (T_from0[i] + C[i][j]) * x[i][j]);
        }

        IloCplex cplex(model);
        if (mipGap >= 0)
            cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, mipGap);
        if (timeLimitSec > 0)
            cplex.setParam(IloCplex::Param::TimeLimit, timeLimitSec);
        // cplex.setParam(IloCplex::Param::Threads, 1); // if desired


        // ---------- Optional MIP starts (aggregate) ----------
        if (!initialRoutes.empty()) {
            IloNumVarArray startVars(env);
            IloNumArray    startVals(env);

            auto push = [&](IloNumVar v, double val){ startVars.add(v); startVals.add(val); };

            for (const auto& r : initialRoutes) {
                if (r.empty() || r.front() != s || r.back() != t) continue;
                double acc = 0.0;
                for (size_t p = 0; p + 1 < r.size(); ++p) {
                    int i = r[p], j = r[p+1];
                    if (i == j) continue;
                    push(x[i][j], 1.0);
                    if (i != s && i != t) push(y[i], 1.0);
                    double zij = (i == s) ? C[s][j] : (acc + C[i][j]);
                    push(z[i][j], zij);
                    acc = zij;
                }
                for (size_t p = 1; p + 1 < r.size(); ++p) {
                    int v = r[p];
                    if (v != s && v != t) push(y[v], 1.0);
                }
            }
            if (startVars.getSize() > 0) cplex.addMIPStart(startVars, startVals);
        }

        // Tell CPLEX we’ll use classic callbacks (safer with user/lazy cuts)
        cplex.setParam(IloCplex::Param::MIP::Strategy::Search, IloCplex::Traditional);

        // Install callbacks
        cplex.use(new (env) TopUserCutCB(env, x, y, N, s, t));
        cplex.use(new (env) TopLazyCB   (env, x, y, N, s, t));

        bool ok = cplex.solve();
        if (!ok) {
            std::cerr << "CPLEX failed to optimize.\n";
            env.end();
            return {routes, 0.0};
        }

        bestBound = cplex.getBestObjValue();

        // ---------- Extract m routes ----------
        std::vector<std::vector<char>> used(N, std::vector<char>(N, 0));
        for (int rIdx = 0; rIdx < m; ++rIdx) {
            std::vector<int> path; path.reserve(N);
            int cur = s; path.push_back(cur);
            while (cur != t) {
                bool found = false;
                for (int j = 0; j < N; ++j) if (j != cur && !used[cur][j]) {
                    if (cplex.getValue(x[cur][j]) > 0.5) {
                        used[cur][j] = 1;
                        cur = j;
                        path.push_back(cur);
                        found = true;
                        break;
                    }
                }
                if (!found) break;
            }
            if (path.back() == t) routes.push_back(std::move(path));
        }
    } catch (const IloException& e) {
        std::cerr << "CPLEX exception: " << e << "\n";
    } catch (const std::exception& e) {
        std::cerr << "std::exception: " << e.what() << "\n";
    } catch (...) {
        std::cerr << "Unknown error.\n";
    }
    env.end();
    return {routes, bestBound};
}


// wrapper
std::pair<std::vector<std::vector<int>>, double>
cplexSolveTeamST(const std::vector<std::vector<double>>& Cost,
                 const std::vector<double>& Prize,
                 int start, int end, double BudgetPerTeam,
                 int mTeams,
                 double mipGap, double timeLimitSec,
                 const std::vector<std::vector<int>>& initialRoutes){
    const int N = static_cast<int>(Prize.size());
    ProblemDataTeam data(N, mTeams, start, end, BudgetPerTeam, Cost, Prize);

    // Solve
    auto [paths, bound] = solveTeamST_TOP_CPLEX(data, mipGap, timeLimitSec, initialRoutes);
    double totalPrize = 0.0;
    std::vector<char> seen(N, 0);
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
cplexSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                 const std::vector<double>& Prize,
                 int start, int end, double BudgetPerTeam,
                 int mTeams,
                 double mipGap, double timeLimitSec,
                 const std::vector<std::vector<int>>& initialRoutes) {
    ProblemDataTeam data((int)Prize.size(), mTeams, start, end, BudgetPerTeam, Cost, Prize);
    // return solveTeamST_TOP_CPLEX(data, mipGap, timeLimitSec, initialRoutes);
    return solveTeamST_TOP_2index_CPX(data, mipGap, timeLimitSec, initialRoutes);
}

