#pragma once
#include <functional>

// // ---------- generic PSO that calls a splitter directly ----------
// using SplitterFn = std::function<
//     std::pair<std::vector<std::vector<int>>, double>(
//         const Instance&, double /*B*/, const std::vector<int>& /*perm*/, int /*m_routes*/
//     )
// >;

// std::pair<std::vector<std::vector<int>>, double>
// run_pso_team_op(const Instance& I, double B, int m_routes,
//                 int swarm_size, int iters,
//                 double vmax, unsigned seed,
//                 const SplitterFn& splitter);


std::pair<std::vector<std::vector<int>>, double>
run_pso_team_op_multi(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    double budget,
                    int m_routes,
                    const std::vector<int>& starts,
                    int t, 
                    int swarm_size=1000, int iters=1000,
                    double vmax=5, unsigned seed=42);