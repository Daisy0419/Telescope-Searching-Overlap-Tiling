#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>
#include <iostream>

#include "SimulatedAnnealing.h"
#include "helplers.h"
#include "LocalSearch.h"
#include "PathOperations.h"

inline std::vector<std::vector<int>>
construct_paths(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                int start, int end, int nPath,
                double budgetT, double greediness) {

    const double EPS = 1e-12;
    const int n = (int)costs.size();

    if (costs[start][end] > budgetT + EPS) {
        std::vector<std::vector<int>> triv(nPath, std::vector<int>{start, end});
        return triv;
    }

    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double alpha = (greediness < 0.0 ? U01(rng) : std::clamp(greediness, 0.0, 1.0));

    // init routes as [s, t]
    std::vector<std::vector<int>> paths(nPath, std::vector<int>{start, end});

    // per-route costs: c(s,t)
    std::vector<double> pathCosts(nPath, costs[start][end]);

    std::vector<char> used(n, 0);
    used[start] = 1; used[end] = 1;

    struct Cand { int v, r, pos; double dt, h; };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue;           
            for (int r = 0; r < nPath; ++r) {
                auto [pos, dt] = cheapestInsertionOnePath(costs, paths[r], v);
                if (pos == -1) continue;
                if (pathCosts[r] + dt <= budgetT + EPS) {
                    double denom = std::max(dt, EPS);
                    double h = prizes[v] / denom; 
                    C.push_back({v, r, pos, dt, h});
                }
            }
        }
    };

    int steps = 0;
    while (steps++ < 10000) {
        std::vector<Cand> C; buildCandidates(C);
        if (C.empty()) break;

        // Value-based RCL: threshold = minH + alpha * (maxH - minH)
        double minH = C.front().h, maxH = C.front().h;
        for (const auto& c : C) { minH = std::min(minH, c.h); maxH = std::max(maxH, c.h); }
        const double thr = minH + alpha * (maxH - minH);

        std::vector<int> RCL; RCL.reserve(C.size());
        for (int i = 0; i < (int)C.size(); ++i)
            if (C[i].h + EPS >= thr) RCL.push_back(i);
        if (RCL.empty()) break;

        std::uniform_int_distribution<int> pick(0, (int)RCL.size() - 1);
        const Cand& mv = C[RCL[pick(rng)]];

        auto& route = paths[mv.r];
        route.insert(route.begin() + mv.pos, mv.v);
        pathCosts[mv.r] += mv.dt;
        used[mv.v] = 1;

    }

    return paths;
}


// Multi-start constructor: start_indices[r] is the start for route r, shared end.
inline std::vector<std::vector<int>>
construct_paths_multiS(const std::vector<std::vector<double>>& costs,
                       const std::vector<double>& prizes,
                       const std::vector<int>& start_indices,
                       int end,
                       int nPath,               
                       double budgetT,
                       double greediness) {
    const double EPS = 1e-12;
    const int n = (int)costs.size();
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    if ((int)start_indices.size() != nPath) throw std::runtime_error("start_indices.size() != nPath");
    if (end < 0 || end >= n) throw std::runtime_error("bad end index");

    {
        std::vector<char> seen(n, 0);
        for (int s : start_indices) {
            if (s < 0 || s >= n) throw std::runtime_error("bad start index");
            if (s == end) throw std::runtime_error("start cannot equal end");
            if (seen[s]) throw std::runtime_error("duplicate start in start_indices");
            seen[s] = 1;
        }
    }

    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double alpha = (greediness < 0.0 ? U01(rng) : std::clamp(greediness, 0.0, 1.0));

    std::vector<std::vector<int>> paths(nPath);
    std::vector<double> pathCosts(nPath, 0.0);
    for (int r = 0; r < nPath; ++r) {
        int s_r = start_indices[r];
        paths[r] = { s_r, end };
        pathCosts[r] = costs[s_r][end]; 
    }

    std::vector<char> used(n, 0);
    for (int s_r : start_indices) used[s_r] = 1;
    used[end] = 1;

    struct Cand { int v, r, pos; double dt, h; };

    auto cheapestInsertionOnePath = [&](const std::vector<int>& route, int v) -> std::pair<int,double> {
        int best_pos = -1; double best_dt = std::numeric_limits<double>::infinity();
        if (route.size() < 2) return { -1, best_dt };
        for (std::size_t j = 0; j + 1 < route.size(); ++j) {
            int a = route[j], b = route[j+1];
            double dt = costs[a][v] + costs[v][b] - costs[a][b];
            if (std::isnan(dt)) continue;
            if (dt < best_dt) { best_dt = dt; best_pos = (int)j + 1; }
        }
        return { best_pos, best_dt };
    };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue; // only unvisited
            for (int r = 0; r < nPath; ++r) {
                auto [pos, dt] = cheapestInsertionOnePath(paths[r], v);
                if (pos == -1) continue;
                if (pathCosts[r] + dt <= budgetT + EPS) {
                    double denom = std::max(dt, EPS);
                    double h = prizes[v] / denom; // value-per-extra-cost
                    C.push_back({v, r, pos, dt, h});
                }
            }
        }
    };

    int steps = 0;
    while (steps++ < 10000) {
        std::vector<Cand> C; C.reserve(256);
        buildCandidates(C);
        if (C.empty()) break;

        // Value-based RCL: threshold = minH + alpha * (maxH - minH)
        double minH = C.front().h, maxH = C.front().h;
        for (const auto& c : C) { minH = std::min(minH, c.h); maxH = std::max(maxH, c.h); }
        const double thr = minH + alpha * (maxH - minH);

        std::vector<int> RCL; RCL.reserve(C.size());
        for (int i = 0; i < (int)C.size(); ++i)
            if (C[i].h + EPS >= thr) RCL.push_back(i);
        if (RCL.empty()) break;

        std::uniform_int_distribution<int> pick(0, (int)RCL.size() - 1);
        const Cand& mv = C[RCL[pick(rng)]];

        auto& route = paths[mv.r];
        route.insert(route.begin() + mv.pos, mv.v);
        pathCosts[mv.r] += mv.dt;
        used[mv.v] = 1;
    }

    return paths;
}


void top_neighbor(const std::vector<std::vector<int>>& paths,
                  const std::vector<std::vector<double>>& costs,
                  const std::vector<double>& prizes,
                  double budget,
                  std::vector<std::vector<int>>& candidate) {

    thread_local std::mt19937_64 rng{std::random_device{}()};

    candidate = paths;
    std::uniform_int_distribution<int> pick_op(0, 7);

    int op = pick_op(rng);
    auto backup = candidate;
    bool changed = false;

    switch (op) {
        case 0:
            changed = relocate_k(candidate, costs, prizes, 1);
            break;
        case 1:
            changed = swap_1_1(candidate, costs, prizes, true);
            break;
        case 2:
            changed = cross_exchange(candidate, costs, prizes);
            break;
        case 3:
            changed = ejection_chain(candidate, costs, prizes, std::max<int>(1, (int)paths.size()-1));
            break;
        case 4: {
            std::vector<double> pathCosts = multiPathsCost(costs, candidate);
            changed = swapNodesBetweenPaths(costs, candidate, pathCosts, budget);
            break;
        }
        case 5: {
            std::vector<double> pathCosts = multiPathsCost(costs, candidate);
            changed = replaceNodesInPaths(costs, candidate, pathCosts, prizes, budget);
            break;
        }
        case 6: {
            std::vector<double> pathCosts = multiPathsCost(costs, candidate);
            changed = insertNodesToPaths(costs, candidate, pathCosts, prizes, budget);
            break;
        }
        case 7:
            changed = merge_split_greedy_budget(candidate, costs, prizes, budget);
            break;
    }

    if (!changed) {
        candidate = std::move(backup); 
        return;
    }

    {
        std::vector<double> pathCosts = multiPathsCost(costs, candidate);
        twoOptAllRoutes(costs, candidate, pathCosts, budget);
    }
}

std::vector<std::vector<int>> top_annealing(const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& prizes,
                                     double budget, std::vector<int> starts, int end,
                                     std::vector<std::vector<int>> init_solution,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate) {

    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    auto current_solution = std::move(init_solution);
    auto best_solution = current_solution;

    double cur_prize  = multiPathsSumPrize(prizes, current_solution);
    double best_prize = cur_prize;

    double temp = initial_temp;

    for (int iter = 0; iter < iterations; ++iter) {
        std::vector<std::vector<int>> candidate_solution;
        top_neighbor(current_solution, costs, prizes, budget, candidate_solution);

        if (multiPathsMaxCost(costs, candidate_solution) > budget + 1e-12) {
            temp *= cooling_rate;
            if (temp < 1e-6) break;
            continue;
        }

        double cand_prize = multiPathsSumPrize(prizes, candidate_solution);
        double delta = cand_prize - cur_prize;

        if (delta >= 0.0 || std::exp(delta / std::max(1e-12, temp)) > U01(rng)) {
            current_solution = std::move(candidate_solution);
            cur_prize = cand_prize;                      // <— update once here

            if (cand_prize > best_prize + 1e-12) {
                best_prize = cand_prize;
                best_solution = current_solution;
            }
        }

        temp *= cooling_rate;
        if (temp < 1e-6) break;
    }
    return best_solution;
}


std::pair<std::vector<std::vector<int>>, double> TopAnnealing(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& prizes,
        double budget, int nPath,
        const std::vector<int>& starts, int end) {

    const int iterations = 2000;
    const double initial_temp = 100.0;
    const double cooling_rate = 0.99;
    const int num_threads = 40;

    std::vector<std::vector<int>> best_solution_global;
    double best_prize_global = -std::numeric_limits<double>::infinity();

    #pragma omp parallel num_threads(num_threads)
    {
        thread_local std::mt19937_64 rng{std::random_device{}()};
        std::uniform_real_distribution<double> U01(0.0, 1.0);

        const double greediness = U01(rng);

        auto init_solution = construct_paths_multiS(costs, prizes, starts, end, nPath, budget, greediness);

        auto sol = top_annealing(costs, prizes, budget, starts, end,
                                 std::move(init_solution),
                                 iterations, initial_temp, cooling_rate);

        double prize = multiPathsSumPrize(prizes, sol);

        #pragma omp critical
        {
            if (prize > best_prize_global + 1e-12) {
                best_prize_global = prize;
                best_solution_global = std::move(sol);
            }
        }
    }

    return { best_solution_global, best_prize_global };
}
