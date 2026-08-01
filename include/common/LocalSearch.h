#pragma once

#include <vector>
#include <cstddef>
#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <stdexcept>

// // std::mt19937 rng;
// inline std::mt19937 rng{ std::random_device{}() };

// Each OpenMP worker needs its own engine; sharing std::mt19937 is a data race.
inline thread_local std::mt19937 rng{ std::random_device{}() };

/*****one-pass local improvement, 2-opt, 3-opt, or-opt*****/
// 2-opt, ensures that the path remains a valid open path from fixed start to fixed end
bool two_opt_st(std::vector<int>& path, const std::vector<std::vector<double>>& costs);

// 3-opt neighbourhood for an open s-t path (2-opt case not included )
bool three_opt_st(std::vector<int>& path, const std::vector<std::vector<double>>& cost);

bool or_opt_improvement(std::vector<int>& path,
                             const std::vector<std::vector<double>>& dist,
                             int maxLen = 5);


/*****Intra-path one-step local search, 2-opt, 3-opt, or-opt*****/
// std::mt19937 rng(std::random_device{}());
std::vector<int> two_opt(const std::vector<int>& path);

std::vector<int> three_opt_one_step(const std::vector<int>& path);

/// One random double-bridge move (duplicated with a case in 3-opt in s-t setting)
std::vector<int> double_bridge(const std::vector<int>& path);

std::vector<int> swap_tiles(const std::vector<int>& path);

std::vector<int> shift_segment(const std::vector<int>& path);

std::vector<int> shuffle_segment(const std::vector<int>& path);

//helpler
void insert_high_value_cheapest(std::vector<int>& path,
                                const std::vector<std::vector<double>>& cost,
                                const std::vector<double>& prize,
                                double  budget,
                                double& current_cost,
                                std::vector<bool>& visited);

//helpler
void insert_high_value_end(std::vector<int>& path,
                                const std::vector<std::vector<double>>& cost,
                                const std::vector<double>& prize,
                                double  budget,
                                double& current_cost,
                                std::vector<bool>& visited);

//randomly remove a tile and insert_high_value_cheapest
std::vector<int> change_tile(const std::vector<int>& path,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                double budget,
                std::vector<bool>& visited);
//overload
std::vector<int> change_tile(const std::vector<int>& path,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                double budget);

std::vector<int> change_tile(
                const std::vector<int>& path,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& pixel_probs,
                const std::vector<std::vector<int>>& member_pixels,
                double budget,
                std::vector<bool>& visited);

std::vector<int> change_tile(
                const std::vector<int>& path,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& pixel_probs,
                const std::vector<std::vector<int>>& member_pixels,
                double budget);


//remove worst prize/Δ‑cost node and insert_high_value_cheapest
std::vector<int> swap_by_ratio(std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes,
              double budget, std::vector<bool>& visited);
//overload
std::vector<int> swap_by_ratio(std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes, double budget);


std::vector<int> swap_by_ratio(
              std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& pixel_probs,
              const std::vector<std::vector<int>>& member_pixels,
              double budget,
              std::vector<bool>& visited);

std::vector<int> swap_by_ratio(
              std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& pixel_probs,
              const std::vector<std::vector<int>>& member_pixels,
              double budget);

/*****Inter-path one-step local search, 2-opt, 3-opt, or-opt*****/

// static double path_cost(const std::vector<int>& path,
//                         const std::vector<std::vector<double>>& costs)
// {
//     double c = 0.0;
//     for (size_t i = 1; i < path.size(); ++i)
//         c += costs[path[i - 1]][path[i]];
//     return c;
// }

// static double solution_cost(const std::vector<std::vector<int>>& paths,
//                             const std::vector<std::vector<double>>& costs)
// {
//     double sum = 0.0;
//     for (const auto& p : paths) sum += path_cost(p, costs);
//     return sum;
// }

bool twoOptAllRoutes(const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<int>>& paths,
                    std::vector<double>& pathCosts,
                    double budgetT);

bool swapNodesBetweenPaths(const std::vector<std::vector<double>>& costs,
                           std::vector<std::vector<int>>& paths,
                           std::vector<double>& pathCosts,
                           double budget);

bool replaceNodesInPaths(const std::vector<std::vector<double>>& costs,
                                   std::vector<std::vector<int>>& paths,
                                   std::vector<double>& pathCosts,
                                   const std::vector<double>& prizes,
                                   double budget);

bool insertNodesToPaths(const std::vector<std::vector<double>>& costs,
                         std::vector<std::vector<int>>& paths,
                         std::vector<double>& pathCosts,
                         const std::vector<double>& prizes,
                         double budget);

bool relocate_k(std::vector<std::vector<int>>& paths,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                int k = 1);

bool swap_1_1(std::vector<std::vector<int>>& paths,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes,
              bool accept_if_cost_reduced = true);

bool cross_exchange(std::vector<std::vector<int>>& paths,
                    const std::vector<std::vector<double>>& costs,
                    const std::vector<double>&  prizes);

bool merge_split_greedy_budget(std::vector<std::vector<int>>& paths,
                                const std::vector<std::vector<double>>& costs,
                                const std::vector<double>& prizes,
                                double budget);

bool ejection_chain(std::vector<std::vector<int>>& paths,
                    const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes,
                    int depth = 1);



