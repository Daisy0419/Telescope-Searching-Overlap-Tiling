#pragma once

#include <cstdint>
#include <vector>

// Genetic heuristic for one no-revisit, budget-feasible s-t path. 
//Population evaluation and child construction use OpenMP when it is available.
std::vector<int> genetic_optimization_st_pixel_overlap( const std::vector<std::vector<double>>& costs,
                                                const std::vector<double>& pixel_probabilities,
                                                const std::vector<std::vector<int>>& member_pixels,
                                                double budget,
                                                int start_node,
                                                int end_node,
                                                std::vector<std::vector<int>> initial_paths = {},
                                                int population_size = 40,
                                                int generation_count = 1000,
                                                std::uint64_t random_seed = 5489U,
                                                int thread_count = 0);