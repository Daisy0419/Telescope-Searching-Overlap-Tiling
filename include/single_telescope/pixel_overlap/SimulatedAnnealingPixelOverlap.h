#pragma once

#include <cstdint>
#include <vector>

// Run independent simulated-annealing restarts in parallel and return the best
// no-revisit, budget-feasible path. A thread count of zero uses the OpenMP
// runtime default; without OpenMP, the same code runs sequentially.
std::vector<int> simulated_annealing_optimization_pixel_overlap(
                                            const std::vector<std::vector<double>>& costs,
                                            const std::vector<double>& pixel_probabilities,
                                            const std::vector<std::vector<int>>& member_pixels,
                                            double budget,
                                            int start_node,
                                            int end_node,
                                            std::vector<int> initial_path = {},
                                            int iteration_count = 2000,
                                            double initial_temperature = 100.0,
                                            double cooling_rate = 0.99,
                                            int restart_count = 8,
                                            std::uint64_t random_seed = 5489U,
                                            int thread_count = 0);


std::vector<int> simulated_annealing_parallel_pixel_overlap(
                                            const std::vector<std::vector<double>>& costs,
                                            const std::vector<double>& pixel_probabilities,
                                            const std::vector<std::vector<int>>& member_pixels,
                                            double budget,
                                            int start_node,
                                            int end_node,
                                            std::vector<std::vector<int>> initial_paths,
                                            int iteration_count,
                                            double initial_temperature,
                                            double cooling_rate,
                                            int num_threads);
