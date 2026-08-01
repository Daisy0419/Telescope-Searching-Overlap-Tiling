#pragma once

#include <vector>

// Every public name has a pixel_overlap suffix so this implementation can be
// linked alongside the legacy scalar-prize simulated-annealing code.
double sa_path_cost_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path);

double sa_path_prize_pixel_overlap(
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<int>& path);

std::vector<int> sa_initial_solution_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    double budget, int start_city, int end_city);

std::vector<int> sa_neighbor_pixel_overlap(
    const std::vector<int>& path,
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    std::vector<bool>& visited);

std::vector<int> simulated_annealing_pixel_overlap(
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& pixel_probs,
                          const std::vector<std::vector<int>>& member_pixels,
                          double budget, int start_city, int end_city,
                          std::vector<int> init_path, int iterations, 
                          double initial_temp, double cooling_rate);

std::vector<int> simulated_annealing_parallel2_pixel_overlap(
                                     const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& pixel_probs,
                                     const std::vector<std::vector<int>>& member_pixels,
                                     double budget, int start_tile, int end_tile,
                                     std::vector<int> init_path,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate,
                                     int num_threads); 

std::vector<int> simulated_annealing_parallel_pixel_overlap(
                                              const std::vector<std::vector<double>>& costs,
                                              const std::vector<double>& pixel_probs,
                                              const std::vector<std::vector<int>>& member_pixels,
                                              double budget, int start_city,int end_city,
                                              std::vector<std::vector<int>> init_path,
                                              int iterations, double initial_temp,
                                              double cooling_rate, int num_threads);

std::vector<int> simulated_annealing_optimization_pixel_overlap(
                                   const std::vector<std::vector<double>>& costs,
                                   const std::vector<double>& pixel_probs,
                                   const std::vector<std::vector<int>>& member_pixels,
                                   double budget, int start_city, int end_city,
                                   std::vector<int> init_path={});
