#pragma once

#include <vector>

double path_cost(const std::vector<std::vector<double>>& costs, const std::vector<int>& path);
double path_prize(const std::vector<double>& prizes, const std::vector<int>& path);
std::vector<int> initial_solution(const std::vector<std::vector<double>> &costs, 
                                    int budget, int start_city);


// std::vector<int> two_opt(const std::vector<int>& path);
// std::vector<int> three_opt(const std::vector<int>& path);
// std::vector<int> swap_cities(const std::vector<int>& path);
// std::vector<int> shift_segment(const std::vector<int>& path);
// std::vector<int> shuffle_segment(const std::vector<int>& path);

// std::vector<int> change_city(const std::vector<int>& path,
//                              const std::vector<std::vector<double>>& costs,
//                              const std::vector<double>& prizes,
//                              double budget,
//                              std::vector<bool>& visited);

std::vector<int> add_city(const std::vector<int>& path,
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double budget,
                          std::vector<bool>& visited);

std::vector<int> neighbor(const std::vector<int>& path,
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double budget,
                          std::vector<bool>& visited);

std::vector<int> simulated_annealing(const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double budget, int start_city, int end_city,
                          std::vector<int> init_path, int iterations, 
                          double initial_temp, double cooling_rate);

std::vector<int> simulated_annealing_parallel2(const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& prizes,
                                     double budget, int start_tile, int end_tile,
                                     std::vector<int> init_path,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate,
                                     int num_threads); 

std::vector<int> simulated_annealing_parallel(const std::vector<std::vector<double>>& costs,
                                              const std::vector<double>& prizes,
                                              double budget, int start_city,int end_city,
                                              std::vector<std::vector<int>> init_path,
                                              int iterations, double initial_temp,
                                              double cooling_rate, int num_threads);

std::vector<int> simulated_annealing_optimization(const std::vector<std::vector<double>>& costs,
                                   const std::vector<double>& prizes,
                                   double budget, int start_city, int end_city,
                                   std::vector<int> init_path={});

