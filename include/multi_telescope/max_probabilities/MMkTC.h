#pragma once
#include <vector>

double binarySearchBestCover(const std::vector<std::vector<double>>& costs, 
                           const std::vector<double>& prizes,
                           const std::vector<int>& indices, 
                           const std::vector<int>& local_roots,
                           int t,
                           std::vector<std::vector<int>>& best_paths, 
                           double budget);

std::pair<std::vector<std::vector<int>>, double> mmktc_wraper(const std::vector<std::vector<double>>& costs, 
                           const std::vector<double>& prizes,
                           double budget,
                           const std::vector<int>& start_indices,
                           int t, int m);

std::pair<std::vector<std::vector<int>>, double> bestPrizeMMkTC(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m);