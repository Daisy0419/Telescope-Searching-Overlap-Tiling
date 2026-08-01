#pragma once
#include <vector>

void sortWithIndices(std::vector<double>& values, std::vector<int>& indices);

std::vector<int> prizeGreedyPathOneFixed(const std::vector<std::vector<double>>& costs,
                                   const std::vector<double>& probability,
                                   int init_tile, double budget);

std::vector<int> prizeRatioGreedyPathOneFixed(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   int init_tile, double budget);

std::vector<int> costGreedyPathOneFixed(const std::vector<std::vector<double>>& distances, int current_city, double budget);

std::vector<int> prizeRatioGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   double budget, int start_tile, int end_tile);

std::vector<int> prizeGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                         const std::vector<double>& probability, 
                                         double budget, int start_tile, int end_tile);