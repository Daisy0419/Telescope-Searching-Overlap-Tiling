#pragma once
#include <vector>


std::pair<std::vector<std::vector<int>>, double> PROptimization(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               int start, int end);

std::pair<std::vector<std::vector<int>>, double> PROptimization_parallel(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               int start, int end);

std::pair<std::vector<std::vector<int>>, double> PROptimization_parallel(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               const std::vector<int>& starts, int end);


