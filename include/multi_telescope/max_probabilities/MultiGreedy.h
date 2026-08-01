#pragma once
#include <vector>

std::pair<std::vector<std::vector<int>>, double> MultiGreedy(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& prizes,
        double budget, int nPath,
        const std::vector<int>& starts, int end);