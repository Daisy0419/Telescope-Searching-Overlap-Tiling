#pragma once

std::pair<std::vector<std::vector<int>>, double> TopAnnealing(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& prizes,
        double budget, int nPath,
        const std::vector<int>& starts, int end);