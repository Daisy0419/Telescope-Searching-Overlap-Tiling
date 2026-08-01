#pragma once
#include <vector>

int UnrootedMinKCover(const std::vector<std::vector<double>>& slew_costs,
                      const std::vector<double>& dwell_costs,
                      double budget,
                      std::vector<std::vector<int>>& best_paths);



int RootedOpenMinKCover(const std::vector<std::vector<double>>& slew_costs,
                        const std::vector<double>& dwell_costs,
                        const std::vector<int>& starts,
                        double budget,
                        std::vector<std::vector<int>>& best_paths);