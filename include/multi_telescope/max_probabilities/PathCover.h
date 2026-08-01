#pragma once

#include <vector>

std::pair<std::vector<std::vector<int>>, double> PathCoverWraper(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m);