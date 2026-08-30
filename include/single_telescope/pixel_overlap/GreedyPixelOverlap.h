#pragma once

#include <vector>

// Greedy by residual pixel probability divided by added path cost
std::vector<int> prizeRatioGreedyPathTwoFixedPixelOverlap(
                        const std::vector<std::vector<double>>& costs,
                        const std::vector<double>& pixel_probabilities,
                        const std::vector<std::vector<int>>& member_pixels,
                        double budget, int start_node, int end_node);

// Greedy by residual pixel probability
std::vector<int> prizeGreedyPathTwoFixedPixelOverlap(
                        const std::vector<std::vector<double>>& costs,
                        const std::vector<double>& pixel_probabilities,
                        const std::vector<std::vector<int>>& member_pixels,
                        double budget, int start_node, int end_node);