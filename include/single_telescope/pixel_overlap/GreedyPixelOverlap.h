#pragma once
#include <vector>

// Pixel-overlap-aware s-t variants. These deliberately use names distinct
// from the legacy scalar-prize greedy functions, so both implementations can
// be linked into the same executable.
std::vector<int> prizeRatioGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_tile,
    int end_tile);

std::vector<int> prizeGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_tile,
    int end_tile);
