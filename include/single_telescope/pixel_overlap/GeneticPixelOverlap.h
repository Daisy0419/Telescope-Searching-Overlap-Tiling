#pragma once

#include <vector>

// Optimize an s-t path on the already padded edge-cost graph. The objective
// is the union probability of member_pixels, so a HEALPixel contributes once
// even when several selected tiles cover it. Its distinct name allows it to
// coexist with the legacy genetic_optimization_st implementation.
std::vector<int> genetic_optimization_st_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_city,
    int end_city,
    std::vector<std::vector<int>> init_paths = {});
