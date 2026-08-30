#pragma once

#include <vector>

// Controls the RG-QP search. A candidate_limit of zero considers every
// midpoint and is required for the paper's approximation guarantee.
struct RecursiveGreedyParams {
    double time_quantum = 1.0;
    long long reward_levels = 1000;
    int max_depth = 1;
    double time_limit = 0.0;
    int candidate_limit = 0;
    bool iterative_deepening = false;
    bool verbose = false;
};

// Complete result and diagnostics returned by the telescope RG-QP solver.
struct RecursiveGreedyResult {
    bool metric_ok = false;
    int guaranteed_arcs = 0;
    double ratio_bound = 0.0;
    std::vector<int> path;
    int completed_depth = 0;
    long long calls = 0;
    bool timed_out = false;
    long long reward_units = 0;
    double reward = 0.0;
    double cost = 0.0;
    double true_duration = 0.0;
    bool guarantee_holds = false;
};

// Solve the single-telescope, no-revisit, pixel-overlap orienteering problem.
//
// costs             Symmetric metric graph cost matrix.
// pixel_probs       Probability of each compact HEALPixel.
// member_pixels[v]  Compact pixels covered by graph node v.
// start/end         Dummy start and end node indices.
// padded_budget     Observation duration plus endpoint padding.
// padding           The endpoint padding included in padded_budget.
RecursiveGreedyResult recursiveGreedyOrienteering(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double padded_budget,
    double padding,
    const RecursiveGreedyParams& params = RecursiveGreedyParams());
