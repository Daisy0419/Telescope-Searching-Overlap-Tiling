#pragma once

#include <limits>
#include <string>
#include <vector>

// Complete fixed-start/fixed-end graph passed to the existing solvers.
// The last two nodes are zero-reward dummy start and end nodes.
struct OrienteeringGraph {
    std::vector<std::vector<double>> costs;
    std::vector<std::vector<int>> member_pixels;

    // Convert graph-node indices back to physical M4OPT fields.
    // Dummy start and end nodes have tile index -1 and exposure time 0.
    std::vector<int> tile_of_node;
    std::vector<double> exposure_of_node;

    int start_node = -1;
    int end_node = -1;
    int physical_tile_count = 0;

    double physical_budget = 0.0;
    double padding = 0.0;
    double solver_budget = 0.0;
};

// Small result returned to Python and M4OPT.
struct M4OptResult {
    std::string algorithm;
    std::vector<int> tile_indices;
    std::vector<double> exposure_seconds;
    std::vector<double> start_seconds;

    double coverage = 0.0;
    double duration_seconds = 0.0;
    double runtime_seconds = 0.0;
    double best_bound = std::numeric_limits<double>::quiet_NaN();
    int completed_depth = 0;
    bool has_path = false;
    bool timed_out = false;
};

// Build the P1 graph from pure tile-to-tile slew times, tile footprints,
// fixed per-tile exposure times, and the physical observing budget.
OrienteeringGraph buildP1Graph(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& dwell_seconds,
    double physical_budget_seconds);

// Reduce deterministic variable-exposure P2 to a P1 mode graph. One mode is
// created for every useful exposure threshold of each tile.
OrienteeringGraph buildP2Graph(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& required_exposure_seconds,
    double minimum_exposure_seconds,
    double maximum_exposure_seconds,
    double physical_budget_seconds);

// Solve fixed-exposure, one-visit, no-window P1 with the named algorithm.
M4OptResult solveP1(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& dwell_seconds,
    double physical_budget_seconds,
    const std::string& algorithm,
    int max_depth = 0,
    long long reward_levels = 100,
    double time_quantum = 0.01,
    double time_limit_seconds = 0.0,
    int thread_count = 0,
    double mip_gap = 0.001);

// Solve deterministic variable-exposure, one-visit, no-window P2 using the
// threshold-mode reduction followed by the named P1 algorithm.
M4OptResult solveP2(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& required_exposure_seconds,
    double minimum_exposure_seconds,
    double maximum_exposure_seconds,
    double physical_budget_seconds,
    const std::string& algorithm,
    int max_depth = 0,
    long long reward_levels = 100,
    double time_quantum = 0.01,
    double time_limit_seconds = 0.0,
    int thread_count = 0,
    double mip_gap = 0.001);

// Return every accepted algorithm name.
std::vector<std::string> availableAlgorithms();
