#pragma once

#include <vector>

// Solve an s-t orienteering problem whose reward is the probability of the
// union of covered HEALPixels. A pixel contributes at most once even when it
// belongs to several selected tiles.
//
// Graph contract:
//   * costs is a square all-in edge-cost matrix. It already contains slew,
//     acceleration/deceleration, settle, dwell, and terminal padding;
//   * member_pixels has one compact HEALPixel list per graph node;
//   * pixel_probs is indexed by those compact HEALPixel ids;
//   * start and end are distinct dummy nodes with empty pixel lists;
//   * budget already includes terminal padding. Do not add padding here;
//   * initial_route, when supplied, must be a simple feasible start-to-end
//     path and is used only as a MIP warm start.
//
// mip_gap < 0 leaves Gurobi's default unchanged. time_limit <= 0 leaves the
// time limit unset. An empty return value means there was no incumbent or the
// direct start-to-end path was infeasible.
std::vector<int> gurobiPathCoverPixelOverlapMTZ(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap = -1.0,
    double time_limit = -1.0,
    const std::vector<int>& initial_route = {});

// Equivalent objective and feasible paths, using cumulative-time variables in
// addition to ordering constraints. This is useful as an alternative
// formulation/check; it is not expected to return the identical path when
// several optimal paths exist or either solve stops at a time limit.
std::vector<int> gurobiPathCoverPixelOverlapContinuousTime(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap = -1.0,
    double time_limit = -1.0,
    const std::vector<int>& initial_route = {});