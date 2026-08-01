#pragma once

#include <utility>
#include <vector>

// Result format used by both independent formulations:
//   first  = zero or one complete s-t path;
//   second = Gurobi's objective bound (NaN only when unavailable).
//
// The outer vector is retained to match the requested/common multi-path
// result form. This single-telescope solver returns at most one inner path.
//
// Shared graph contract:
//   * costs is the square all-in edge-time matrix returned by the graph
//     builder and already includes slew, dwell, settle, and terminal padding;
//   * member_pixels has one compact HEALPixel list per graph node;
//   * pixel_probs is indexed by those compact pixel ids;
//   * start and end are distinct dummy nodes with empty pixel lists;
//   * budget already includes terminal padding;
//   * initial_route, when supplied, is a simple feasible s-t warm start.
//
// mip_gap < 0 and time_limit <= 0 leave the corresponding Gurobi parameter
// unset.
// Uses only MTZ ordering variables/constraints for subtour elimination.
std::pair<std::vector<std::vector<int>>, double>
gurobiPathCoverPixelOverlapMTZ(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap = -1.0,
    double time_limit = -1.0,
    const std::vector<int>& initial_route = {});

// Uses only actual cumulative edge-time propagation for subtour elimination;
// it does not create MTZ/order variables. Because time cannot eliminate an
// all-zero-time cycle, the implementation detects that exceptional graph and
// throws with instructions to use MTZ or positive tile dwell/edge times.
std::pair<std::vector<std::vector<int>>, double>
gurobiPathCoverPixelOverlapContinuousTime(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap = -1.0,
    double time_limit = -1.0,
    const std::vector<int>& initial_route = {});
