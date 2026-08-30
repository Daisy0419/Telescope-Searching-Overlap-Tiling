#pragma once

#include "SubmodularSelection.h"

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

// Validate the common telescope-orienteering inputs used by every solver.
void validatePixelOverlapProblem(
    const std::vector<std::vector<double>>& costs,
    const std::vector<std::vector<int>>& member_pixels,
    const CoverageOracle& oracle,
    double budget,
    int start_node,
    int end_node);

// Sum the edge costs along a path.
double pixelOverlapPathCost(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path);

// Check endpoints, node ranges, and the no-revisit requirement.
bool pixelOverlapIsSimpleSTPath(const std::vector<int>& path,
                               int node_count,
                               int start_node,
                               int end_node);

// Check the complete P1 path constraints, including the padded budget.
bool pixelOverlapIsFeasibleSTPath(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node);

// Return {insertion position, added cost}. Position -1 means unavailable.
std::pair<int, double> pixelOverlapCheapestInsertion(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path,
    int node);

// Cost-improving 2-opt for the verified symmetric graph. Endpoints and the
// selected node set are preserved.
void pixelOverlapTwoOpt(const std::vector<std::vector<double>>& costs,
                        std::vector<int>& path);

// Produce an independent standard C++ random engine for one logical worker.
std::mt19937_64 makePixelOverlapEngine(std::uint64_t base_seed,
                                      std::uint64_t stream_number);

// Insert nodes from an order whenever their residual reward is positive and
// their cheapest insertion respects the budget.
std::vector<int> pixelOverlapPathFromOrder(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& order,
    double budget,
    int start_node,
    int end_node);

// Construct a feasible randomized greedy path using marginal reward divided
// by added path cost.
std::vector<int> pixelOverlapRandomizedGreedyPath(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    double budget,
    int start_node,
    int end_node,
    double greediness,
    std::mt19937_64& engine);

// Coverage-aware equivalent of the original swap_by_ratio move. Remove the
// selected tile with the smallest exclusive-reward/cost-saving ratio, then
// insert the unselected tile with the best marginal-reward/added-cost ratio.
std::vector<int> pixelOverlapSwapByRatio(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node);

// Coverage-aware equivalent of the original change_tile move. Randomly
// replace one selected interior tile with a feasible unselected tile.
std::vector<int> pixelOverlapChangeTile(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node,
    std::mt19937_64& engine);

// Move one contiguous interior segment to another position.
std::vector<int> pixelOverlapShiftSegment(
    const std::vector<int>& path,
    std::mt19937_64& engine);

// Randomly permute a contiguous interior segment.
std::vector<int> pixelOverlapShuffleSegment(
    const std::vector<int>& path,
    std::mt19937_64& engine);

// Apply one random three-edge reconnection while preserving both endpoints.
std::vector<int> pixelOverlapThreeOptOneStep(
    const std::vector<int>& path,
    std::mt19937_64& engine);

// Coverage-aware equivalent of the original repairPath. It removes duplicate
// interior nodes, prunes low-exclusive-reward tiles if over budget, applies
// 2-opt, and fills available budget by marginal-reward/added-cost insertion.
void pixelOverlapRepairPath(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    std::vector<int>& path,
    double budget,
    int start_node,
    int end_node);

// Generate a neighbor with the same active strategy mix as the original
// non-overlap implementation:
//   20% ratio exchange, 20% tile replacement, 10% segment shift,
//   10% segment shuffle, and 40% one-step 3-opt.
// Every move is followed by 2-opt and coverage-aware repair.
std::vector<int> pixelOverlapRandomNeighbor(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node,
    std::mt19937_64& engine);
