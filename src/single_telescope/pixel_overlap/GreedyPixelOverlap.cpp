// Pixel-overlap-aware s-t greedy implementation.
#include "GreedyPixelOverlap.h"
#include "PathOperations.h"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

std::vector<int> greedyPathTwoFixedPixelOverlapImpl(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_tile,
    int end_tile,
    bool use_ratio
) {
    const int node_count = static_cast<int>(costs.size());
    if (node_count < 2 || start_tile < 0 || start_tile >= node_count ||
        end_tile < 0 || end_tile >= node_count || start_tile == end_tile) {
        throw std::invalid_argument(
            "pixel greedy: invalid graph or start/end node.");
    }
    if (!std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "pixel greedy: budget must be finite and nonnegative.");
    }
    for (const auto& row : costs) {
        if (static_cast<int>(row.size()) != node_count) {
            throw std::invalid_argument(
                "pixel greedy: cost matrix must be square.");
        }
    }
    validatePixelCoverageData(
        costs.size(), pixel_probs, member_pixels);
    if (!member_pixels[start_tile].empty() ||
        !member_pixels[end_tile].empty()) {
        throw std::invalid_argument(
            "pixel greedy: dummy start/end pixel lists must be empty.");
    }

    constexpr double tolerance = 1e-12;
    if (costs[start_tile][end_tile] > budget + tolerance) {
        return {};
    }

    std::vector<int> path{start_tile, end_tile};
    std::vector<char> visited(node_count, 0);
    visited[start_tile] = 1;
    visited[end_tile] = 1;
    auto cover_counts = pixelCoverCounts(
        path, member_pixels, pixel_probs.size());

    int current = start_tile;
    double total_cost = costs[start_tile][end_tile];

    while (true) {
        int best_tile = -1;
        double best_score = -std::numeric_limits<double>::infinity();
        double best_gain = 0.0;
        double best_delta = 0.0;

        for (int tile = 0; tile < node_count; ++tile) {
            if (visited[tile]) {
                continue;
            }

            const double delta =
                costs[current][tile] + costs[tile][end_tile] -
                costs[current][end_tile];
            if (total_cost + delta > budget + tolerance) {
                continue;
            }

            const double gain = marginalPixelPrize(
                tile, pixel_probs, member_pixels, cover_counts);
            if (gain <= tolerance) {
                continue;
            }

            const double score = use_ratio
                ? (delta <= tolerance
                       ? std::numeric_limits<double>::infinity()
                       : gain / delta)
                : gain;
            const bool better =
                score > best_score + tolerance ||
                (std::abs(score - best_score) <= tolerance &&
                 (gain > best_gain + tolerance ||
                  (std::abs(gain - best_gain) <= tolerance &&
                   (delta < best_delta - tolerance || best_tile < 0))));
            if (better) {
                best_tile = tile;
                best_score = score;
                best_gain = gain;
                best_delta = delta;
            }
        }

        if (best_tile < 0) {
            break;
        }
        path.insert(path.end() - 1, best_tile);
        visited[best_tile] = 1;
        total_cost += best_delta;
        current = best_tile;
        updatePixelCoverCounts(
            best_tile, 1, member_pixels, cover_counts);
    }

    return path;
}

}  // namespace

std::vector<int> prizeRatioGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_tile,
    int end_tile
) {
    return greedyPathTwoFixedPixelOverlapImpl(
        costs, pixel_probs, member_pixels,
        budget, start_tile, end_tile, true);
}

std::vector<int> prizeGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_tile,
    int end_tile
) {
    return greedyPathTwoFixedPixelOverlapImpl(
        costs, pixel_probs, member_pixels,
        budget, start_tile, end_tile, false);
}
