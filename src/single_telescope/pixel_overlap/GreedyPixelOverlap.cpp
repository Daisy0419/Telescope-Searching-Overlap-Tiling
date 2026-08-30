#include "GreedyPixelOverlap.h"
#include "PixelOverlapPath.h"
#include "SubmodularSelection.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Shared implementation for the marginal and marginal/cost variants.
static std::vector<int> greedyPixelOverlapPath(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node,
    bool divide_by_added_cost) {

    CoverageOracle oracle(pixel_probabilities, member_pixels);
    validatePixelOverlapProblem(
        costs, member_pixels, oracle, budget, start_node, end_node);

    if (costs[start_node][end_node] > budget + 1.0e-12) {
        return {};
    }

    const int node_count = static_cast<int>(costs.size());
    std::vector<int> path{start_node, end_node};
    std::vector<char> used(static_cast<std::size_t>(node_count), 0);
    used[start_node] = 1;
    used[end_node] = 1;
    double current_cost = costs[start_node][end_node];

    oracle.reset();
    oracle.commitAll(path);

    while (true) {
        std::vector<char> eligible(static_cast<std::size_t>(node_count), 0);
        std::vector<int> positions(static_cast<std::size_t>(node_count), -1);
        std::vector<double> gains(static_cast<std::size_t>(node_count), 0.0);
        std::vector<double> added_costs(
            static_cast<std::size_t>(node_count), 0.0);
        std::vector<double> scores(
            static_cast<std::size_t>(node_count),
            -std::numeric_limits<double>::infinity());

        for (int node = 0; node < node_count; ++node) {
            if (used[node] != 0) {
                continue;
            }

            const double gain = oracle.marginal(node);
            if (gain <= 1.0e-15) {
                continue;
            }

            const std::pair<int, double> insertion =
                pixelOverlapCheapestInsertion(costs, path, node);
            if (insertion.first < 0 ||
                current_cost + insertion.second > budget + 1.0e-12) {
                continue;
            }

            eligible[node] = 1;
            positions[node] = insertion.first;
            gains[node] = gain;
            added_costs[node] = insertion.second;
            scores[node] = divide_by_added_cost
                ? gain / std::max(insertion.second, 1.0e-12)
                : gain;
        }

        int best_node = -1;
        for (int node = 0; node < node_count; ++node) {
            if (eligible[node] == 0) {
                continue;
            }

            bool better = best_node < 0;
            if (!better && scores[node] > scores[best_node] + 1.0e-12) {
                better = true;
            } else if (!better &&
                       std::abs(scores[node] - scores[best_node]) <= 1.0e-12) {
                if (gains[node] > gains[best_node] + 1.0e-12) {
                    better = true;
                } else if (std::abs(gains[node] - gains[best_node]) <= 1.0e-12 &&
                           added_costs[node] <
                               added_costs[best_node] - 1.0e-12) {
                    better = true;
                }
            }

            if (better) {
                best_node = node;
            }
        }

        if (best_node < 0) {
            break;
        }

        path.insert(path.begin() + positions[best_node], best_node);
        current_cost += added_costs[best_node];
        used[best_node] = 1;
        oracle.commit(best_node);
    }

    // pixelOverlapTwoOpt(costs, path);
    if (!pixelOverlapIsFeasibleSTPath(
            costs, path, budget, start_node, end_node)) {
        return {};
    }
    return path;
}

std::vector<int> prizeRatioGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node) {

    return greedyPixelOverlapPath(
        costs, pixel_probabilities, member_pixels,
        budget, start_node, end_node, true);
}

std::vector<int> prizeGreedyPathTwoFixedPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node) {

    return greedyPixelOverlapPath(
        costs, pixel_probabilities, member_pixels,
        budget, start_node, end_node, false);
}