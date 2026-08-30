#include "PixelOverlapPath.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <stdexcept>

void validatePixelOverlapProblem(
    const std::vector<std::vector<double>>& costs,
    const std::vector<std::vector<int>>& member_pixels,
    const CoverageOracle& oracle,
    double budget,
    int start_node,
    int end_node) {

    const int node_count = static_cast<int>(costs.size());
    if (node_count < 2) {
        throw std::invalid_argument(
            "Pixel-overlap problem: at least start and end nodes are required.");
    }
    if (!std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "Pixel-overlap problem: budget must be finite and nonnegative.");
    }
    for (const std::vector<double>& row : costs) {
        if (static_cast<int>(row.size()) != node_count) {
            throw std::invalid_argument(
                "Pixel-overlap problem: cost matrix must be square.");
        }
        for (double cost : row) {
            if (!std::isfinite(cost) || cost < 0.0) {
                throw std::invalid_argument(
                    "Pixel-overlap problem: costs must be finite and nonnegative.");
            }
        }
    }
    if (oracle.nodeCount() != node_count ||
        static_cast<int>(member_pixels.size()) != node_count) {
        throw std::invalid_argument(
            "Pixel-overlap problem: one pixel list is required per graph node.");
    }
    if (start_node < 0 || start_node >= node_count ||
        end_node < 0 || end_node >= node_count ||
        start_node == end_node) {
        throw std::invalid_argument(
            "Pixel-overlap problem: invalid start or end node.");
    }
    if (!member_pixels[start_node].empty() ||
        !member_pixels[end_node].empty()) {
        throw std::invalid_argument(
            "Pixel-overlap problem: dummy start/end pixel lists must be empty.");
    }
}

double pixelOverlapPathCost(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path) {

    double total = 0.0;
    for (std::size_t position = 1; position < path.size(); ++position) {
        total += costs[path[position - 1]][path[position]];
    }
    return total;
}

bool pixelOverlapIsSimpleSTPath(const std::vector<int>& path,
                               int node_count,
                               int start_node,
                               int end_node) {
    if (path.size() < 2 || path.front() != start_node ||
        path.back() != end_node) {
        return false;
    }

    std::vector<char> seen(static_cast<std::size_t>(node_count), 0);
    for (int node : path) {
        if (node < 0 || node >= node_count || seen[node] != 0) {
            return false;
        }
        seen[node] = 1;
    }
    return true;
}

bool pixelOverlapIsFeasibleSTPath(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node) {

    if (!pixelOverlapIsSimpleSTPath(
            path, static_cast<int>(costs.size()), start_node, end_node)) {
        return false;
    }
    return pixelOverlapPathCost(costs, path) <= budget + 1.0e-12;
}

std::pair<int, double> pixelOverlapCheapestInsertion(
    const std::vector<std::vector<double>>& costs,
    const std::vector<int>& path,
    int node) {

    if (path.size() < 2) {
        return {-1, std::numeric_limits<double>::infinity()};
    }

    int best_position = -1;
    double best_added_cost = std::numeric_limits<double>::infinity();
    for (std::size_t position = 1; position < path.size(); ++position) {
        const int previous = path[position - 1];
        const int next = path[position];
        const double added_cost =
            costs[previous][node] + costs[node][next] - costs[previous][next];
        if (added_cost < best_added_cost) {
            best_added_cost = added_cost;
            best_position = static_cast<int>(position);
        }
    }
    return {best_position, best_added_cost};
}

void pixelOverlapTwoOpt(const std::vector<std::vector<double>>& costs,
                        std::vector<int>& path) {
    if (path.size() < 4) {
        return;
    }

    bool improved = true;
    while (improved) {
        improved = false;
        for (std::size_t first = 1; first + 1 < path.size(); ++first) {
            for (std::size_t last = first + 1;
                 last + 1 < path.size(); ++last) {
                const int before = path[first - 1];
                const int first_node = path[first];
                const int last_node = path[last];
                const int after = path[last + 1];

                const double old_boundary_cost =
                    costs[before][first_node] + costs[last_node][after];
                const double new_boundary_cost =
                    costs[before][last_node] + costs[first_node][after];

                if (new_boundary_cost + 1.0e-12 < old_boundary_cost) {
                    std::reverse(path.begin() + static_cast<std::ptrdiff_t>(first),
                                 path.begin() + static_cast<std::ptrdiff_t>(last + 1));
                    improved = true;
                }
            }
        }
    }
}

std::mt19937_64 makePixelOverlapEngine(std::uint64_t base_seed,
                                      std::uint64_t stream_number) {
    const std::uint32_t base_low = static_cast<std::uint32_t>(base_seed);
    const std::uint32_t base_high =
        static_cast<std::uint32_t>(base_seed >> 32U);
    const std::uint32_t stream_low = static_cast<std::uint32_t>(stream_number);
    const std::uint32_t stream_high =
        static_cast<std::uint32_t>(stream_number >> 32U);

    std::seed_seq seeds{base_low, base_high, stream_low, stream_high};
    return std::mt19937_64(seeds);
}

std::vector<int> pixelOverlapPathFromOrder(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& order,
    double budget,
    int start_node,
    int end_node) {

    std::vector<int> path{start_node, end_node};
    double current_cost = costs[start_node][end_node];
    if (current_cost > budget + 1.0e-12) {
        return {};
    }

    std::vector<char> used(costs.size(), 0);
    used[start_node] = 1;
    used[end_node] = 1;

    oracle.reset();
    oracle.commitAll(path);

    for (int node : order) {
        if (node < 0 || node >= static_cast<int>(costs.size()) ||
            used[node] != 0 || oracle.marginal(node) <= 1.0e-15) {
            continue;
        }

        const std::pair<int, double> insertion =
            pixelOverlapCheapestInsertion(costs, path, node);
        if (insertion.first < 0 ||
            current_cost + insertion.second > budget + 1.0e-12) {
            continue;
        }

        path.insert(path.begin() + insertion.first, node);
        current_cost += insertion.second;
        used[node] = 1;
        oracle.commit(node);
    }

    return path;
}

std::vector<int> pixelOverlapRandomizedGreedyPath(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    double budget,
    int start_node,
    int end_node,
    double greediness,
    std::mt19937_64& engine) {

    std::vector<int> path{start_node, end_node};
    double current_cost = costs[start_node][end_node];
    if (current_cost > budget + 1.0e-12) {
        return {};
    }

    const int node_count = static_cast<int>(costs.size());
    std::vector<char> used(static_cast<std::size_t>(node_count), 0);
    used[start_node] = 1;
    used[end_node] = 1;

    oracle.reset();
    oracle.commitAll(path);
    greediness = std::max(0.0, std::min(1.0, greediness));

    while (true) {
        std::vector<int> candidates;
        std::vector<int> positions;
        std::vector<double> added_costs;
        std::vector<double> scores;

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

            candidates.push_back(node);
            positions.push_back(insertion.first);
            added_costs.push_back(insertion.second);
            scores.push_back(gain / std::max(insertion.second, 1.0e-12));
        }

        if (candidates.empty()) {
            break;
        }

        const std::pair<std::vector<double>::iterator,
                        std::vector<double>::iterator> limits =
            std::minmax_element(scores.begin(), scores.end());
        const double threshold =
            *limits.first + greediness * (*limits.second - *limits.first);

        std::vector<int> restricted_candidates;
        for (int index = 0; index < static_cast<int>(scores.size()); ++index) {
            if (scores[index] + 1.0e-12 >= threshold) {
                restricted_candidates.push_back(index);
            }
        }

        std::uniform_int_distribution<int> choose(
            0, static_cast<int>(restricted_candidates.size()) - 1);
        const int selected_index = restricted_candidates[choose(engine)];
        const int selected_node = candidates[selected_index];

        path.insert(path.begin() + positions[selected_index], selected_node);
        current_cost += added_costs[selected_index];
        used[selected_node] = 1;
        oracle.commit(selected_node);
    }

    return path;
}

static std::vector<int> pixelOverlapPathWithoutPosition(
    const std::vector<int>& path,
    std::size_t position) {

    std::vector<int> result = path;
    result.erase(result.begin() + static_cast<std::ptrdiff_t>(position));
    return result;
}

static void pixelOverlapRemoveDuplicateInteriorNodes(
    std::vector<int>& path,
    int node_count,
    int start_node,
    int end_node) {

    if (path.size() < 2 || path.front() != start_node ||
        path.back() != end_node) {
        path = {start_node, end_node};
        return;
    }

    std::vector<char> used(static_cast<std::size_t>(node_count), 0);
    used[start_node] = 1;
    used[end_node] = 1;

    std::vector<int> cleaned;
    cleaned.reserve(path.size());
    cleaned.push_back(start_node);
    for (std::size_t position = 1;
         position + 1 < path.size(); ++position) {
        const int node = path[position];
        if (node >= 0 && node < node_count && used[node] == 0) {
            cleaned.push_back(node);
            used[node] = 1;
        }
    }
    cleaned.push_back(end_node);
    path.swap(cleaned);
}

void pixelOverlapRepairPath(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    std::vector<int>& path,
    double budget,
    int start_node,
    int end_node) {

    const double tolerance = 1.0e-12;
    const int node_count = static_cast<int>(costs.size());
    if (costs[start_node][end_node] > budget + tolerance) {
        path.clear();
        return;
    }

    pixelOverlapRemoveDuplicateInteriorNodes(
        path, node_count, start_node, end_node);
    double current_cost = pixelOverlapPathCost(costs, path);

    // If an order-changing move violates the budget, remove the tile with the
    // least exclusive probability per unit of saved path cost.
    while (current_cost > budget + tolerance && path.size() > 2) {
        const double current_reward = oracle.rewardOf(path);
        std::size_t best_position = path.size();
        double best_ratio = std::numeric_limits<double>::infinity();
        double best_saving = -std::numeric_limits<double>::infinity();
        double best_candidate_cost = std::numeric_limits<double>::infinity();

        for (std::size_t position = 1;
             position + 1 < path.size(); ++position) {
            const std::vector<int> candidate =
                pixelOverlapPathWithoutPosition(path, position);
            const double candidate_cost =
                pixelOverlapPathCost(costs, candidate);
            const double saving = current_cost - candidate_cost;
            const double loss = std::max(
                0.0, current_reward - oracle.rewardOf(candidate));
            const double ratio = saving > tolerance
                ? loss / saving
                : std::numeric_limits<double>::infinity();

            bool better = best_position == path.size() ||
                ratio < best_ratio - tolerance;
            if (!better && std::abs(ratio - best_ratio) <= tolerance) {
                better = saving > best_saving + tolerance;
            }
            if (!better && !std::isfinite(ratio) &&
                !std::isfinite(best_ratio)) {
                // Defensive fallback for a non-metric input with no
                // positive-saving removal.
                better = candidate_cost < best_candidate_cost - tolerance;
            }
            if (better) {
                best_position = position;
                best_ratio = ratio;
                best_saving = saving;
                best_candidate_cost = candidate_cost;
            }
        }

        if (best_position == path.size()) {
            path = {start_node, end_node};
            current_cost = costs[start_node][end_node];
            break;
        }
        path.erase(
            path.begin() + static_cast<std::ptrdiff_t>(best_position));
        current_cost = pixelOverlapPathCost(costs, path);
    }

    // Repeatedly shorten the order and refill the freed budget. An insertion
    // changes the marginal gains, so the oracle is updated after every choice.
    while (true) {
        pixelOverlapTwoOpt(costs, path);
        current_cost = pixelOverlapPathCost(costs, path);

        std::vector<char> used(static_cast<std::size_t>(node_count), 0);
        for (int node : path) {
            used[node] = 1;
        }

        oracle.reset();
        oracle.commitAll(path);
        const std::size_t size_before_insertion = path.size();

        while (true) {
            int best_node = -1;
            int best_position = -1;
            double best_added_cost = 0.0;
            double best_gain = 0.0;
            double best_score = -std::numeric_limits<double>::infinity();

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
                    current_cost + insertion.second > budget + tolerance) {
                    continue;
                }

                const double score = insertion.second <= tolerance
                    ? std::numeric_limits<double>::infinity()
                    : gain / insertion.second;
                bool better = best_node < 0 ||
                    score > best_score + tolerance;
                const bool scores_tied =
                    (std::isinf(score) && std::isinf(best_score)) ||
                    std::abs(score - best_score) <= tolerance;
                if (!better && scores_tied) {
                    better = gain > best_gain + tolerance ||
                        (std::abs(gain - best_gain) <= tolerance &&
                         (best_node < 0 ||
                          insertion.second < best_added_cost - tolerance));
                }
                if (better) {
                    best_node = node;
                    best_position = insertion.first;
                    best_added_cost = insertion.second;
                    best_gain = gain;
                    best_score = score;
                }
            }

            if (best_node < 0) {
                break;
            }
            path.insert(path.begin() + best_position, best_node);
            current_cost += best_added_cost;
            used[best_node] = 1;
            oracle.commit(best_node);
        }

        if (path.size() == size_before_insertion) {
            break;
        }
    }

    pixelOverlapTwoOpt(costs, path);
    if (!pixelOverlapIsFeasibleSTPath(
            costs, path, budget, start_node, end_node)) {
        path = {start_node, end_node};
    }
}

std::vector<int> pixelOverlapSwapByRatio(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node) {

    if (!pixelOverlapIsFeasibleSTPath(
            costs, path, budget, start_node, end_node) ||
        path.size() <= 2) {
        return path;
    }

    const double tolerance = 1.0e-12;
    const double current_cost = pixelOverlapPathCost(costs, path);
    const double current_reward = oracle.rewardOf(path);
    std::size_t remove_position = path.size();
    double smallest_ratio = std::numeric_limits<double>::infinity();
    double largest_saving = -std::numeric_limits<double>::infinity();

    for (std::size_t position = 1;
         position + 1 < path.size(); ++position) {
        const std::vector<int> reduced =
            pixelOverlapPathWithoutPosition(path, position);
        const double saving =
            current_cost - pixelOverlapPathCost(costs, reduced);
        const double loss = std::max(
            0.0, current_reward - oracle.rewardOf(reduced));
        const double ratio = loss / std::max(saving, tolerance);

        if (ratio < smallest_ratio - tolerance ||
            (std::abs(ratio - smallest_ratio) <= tolerance &&
             saving > largest_saving + tolerance)) {
            remove_position = position;
            smallest_ratio = ratio;
            largest_saving = saving;
        }
    }

    if (remove_position == path.size()) {
        return path;
    }

    const int removed_node = path[remove_position];
    std::vector<int> candidate =
        pixelOverlapPathWithoutPosition(path, remove_position);
    double candidate_cost = pixelOverlapPathCost(costs, candidate);

    std::vector<char> used(costs.size(), 0);
    for (int node : candidate) {
        used[node] = 1;
    }
    oracle.reset();
    oracle.commitAll(candidate);

    int best_node = -1;
    int best_position = -1;
    double best_gain = 0.0;
    double best_added_cost = 0.0;
    double best_score = -std::numeric_limits<double>::infinity();
    for (int node = 0; node < static_cast<int>(costs.size()); ++node) {
        if (used[node] != 0 || node == removed_node) {
            continue;
        }
        const double gain = oracle.marginal(node);
        if (gain <= 1.0e-15) {
            continue;
        }
        const std::pair<int, double> insertion =
            pixelOverlapCheapestInsertion(costs, candidate, node);
        if (insertion.first < 0 ||
            candidate_cost + insertion.second > budget + tolerance) {
            continue;
        }

        const double score = gain /
            std::max(insertion.second, tolerance);
        const bool scores_tied =
            (std::isinf(score) && std::isinf(best_score)) ||
            std::abs(score - best_score) <= tolerance;
        if (best_node < 0 || score > best_score + tolerance ||
            (scores_tied &&
             (gain > best_gain + tolerance ||
              (std::abs(gain - best_gain) <= tolerance &&
               (best_node < 0 ||
                insertion.second < best_added_cost - tolerance))))) {
            best_node = node;
            best_position = insertion.first;
            best_gain = gain;
            best_added_cost = insertion.second;
            best_score = score;
        }
    }

    if (best_node < 0) {
        return path;
    }
    candidate.insert(candidate.begin() + best_position, best_node);
    return candidate;
}

std::vector<int> pixelOverlapChangeTile(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node,
    std::mt19937_64& engine) {

    if (!pixelOverlapIsFeasibleSTPath(
            costs, path, budget, start_node, end_node) ||
        path.size() <= 2) {
        return path;
    }

    std::vector<int> removal_positions;
    for (int position = 1;
         position + 1 < static_cast<int>(path.size()); ++position) {
        removal_positions.push_back(position);
    }
    std::shuffle(removal_positions.begin(), removal_positions.end(), engine);

    for (int remove_position : removal_positions) {
        const int removed_node = path[remove_position];
        std::vector<int> candidate = pixelOverlapPathWithoutPosition(
            path, static_cast<std::size_t>(remove_position));
        const double candidate_cost = pixelOverlapPathCost(costs, candidate);

        std::vector<char> used(costs.size(), 0);
        for (int node : candidate) {
            used[node] = 1;
        }
        oracle.reset();
        oracle.commitAll(candidate);

        std::vector<int> nodes;
        std::vector<int> positions;
        std::vector<double> weights;
        for (int node = 0; node < static_cast<int>(costs.size()); ++node) {
            if (used[node] != 0 || node == removed_node) {
                continue;
            }
            const double gain = oracle.marginal(node);
            if (gain <= 1.0e-15) {
                continue;
            }
            const std::pair<int, double> insertion =
                pixelOverlapCheapestInsertion(costs, candidate, node);
            if (insertion.first < 0 ||
                candidate_cost + insertion.second > budget + 1.0e-12) {
                continue;
            }

            nodes.push_back(node);
            positions.push_back(insertion.first);
            weights.push_back(
                gain / (1.0 + std::max(0.0, insertion.second)));
        }

        if (!nodes.empty()) {
            std::discrete_distribution<int> choose(
                weights.begin(), weights.end());
            const int selected = choose(engine);
            candidate.insert(
                candidate.begin() + positions[selected], nodes[selected]);
            return candidate;
        }
    }

    return path;
}

std::vector<int> pixelOverlapShiftSegment(
    const std::vector<int>& path,
    std::mt19937_64& engine) {

    const int interior_count = static_cast<int>(path.size()) - 2;
    if (interior_count < 2) {
        return path;
    }

    std::uniform_int_distribution<int> choose_interior(1, interior_count);
    for (int attempt = 0; attempt < 12; ++attempt) {
        int first = choose_interior(engine);
        int last = choose_interior(engine);
        if (first > last) {
            std::swap(first, last);
        }
        if (first == 1 && last == interior_count) {
            continue;
        }

        const std::vector<int> segment(
            path.begin() + first, path.begin() + last + 1);
        std::vector<int> candidate = path;
        candidate.erase(candidate.begin() + first,
                        candidate.begin() + last + 1);

        std::uniform_int_distribution<int> choose_position(
            1, static_cast<int>(candidate.size()) - 1);
        const int new_position = choose_position(engine);
        if (new_position == first) {
            continue;
        }
        candidate.insert(candidate.begin() + new_position,
                         segment.begin(), segment.end());
        if (candidate != path) {
            return candidate;
        }
    }
    return path;
}

std::vector<int> pixelOverlapShuffleSegment(
    const std::vector<int>& path,
    std::mt19937_64& engine) {

    const int interior_count = static_cast<int>(path.size()) - 2;
    if (interior_count < 2) {
        return path;
    }

    std::uniform_int_distribution<int> choose_interior(1, interior_count);
    for (int attempt = 0; attempt < 12; ++attempt) {
        int first = choose_interior(engine);
        int last = choose_interior(engine);
        if (first > last) {
            std::swap(first, last);
        }
        if (first == last) {
            continue;
        }

        std::vector<int> candidate = path;
        std::shuffle(candidate.begin() + first,
                     candidate.begin() + last + 1, engine);
        if (candidate != path) {
            return candidate;
        }
    }
    return path;
}

std::vector<int> pixelOverlapThreeOptOneStep(
    const std::vector<int>& path,
    std::mt19937_64& engine) {

    const int path_size = static_cast<int>(path.size());
    if (path_size < 6) {
        return pixelOverlapShiftSegment(path, engine);
    }

    for (int attempt = 0; attempt < 12; ++attempt) {
        std::uniform_int_distribution<int> choose_first(0, path_size - 4);
        const int first = choose_first(engine);
        std::uniform_int_distribution<int> choose_second(
            first + 1, path_size - 3);
        const int second = choose_second(engine);
        std::uniform_int_distribution<int> choose_third(
            second + 1, path_size - 2);
        const int third = choose_third(engine);
        std::uniform_int_distribution<int> choose_variant(0, 3);
        const int variant = choose_variant(engine);

        std::vector<int> candidate;
        candidate.reserve(path.size());
        candidate.insert(candidate.end(),
                         path.begin(), path.begin() + first + 1);

        if (variant == 0) {
            candidate.insert(candidate.end(),
                             std::make_reverse_iterator(path.begin() + second + 1),
                             std::make_reverse_iterator(path.begin() + first + 1));
            candidate.insert(candidate.end(),
                             std::make_reverse_iterator(path.begin() + third + 1),
                             std::make_reverse_iterator(path.begin() + second + 1));
        } else if (variant == 1) {
            candidate.insert(candidate.end(),
                             path.begin() + second + 1,
                             path.begin() + third + 1);
            candidate.insert(candidate.end(),
                             path.begin() + first + 1,
                             path.begin() + second + 1);
        } else if (variant == 2) {
            candidate.insert(candidate.end(),
                             std::make_reverse_iterator(path.begin() + third + 1),
                             std::make_reverse_iterator(path.begin() + second + 1));
            candidate.insert(candidate.end(),
                             path.begin() + first + 1,
                             path.begin() + second + 1);
        } else {
            candidate.insert(candidate.end(),
                             path.begin() + second + 1,
                             path.begin() + third + 1);
            candidate.insert(candidate.end(),
                             std::make_reverse_iterator(path.begin() + second + 1),
                             std::make_reverse_iterator(path.begin() + first + 1));
        }

        candidate.insert(candidate.end(),
                         path.begin() + third + 1, path.end());
        if (candidate != path) {
            return candidate;
        }
    }
    return path;
}

std::vector<int> pixelOverlapRandomNeighbor(
    const std::vector<std::vector<double>>& costs,
    CoverageOracle& oracle,
    const std::vector<int>& path,
    double budget,
    int start_node,
    int end_node,
    std::mt19937_64& engine) {

    if (!pixelOverlapIsFeasibleSTPath(
            costs, path, budget, start_node, end_node)) {
        return {};
    }

    std::uniform_real_distribution<double> choose_method(0.0, 1.0);
    const double method = choose_method(engine);
    std::vector<int> candidate;
    if (method < 0.2) {
        candidate = pixelOverlapSwapByRatio(
            costs, oracle, path, budget, start_node, end_node);
    } else if (method < 0.4) {
        candidate = pixelOverlapChangeTile(
            costs, oracle, path, budget,
            start_node, end_node, engine);
    } else if (method < 0.5) {
        candidate = pixelOverlapShiftSegment(path, engine);
    } else if (method < 0.6) {
        candidate = pixelOverlapShuffleSegment(path, engine);
    } else {
        candidate = pixelOverlapThreeOptOneStep(path, engine);
    }

    if (candidate.empty()) {
        candidate = path;
    }
    pixelOverlapTwoOpt(costs, candidate);
    pixelOverlapRepairPath(
        costs, oracle, candidate, budget, start_node, end_node);

    if (!pixelOverlapIsFeasibleSTPath(
            costs, candidate, budget, start_node, end_node)) {
        return path;
    }
    return candidate;
}
