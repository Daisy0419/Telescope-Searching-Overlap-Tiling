#include "AntColonyPixelOverlap.h"

#include "PixelOverlapOpenMP.h"
#include "PixelOverlapPath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

static void validateAntCostMatrix(
    const std::vector<std::vector<double>>& costs) {
    if (costs.size() < 2) {
        throw std::invalid_argument(
            "Ant colony: at least start and end nodes are required.");
    }
    for (const std::vector<double>& row : costs) {
        if (row.size() != costs.size()) {
            throw std::invalid_argument(
                "Ant colony: cost matrix must be square.");
        }
        for (double cost : row) {
            if (!std::isfinite(cost) || cost < 0.0) {
                throw std::invalid_argument(
                    "Ant colony: costs must be finite and nonnegative.");
            }
        }
    }
}

AntColonyPixelOverlap::AntColonyPixelOverlap(
                            const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& pixel_probabilities,
                            const std::vector<std::vector<int>>& member_pixels,
                            int ant_count,
                            int iteration_count,
                            std::uint64_t random_seed,
                            int thread_count)
                            : costs_(costs),
                            oracle_(pixel_probabilities, member_pixels),
                            node_has_pixels_(member_pixels.size(), 0),
                            node_count_(static_cast<int>(costs.size())),
                            ant_count_(ant_count),
                            iteration_count_(iteration_count),
                            thread_count_(thread_count),
                            random_seed_(random_seed) {

    validateAntCostMatrix(costs_);
    if (oracle_.nodeCount() != node_count_) {
        throw std::invalid_argument(
            "Ant colony: one pixel list is required per graph node.");
    }
    if (iteration_count_ <= 0) {
        throw std::invalid_argument(
            "Ant colony: iteration count must be positive.");
    }
    if (thread_count_ < 0) {
        throw std::invalid_argument(
            "Ant colony: thread count cannot be negative.");
    }
    if (ant_count_ <= 0) {
        ant_count_ = 2 * node_count_;
    }

    for (int node = 0; node < node_count_; ++node) {
        node_has_pixels_[node] = member_pixels[node].empty() ? 0 : 1;
    }

    pheromones_.assign(
        static_cast<std::size_t>(node_count_),
        std::vector<double>(static_cast<std::size_t>(node_count_), 0.5));
}

std::pair<int, double> AntColonyPixelOverlap::chooseNextNode(
                                                int current_node,
                                                int end_node,
                                                const std::vector<char>& used,
                                                const CoverageOracle& coverage,
                                                double remaining_budget,
                                                std::mt19937_64& engine) const {

    std::vector<int> candidates;
    std::vector<double> log_weights;
    double largest_log_weight = -std::numeric_limits<double>::infinity();

    for (int node = 0; node < node_count_; ++node) {
        if (used[node] != 0) {
            continue;
        }
        if (costs_[current_node][node] + costs_[node][end_node] >
            remaining_budget + 1.0e-12) {
            continue;
        }

        const double gain = coverage.marginal(node);
        if (gain <= 1.0e-15) {
            continue;
        }

        const double pheromone =
            std::max(pheromones_[current_node][node], 1.0e-15);
        const double travel_cost =
            std::max(costs_[current_node][node], 1.0e-12);
        const double log_weight =
            pheromone_influence_ * std::log(pheromone) +
            reward_influence_ * std::log(gain) -
            distance_influence_ * std::log(travel_cost);

        candidates.push_back(node);
        log_weights.push_back(log_weight);
        largest_log_weight = std::max(largest_log_weight, log_weight);
    }

    if (candidates.empty()) {
        return {-1, 0.0};
    }

    std::vector<double> weights(log_weights.size(), 0.0);
    double total_weight = 0.0;
    for (std::size_t index = 0; index < log_weights.size(); ++index) {
        weights[index] = std::exp(log_weights[index] - largest_log_weight);
        total_weight += weights[index];
    }

    std::uniform_real_distribution<double> choose(0.0, total_weight);
    const double selected_weight = choose(engine);
    double cumulative_weight = 0.0;
    for (std::size_t index = 0; index < candidates.size(); ++index) {
        cumulative_weight += weights[index];
        if (selected_weight <= cumulative_weight) {
            const int selected_node = candidates[index];
            return {selected_node, costs_[current_node][selected_node]};
        }
    }

    const int selected_node = candidates.back();
    return {selected_node, costs_[current_node][selected_node]};
}

std::vector<int> AntColonyPixelOverlap::constructAntPath(
                                                    double budget,
                                                    int start_node,
                                                    int end_node,
                                                    CoverageOracle& coverage,
                                                    std::mt19937_64& engine) const {

    if (costs_[start_node][end_node] > budget + 1.0e-12) {
        return {};
    }

    std::vector<int> path{start_node};
    std::vector<char> used(static_cast<std::size_t>(node_count_), 0);
    used[start_node] = 1;
    used[end_node] = 1;

    coverage.reset();
    coverage.commit(start_node);

    int current_node = start_node;
    double remaining_budget = budget;
    while (true) {
        const std::pair<int, double> next = chooseNextNode(
            current_node, end_node, used, coverage,
            remaining_budget, engine);
        if (next.first < 0) {
            break;
        }

        path.push_back(next.first);
        used[next.first] = 1;
        coverage.commit(next.first);
        remaining_budget -= next.second;
        current_node = next.first;
    }

    if (costs_[current_node][end_node] > remaining_budget + 1.0e-12) {
        return {};
    }
    path.push_back(end_node);
    pixelOverlapTwoOpt(costs_, path);
    return path;
}

std::vector<int> AntColonyPixelOverlap::ant_colony_optimization_pixel_overlap(
                                                            double budget,
                                                            int start_node,
                                                            int end_node) {

    if (!std::isfinite(budget) || budget < 0.0 ||
        start_node < 0 || start_node >= node_count_ ||
        end_node < 0 || end_node >= node_count_ ||
        start_node == end_node) {
        throw std::invalid_argument(
            "Ant colony: invalid budget or start/end node.");
    }
    if (node_has_pixels_[start_node] != 0 ||
        node_has_pixels_[end_node] != 0) {
        throw std::invalid_argument(
            "Ant colony: dummy start/end nodes must cover no pixels.");
    }
    if (costs_[start_node][end_node] > budget + 1.0e-12) {
        return {};
    }

    for (std::vector<double>& row : pheromones_) {
        std::fill(row.begin(), row.end(), 0.5);
    }

    const int worker_count = pixelOverlapWorkerCount(
        thread_count_, ant_count_);

    // CoverageOracle owns mutable coverage counters and scratch storage.
    // Give each OpenMP worker its own copy instead of sharing oracle_.
    std::vector<CoverageOracle> worker_oracles;
    worker_oracles.reserve(static_cast<std::size_t>(worker_count));
    for (int worker = 0; worker < worker_count; ++worker) {
        worker_oracles.push_back(oracle_);
    }

    std::vector<int> best_path{start_node, end_node};
    double best_reward = oracle_.rewardOf(best_path);
    double best_cost = costs_[start_node][end_node];

    for (int iteration = 0; iteration < iteration_count_; ++iteration) {
        std::vector<std::vector<int>> paths(
            static_cast<std::size_t>(ant_count_));
        std::vector<double> rewards(
            static_cast<std::size_t>(ant_count_), 0.0);
        std::vector<double> path_costs(
            static_cast<std::size_t>(ant_count_), 0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
        for (int ant = 0; ant < ant_count_; ++ant) {
            CoverageOracle& coverage =
                worker_oracles[pixelOverlapThreadNumber()];
            const std::uint64_t stream =
                static_cast<std::uint64_t>(iteration) *
                    static_cast<std::uint64_t>(ant_count_) +
                static_cast<std::uint64_t>(ant);
            std::mt19937_64 engine =
                makePixelOverlapEngine(random_seed_, stream);

            paths[ant] = constructAntPath(
                budget, start_node, end_node, coverage, engine);
            if (pixelOverlapIsFeasibleSTPath(
                    costs_, paths[ant], budget, start_node, end_node)) {
                rewards[ant] = coverage.rewardOf(paths[ant]);
                path_costs[ant] = pixelOverlapPathCost(costs_, paths[ant]);
            } else {
                paths[ant].clear();
                rewards[ant] = -std::numeric_limits<double>::infinity();
                path_costs[ant] = std::numeric_limits<double>::infinity();
            }
        }

        for (int ant = 0; ant < ant_count_; ++ant) {
            if (paths[ant].empty()) {
                continue;
            }
            if (rewards[ant] > best_reward + 1.0e-12 ||
                (std::abs(rewards[ant] - best_reward) <= 1.0e-12 &&
                 path_costs[ant] < best_cost - 1.0e-12)) {
                best_path = paths[ant];
                best_reward = rewards[ant];
                best_cost = path_costs[ant];
            }
        }

#ifdef _OPENMP
#pragma omp parallel for collapse(2) schedule(static) num_threads(worker_count)
#endif
        for (int first = 0; first < node_count_; ++first) {
            for (int second = 0; second < node_count_; ++second) {
                pheromones_[first][second] *= (1.0 - evaporation_rate_);
                pheromones_[first][second] =
                    std::max(pheromones_[first][second], 1.0e-15);
            }
        }

        // Pheromone updates are sequential and deterministic.
        for (int ant = 0; ant < ant_count_; ++ant) {
            if (paths[ant].size() < 2 ||
                path_costs[ant] <= 1.0e-12 || rewards[ant] <= 0.0) {
                continue;
            }
            const double deposit =
                deposit_scale_ * rewards[ant] / path_costs[ant];
            for (std::size_t position = 1;
                 position < paths[ant].size(); ++position) {
                const int first = paths[ant][position - 1];
                const int second = paths[ant][position];
                pheromones_[first][second] += deposit;
                pheromones_[second][first] += deposit;
            }
        }
    }

    if (!pixelOverlapIsFeasibleSTPath(
            costs_, best_path, budget, start_node, end_node)) {
        return {};
    }
    return best_path;
}