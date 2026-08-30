#include "SimulatedAnnealingPixelOverlap.h"

#include "PixelOverlapOpenMP.h"
#include "PixelOverlapPath.h"
#include "SubmodularSelection.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>

static bool annealingCandidateIsBetter(double first_reward,
                                       double first_cost,
                                       double second_reward,
                                       double second_cost) {
    if (first_reward > second_reward + 1.0e-12) {
        return true;
    }
    return std::abs(first_reward - second_reward) <= 1.0e-12 &&
           first_cost < second_cost - 1.0e-12;
}

static std::vector<int> runParallelAnnealingRestarts(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node,
    std::vector<std::vector<int>> initial_paths,
    int iteration_count,
    double initial_temperature,
    double cooling_rate,
    int restart_count,
    std::uint64_t random_seed,
    int thread_count) {

    CoverageOracle oracle(pixel_probabilities, member_pixels);
    validatePixelOverlapProblem(
        costs, member_pixels, oracle, budget, start_node, end_node);

    if (iteration_count < 0 ||
        !std::isfinite(initial_temperature) || initial_temperature <= 0.0 ||
        !std::isfinite(cooling_rate) ||
        cooling_rate <= 0.0 || cooling_rate > 1.0) {
        throw std::invalid_argument(
            "Simulated annealing: invalid temperature schedule.");
    }
    if (thread_count < 0) {
        throw std::invalid_argument(
            "Simulated annealing: thread count cannot be negative.");
    }
    if (costs[start_node][end_node] > budget + 1.0e-12) {
        return {};
    }
    if (restart_count <= 0) {
        restart_count = 1;
    }

    std::vector<std::vector<int>> valid_initial_paths;
    for (std::vector<int>& path : initial_paths) {
        if (pixelOverlapIsFeasibleSTPath(
                costs, path, budget, start_node, end_node)) {
            pixelOverlapRepairPath(
                costs, oracle, path, budget, start_node, end_node);
            valid_initial_paths.push_back(std::move(path));
        }
    }

    const int worker_count = pixelOverlapWorkerCount(
        thread_count, restart_count);

    // Restarts run independently, including their coverage-oracle state.
    std::vector<CoverageOracle> worker_oracles;
    worker_oracles.reserve(static_cast<std::size_t>(worker_count));
    for (int worker = 0; worker < worker_count; ++worker) {
        worker_oracles.push_back(oracle);
    }

    std::vector<std::vector<int>> restart_best_paths(
        static_cast<std::size_t>(restart_count));
    std::vector<double> restart_best_rewards(
        static_cast<std::size_t>(restart_count),
        -std::numeric_limits<double>::infinity());
    std::vector<double> restart_best_costs(
        static_cast<std::size_t>(restart_count),
        std::numeric_limits<double>::infinity());

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (int restart = 0; restart < restart_count; ++restart) {
        CoverageOracle& worker_oracle =
            worker_oracles[pixelOverlapThreadNumber()];
        std::mt19937_64 engine = makePixelOverlapEngine(
            random_seed, static_cast<std::uint64_t>(restart));
        std::uniform_real_distribution<double> unit_random(0.0, 1.0);

        std::vector<int> current_path;
        if (!valid_initial_paths.empty() && restart % 5 != 4) {
            const std::size_t path_index =
                static_cast<std::size_t>(restart) %
                valid_initial_paths.size();
            current_path = valid_initial_paths[path_index];
        } else {
            current_path = pixelOverlapRandomizedGreedyPath(
                costs, worker_oracle, budget, start_node, end_node,
                unit_random(engine), engine);
            pixelOverlapRepairPath(
                costs, worker_oracle, current_path,
                budget, start_node, end_node);
        }
        if (!pixelOverlapIsFeasibleSTPath(
                costs, current_path, budget, start_node, end_node)) {
            current_path = {start_node, end_node};
        }

        double current_reward = worker_oracle.rewardOf(current_path);
        double current_cost = pixelOverlapPathCost(costs, current_path);
        std::vector<int> best_path = current_path;
        double best_reward = current_reward;
        double best_cost = current_cost;

        // Slightly different starting temperatures improve restart diversity.
        const double temperature_factor =
            0.9 + 0.2 *
                (static_cast<double>(restart) /
                 static_cast<double>(std::max(1, restart_count - 1)));
        double temperature = initial_temperature * temperature_factor;

        for (int iteration = 0; iteration < iteration_count; ++iteration) {
            const std::vector<int> candidate = pixelOverlapRandomNeighbor(
                costs, worker_oracle, current_path, budget,
                start_node, end_node, engine);
            if (candidate.empty()) {
                continue;
            }

            const double candidate_reward = worker_oracle.rewardOf(candidate);
            const double candidate_cost = pixelOverlapPathCost(costs, candidate);
            const double reward_change = candidate_reward - current_reward;
            const bool accept = reward_change >= 0.0 ||
                std::exp(reward_change / temperature) > unit_random(engine);

            if (accept) {
                current_path = candidate;
                current_reward = candidate_reward;
                current_cost = candidate_cost;

                if (annealingCandidateIsBetter(
                        current_reward, current_cost,
                        best_reward, best_cost)) {
                    best_path = current_path;
                    best_reward = current_reward;
                    best_cost = current_cost;
                }
            }

            temperature *= cooling_rate;
            if (temperature < 1.0e-9) {
                break;
            }
        }

        restart_best_paths[restart] = std::move(best_path);
        restart_best_rewards[restart] = best_reward;
        restart_best_costs[restart] = best_cost;
    }

    int best_restart = 0;
    for (int restart = 1; restart < restart_count; ++restart) {
        if (annealingCandidateIsBetter(
                restart_best_rewards[restart], restart_best_costs[restart],
                restart_best_rewards[best_restart],
                restart_best_costs[best_restart])) {
            best_restart = restart;
        }
    }

    if (!pixelOverlapIsFeasibleSTPath(
            costs, restart_best_paths[best_restart], budget,
            start_node, end_node)) {
        return {};
    }
    return restart_best_paths[best_restart];
}

std::vector<int> simulated_annealing_optimization_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node,
    std::vector<int> initial_path,
    int iteration_count,
    double initial_temperature,
    double cooling_rate,
    int restart_count,
    std::uint64_t random_seed,
    int thread_count) {

    std::vector<std::vector<int>> initial_paths;
    if (!initial_path.empty()) {
        initial_paths.push_back(std::move(initial_path));
    }
    return runParallelAnnealingRestarts(
        costs, pixel_probabilities, member_pixels,
        budget, start_node, end_node, std::move(initial_paths),
        iteration_count, initial_temperature, cooling_rate,
        restart_count, random_seed, thread_count);
}

std::vector<int> simulated_annealing_parallel_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node,
    std::vector<std::vector<int>> initial_paths,
    int iteration_count,
    double initial_temperature,
    double cooling_rate,
    int num_threads) {

    if (num_threads < 0) {
        throw std::invalid_argument(
            "Simulated annealing: thread count cannot be negative.");
    }
    int restart_count = std::max(1, num_threads);
    restart_count = std::max(
        restart_count, static_cast<int>(initial_paths.size()));

    return runParallelAnnealingRestarts(
        costs, pixel_probabilities, member_pixels,
        budget, start_node, end_node, std::move(initial_paths),
        iteration_count, initial_temperature, cooling_rate,
        restart_count, 5489U, num_threads);
}
