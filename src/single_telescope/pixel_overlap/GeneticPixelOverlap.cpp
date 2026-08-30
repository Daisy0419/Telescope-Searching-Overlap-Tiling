#include "GeneticPixelOverlap.h"

#include "PixelOverlapOpenMP.h"
#include "PixelOverlapPath.h"
#include "SubmodularSelection.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>

static bool geneticCandidateIsBetter(double first_reward,
                                     double first_cost,
                                     double second_reward,
                                     double second_cost) {
    if (first_reward > second_reward + 1.0e-12) {
        return true;
    }
    return std::abs(first_reward - second_reward) <= 1.0e-12 &&
           first_cost < second_cost - 1.0e-12;
}

static void evaluateGeneticPopulation(
    const std::vector<std::vector<int>>& population,
    const std::vector<std::vector<double>>& costs,
    std::vector<CoverageOracle>& worker_oracles,
    int worker_count,
    double budget,
    int start_node,
    int end_node,
    std::vector<double>& rewards,
    std::vector<double>& path_costs) {

    const int population_size = static_cast<int>(population.size());
    rewards.assign(static_cast<std::size_t>(population_size),
                   -std::numeric_limits<double>::infinity());
    path_costs.assign(static_cast<std::size_t>(population_size),
                      std::numeric_limits<double>::infinity());

#ifndef _OPENMP
    (void)worker_count;
#endif

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (int individual = 0; individual < population_size; ++individual) {
        CoverageOracle& oracle =
            worker_oracles[pixelOverlapThreadNumber()];
        if (pixelOverlapIsFeasibleSTPath(
                costs, population[individual], budget,
                start_node, end_node)) {
            rewards[individual] = oracle.rewardOf(population[individual]);
            path_costs[individual] =
                pixelOverlapPathCost(costs, population[individual]);
        }
    }
}

static std::vector<int> geneticRanking(
    const std::vector<double>& rewards,
    const std::vector<double>& path_costs) {

    std::vector<int> ranking(rewards.size());
    std::iota(ranking.begin(), ranking.end(), 0);
    std::sort(ranking.begin(), ranking.end(),
              [&](int first, int second) {
                  // Exact comparisons keep this a strict weak ordering, as
                  // required by std::sort. Tolerances are used only when
                  // deciding whether a newly found solution is materially
                  // better than the incumbent.
                  if (rewards[first] != rewards[second])
                      return rewards[first] > rewards[second];
                  if (path_costs[first] != path_costs[second])
                      return path_costs[first] < path_costs[second];
                  return first < second;
              });
    return ranking;
}

static int geneticTournament(const std::vector<double>& rewards,
                             const std::vector<double>& path_costs,
                             std::mt19937_64& engine) {
    std::uniform_int_distribution<int> choose(
        0, static_cast<int>(rewards.size()) - 1);
    int winner = choose(engine);
    for (int round = 1; round < 3; ++round) {
        const int challenger = choose(engine);
        if (geneticCandidateIsBetter(
                rewards[challenger], path_costs[challenger],
                rewards[winner], path_costs[winner])) {
            winner = challenger;
        }
    }
    return winner;
}

// Form an order from one random segment of parent A, followed by the unused
// interior nodes of parent B and then the remaining nodes of parent A.
static std::vector<int> geneticCrossoverOrder(
    const std::vector<int>& parent_a,
    const std::vector<int>& parent_b,
    int node_count,
    std::mt19937_64& engine) {

    std::vector<int> order;
    std::vector<char> used(static_cast<std::size_t>(node_count), 0);

    if (parent_a.size() > 2) {
        std::uniform_int_distribution<int> choose(
            1, static_cast<int>(parent_a.size()) - 2);
        int first = choose(engine);
        int last = choose(engine);
        if (first > last) {
            std::swap(first, last);
        }
        for (int position = first; position <= last; ++position) {
            const int node = parent_a[position];
            if (used[node] == 0) {
                order.push_back(node);
                used[node] = 1;
            }
        }
    }

    for (std::size_t position = 1;
         position + 1 < parent_b.size(); ++position) {
        const int node = parent_b[position];
        if (used[node] == 0) {
            order.push_back(node);
            used[node] = 1;
        }
    }
    for (std::size_t position = 1;
         position + 1 < parent_a.size(); ++position) {
        const int node = parent_a[position];
        if (used[node] == 0) {
            order.push_back(node);
            used[node] = 1;
        }
    }

    return order;
}

std::vector<int> genetic_optimization_st_pixel_overlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    int start_node,
    int end_node,
    std::vector<std::vector<int>> initial_paths,
    int population_size,
    int generation_count,
    std::uint64_t random_seed,
    int thread_count) {

    CoverageOracle oracle(pixel_probabilities, member_pixels);
    validatePixelOverlapProblem(
        costs, member_pixels, oracle, budget, start_node, end_node);
    if (population_size < 2 || generation_count <= 0) {
        throw std::invalid_argument(
            "Genetic search: population must be at least two and generations positive.");
    }
    if (thread_count < 0) {
        throw std::invalid_argument(
            "Genetic search: thread count cannot be negative.");
    }
    if (costs[start_node][end_node] > budget + 1.0e-12) {
        return {};
    }

    const int worker_count = pixelOverlapWorkerCount(
        thread_count, population_size);

    // Each worker needs private mutable scratch inside CoverageOracle.
    std::vector<CoverageOracle> worker_oracles;
    worker_oracles.reserve(static_cast<std::size_t>(worker_count));
    for (int worker = 0; worker < worker_count; ++worker) {
        worker_oracles.push_back(oracle);
    }

    std::vector<std::vector<int>> population;
    population.reserve(static_cast<std::size_t>(population_size));
    for (std::vector<int>& path : initial_paths) {
        if (static_cast<int>(population.size()) >= population_size) {
            break;
        }
        if (pixelOverlapIsFeasibleSTPath(
                costs, path, budget, start_node, end_node)) {
            pixelOverlapRepairPath(
                costs, oracle, path, budget, start_node, end_node);
            population.push_back(std::move(path));
        }
    }

    const int generated_count =
        population_size - static_cast<int>(population.size());
    std::vector<std::vector<int>> generated(
        static_cast<std::size_t>(generated_count));

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
    for (int individual = 0; individual < generated_count; ++individual) {
        CoverageOracle& worker_oracle =
            worker_oracles[pixelOverlapThreadNumber()];
        std::mt19937_64 engine = makePixelOverlapEngine(
            random_seed, static_cast<std::uint64_t>(individual));
        std::uniform_real_distribution<double> choose_greediness(0.0, 1.0);
        generated[individual] = pixelOverlapRandomizedGreedyPath(
            costs, worker_oracle, budget, start_node, end_node,
            choose_greediness(engine), engine);
        pixelOverlapRepairPath(
            costs, worker_oracle, generated[individual],
            budget, start_node, end_node);
    }
    population.insert(population.end(),
                      std::make_move_iterator(generated.begin()),
                      std::make_move_iterator(generated.end()));

    std::vector<int> best_path{start_node, end_node};
    double best_reward = oracle.rewardOf(best_path);
    double best_cost = costs[start_node][end_node];
    int generations_without_improvement = 0;
    const int stagnation_limit = std::max(50, generation_count / 4);

    for (int generation = 0; generation < generation_count; ++generation) {
        std::vector<double> rewards;
        std::vector<double> path_costs;
        evaluateGeneticPopulation(
            population, costs, worker_oracles, worker_count,
            budget, start_node, end_node,
            rewards, path_costs);
        const std::vector<int> ranking = geneticRanking(rewards, path_costs);
        const int generation_best = ranking.front();

        if (geneticCandidateIsBetter(
                rewards[generation_best], path_costs[generation_best],
                best_reward, best_cost)) {
            best_path = population[generation_best];
            best_reward = rewards[generation_best];
            best_cost = path_costs[generation_best];
            generations_without_improvement = 0;
        } else {
            ++generations_without_improvement;
        }

        if (generation + 1 == generation_count ||
            generations_without_improvement >= stagnation_limit) {
            break;
        }

        const int elite_count = std::max(1, population_size / 5);
        std::vector<std::vector<int>> next_population(
            static_cast<std::size_t>(population_size));
        for (int elite = 0; elite < elite_count; ++elite) {
            next_population[elite] = population[ranking[elite]];
        }

#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(worker_count)
#endif
        for (int child_index = elite_count;
             child_index < population_size; ++child_index) {
            CoverageOracle& worker_oracle =
                worker_oracles[pixelOverlapThreadNumber()];
            const std::uint64_t stream =
                static_cast<std::uint64_t>(population_size) +
                static_cast<std::uint64_t>(generation) *
                    static_cast<std::uint64_t>(population_size) +
                static_cast<std::uint64_t>(child_index);
            std::mt19937_64 engine =
                makePixelOverlapEngine(random_seed, stream);

            const int first_parent =
                geneticTournament(rewards, path_costs, engine);
            const int second_parent =
                geneticTournament(rewards, path_costs, engine);
            const std::vector<int> order = geneticCrossoverOrder(
                population[first_parent], population[second_parent],
                static_cast<int>(costs.size()), engine);

            std::vector<int> child = pixelOverlapPathFromOrder(
                costs, worker_oracle, order,
                budget, start_node, end_node);

            std::uniform_real_distribution<double> mutation_chance(0.0, 1.0);
            // The original implementation generated approximately 85% of
            // each new population through the rich mutation neighborhood.
            if (mutation_chance(engine) < 0.85) {
                child = pixelOverlapRandomNeighbor(
                    costs, worker_oracle, child, budget,
                    start_node, end_node, engine);
            }
            pixelOverlapRepairPath(
                costs, worker_oracle, child,
                budget, start_node, end_node);

            if (!pixelOverlapIsFeasibleSTPath(
                    costs, child, budget, start_node, end_node)) {
                child = population[first_parent];
            }
            next_population[child_index] = std::move(child);
        }

        population.swap(next_population);
    }

    if (!pixelOverlapIsFeasibleSTPath(
            costs, best_path, budget, start_node, end_node)) {
        return {};
    }
    return best_path;
}
