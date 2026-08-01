// Ant-colony s-t search with dynamic marginal HEALPixel probability.
#include "AntColonyPixelOverlap.h"
#include "helplers.h"
#include "PathSplitter.h"
#include "ReadData.h"
#include "LocalSearch.h"
#include "PathOperations.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <limits>
#include <random>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <omp.h>
#include <stdexcept>

namespace {

void validateCostMatrix(const std::vector<std::vector<double>>& costs) {
    if (costs.size() < 2) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: graph must contain at least start and end nodes.");
    }
    for (const auto& row : costs) {
        if (row.size() != costs.size()) {
            throw std::invalid_argument(
                "AntColonyPixelOverlap: cost matrix must be square.");
        }
        for (double cost : row) {
            if (!std::isfinite(cost) || cost < 0.0) {
                throw std::invalid_argument(
                    "AntColonyPixelOverlap: edge costs must be finite and nonnegative.");
            }
        }
    }
}

}  // namespace

AntColonyPixelOverlap::AntColonyPixelOverlap(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs_,
    const std::vector<std::vector<int>>& member_pixels_
) : distances(costs),
    pixel_probs(pixel_probs_),
    member_pixels(member_pixels_) {
    validateCostMatrix(distances);
    validatePixelCoverageData(
        distances.size(), pixel_probs, member_pixels);

    NUM_CITIES = static_cast<int>(distances.size());
    NUM_ANTS = 2 * NUM_CITIES;
    standalone_prizes.resize(NUM_CITIES, 0.0);
    for (int tile = 0; tile < NUM_CITIES; ++tile) {
        standalone_prizes[tile] = tilePixelPrize(
            tile, pixel_probs, member_pixels);
    }

    initialize_heuristic();
    pheromones.assign(
        NUM_CITIES, std::vector<double>(NUM_CITIES, 0.5));
}

AntColonyPixelOverlap::AntColonyPixelOverlap(
    std::vector<std::vector<double>> costs,
    std::vector<double> prize_) {
    validateCostMatrix(costs);
    if (prize_.size() != costs.size()) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: one prize is required per graph node.");
    }
    for (double prize : prize_) {
        if (!std::isfinite(prize) || prize < 0.0) {
            throw std::invalid_argument(
                "AntColonyPixelOverlap: prizes must be finite and nonnegative.");
        }
    }
    distances = std::move(costs);
    pixel_probs = std::move(prize_);

    NUM_CITIES = static_cast<int>(distances.size());
    NUM_ANTS = 2 * NUM_CITIES;
    member_pixels.resize(NUM_CITIES);
    standalone_prizes.resize(NUM_CITIES, 0.0);
    for (int tile = 0; tile < NUM_CITIES; ++tile) {
        if (pixel_probs[tile] > 0.0) {
            member_pixels[tile].push_back(tile);
        }
        standalone_prizes[tile] = pixel_probs[tile];
    }

    initialize_heuristic();
    pheromones.assign(
        NUM_CITIES, std::vector<double>(NUM_CITIES, 0.5));
}

void AntColonyPixelOverlap::initialize_heuristic() {
    heuristic.resize(NUM_CITIES, std::vector<double>(NUM_CITIES));
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            if (i == j) {
                heuristic[i][j] = 0.0;
                continue;
            }
            const double safe_cost = std::max(distances[i][j], 1e-12);
            heuristic[i][j] = std::pow(1.0 / safe_cost, BETA);
        }
    }
}


// Compute next hop for an ant
std::pair<int, double> AntColonyPixelOverlap::choose_next_city(
                                                   int current_city, int end_city,
                                                   const std::vector<bool>& visited_set,
                                                   const std::vector<int>& cover_counts,
                                                   double budget) {
    double sum = 0.0;
    std::vector<std::pair<int, double>> candidates; 

    // Single pass to compute sum and collect candidates
    for (int j = 0; j < NUM_CITIES; j++) {
        if (!visited_set[j] && distances[current_city][j] + distances[j][end_city] <= budget) {
            const double gain = marginalPixelPrize(
                j, pixel_probs, member_pixels, cover_counts);
            if (gain <= 1e-15) {
                continue;
            }
            double prob = std::pow(pheromones[current_city][j], ALPHA) *
                          std::pow(5.0 * gain, GAMMA) *
                          heuristic[current_city][j];
            if (prob > 0.0) {
                sum += prob;
                candidates.emplace_back(j, prob);
            }
        }
    }

    if (candidates.empty()) {
        return {-1, 0.0};
    }

    // Roulette wheel selection in one pass
    double rand_prob = random_double(0.0, sum);
    double cumulative = 0.0;
    for (const auto& [city, prob] : candidates) {
        cumulative += prob;
        if (rand_prob <= cumulative) {
            return {city, distances[current_city][city]};
        }
    }

    return {candidates.back().first, distances[current_city][candidates.back().first]};
}

// Compute tour of one ant
std::vector<int> AntColonyPixelOverlap::ant_tour(
    double budget, int start_city, int end_city) {
    if (distances[start_city][end_city] > budget + 1e-12) {
        return {};
    }
    std::vector<int> visited;
    std::vector<bool> visited_set(NUM_CITIES, false);
    
    visited.push_back(start_city);
    visited_set[start_city] = true;
    visited_set[end_city] = true;

    int current_city = start_city;
    double remaining_budget = budget;
    auto cover_counts = pixelCoverCounts(
        visited, member_pixels, pixel_probs.size());

    while (remaining_budget > 0) {
        auto [next_city, cost] = choose_next_city(
            current_city, end_city, visited_set,
            cover_counts, remaining_budget);
        if (next_city == -1) break;

        visited.push_back(next_city);
        visited_set[next_city] = true;
        updatePixelCoverCounts(
            next_city, 1, member_pixels, cover_counts);

        current_city = next_city;
        remaining_budget -= cost;
    }

    if (distances[current_city][end_city] <= remaining_budget) {
        visited.push_back(end_city);
        visited_set[end_city] = true;
    }

    return visited;
}


double AntColonyPixelOverlap::simulate_ants(
                                std::vector<int>& best_tour, double &best_prize,
                                double initial_budget, int start_city, int end_city) {

    std::vector<std::vector<int>> visited(NUM_ANTS);

    // Simulate each ant's tour
    #pragma omp parallel for schedule(dynamic, 5)
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        visited[ant] = std::move(ant_tour(initial_budget, start_city, end_city));
    }

    //fix-cross(conducting one pass 2-opt)
    #pragma omp parallel for schedule(dynamic, 5)
    for (size_t ant = 0; ant < visited.size(); ant++) {
        if (!visited[ant].empty()) {
            two_opt_st(visited[ant], distances);
        }
    }

    // Evaporate pheromones
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
        }
    }

    // Update pheromones
    #pragma omp parallel for schedule(dynamic, 5)
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        if (visited[ant].size() < 2) continue;

        const double tour_cost = pathCost(distances, visited[ant]);
        const double tour_prize = pathPixelPrize(
            pixel_probs, member_pixels, visited[ant]);

        #pragma omp critical
        if (best_tour.empty() || tour_prize > best_prize) {
            best_prize = tour_prize;
            best_tour = visited[ant];
        }

        if (tour_cost <= 1e-12) {
            continue;
        }
        double contribution = Q * tour_prize / tour_cost;
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            #pragma omp atomic
            pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
        }
    }

    return best_prize;
}

std::vector<int>
AntColonyPixelOverlap::ant_colony_optimization_pixel_overlap(
    double budget, int start_city, int end_city) {
    if (!std::isfinite(budget) || budget < 0.0 ||
        start_city < 0 || end_city < 0 ||
        start_city >= NUM_CITIES || end_city >= NUM_CITIES ||
        start_city == end_city) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: invalid budget or start/end node.");
    }
    if (!member_pixels[start_city].empty() ||
        !member_pixels[end_city].empty()) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: dummy start/end pixel lists must be empty.");
    }
    if (distances[start_city][end_city] > budget + 1e-12) {
        return {};
    }

    std::vector<int> best_tour;
    double best_prize = -std::numeric_limits<double>::infinity();

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        simulate_ants(
            best_tour, best_prize, budget, start_city, end_city);
    }
    if (best_tour.empty()) {
        return {start_city, end_city};
    }
    return best_tour;

}


double AntColonyPixelOverlap::simulate_ants_top(
                                std::vector<int>& best_tour, double &best_prize,
                                double initial_budget, int start_city, int end_city, 
                                const PixelOverlapSplitterFn& splitter,
                                int m_routes) {

    std::vector<std::vector<int>> visited(NUM_ANTS);

    auto decode = [&](const std::vector<int>& path)
        -> std::pair<std::vector<std::vector<int>>, double>
    {

        Instance I(distances, standalone_prizes, start_city, end_city);
        return splitter(I, initial_budget, path, m_routes); 
    };

    auto concat = [&](const std::vector<std::vector<int>>& paths) {
        if (paths.empty()) return std::vector<int>{};

        // // Optional: ensure shared endpoints in debug builds
        // #ifndef NDEBUG
        // for (const auto& p : paths) {
        //     assert(p.size() >= 2);
        //     assert(p.front() == paths.front().front());
        //     assert(p.back()  == paths.front().back());
        // }
        // #endif

        const int num_paths = static_cast<int>(paths.size());

        size_t total_size = 0;
        for (const auto& p : paths) total_size += p.size();

        // keep s and t once: sum|p| - 2*num_paths + 2
        const size_t final_size = total_size - 2ull * num_paths + 2ull;

        std::vector<int> result;
        result.reserve(final_size);
        result.push_back(paths.front().front()); // s

        for (const auto& p : paths) {
            if (p.size() > 2) {
                result.insert(result.end(), p.begin() + 1, p.end() - 1);
            }
        }

        result.push_back(paths.front().back()); // t
        return result;
    };


    auto maxCost = [&](const std::vector<std::vector<int>>& paths) {
        double max_cost = 0.0;
        for (const auto& path : paths) {                  // const ref
            double C = 0.0;
            for (size_t i = 0; i + 1 < path.size(); ++i)
                C += distances[path[i]][path[i+1]];
            if (C > max_cost) max_cost = C;
        }
        return max_cost;
    };

    // Simulate each ant's tour
    std::vector<double> team_prizes(NUM_ANTS);
    std::vector<double> team_costs(NUM_ANTS);
    #pragma omp parallel for
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        visited[ant] = std::move(ant_tour(initial_budget, start_city, end_city));
        two_opt_st(visited[ant], distances);

        auto [team_path, splitter_prize] = decode(visited[ant]);
        (void)splitter_prize;
        // Rescore the returned team by the union of all covered pixels. A
        // shared HEALPixel contributes once across all routes.
        team_prizes[ant] = pathsPixelPrize(
            pixel_probs, member_pixels, team_path);
        team_costs[ant] = maxCost(team_path);
        visited[ant] = concat(team_path);
    }

    // Evaporate pheromones
    #pragma omp parallel for
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            pheromones[i][j] *= (1.0 - EVAPORATION_RATE);
        }
    }

    // Update pheromones
    #pragma omp parallel for
    for (int ant = 0; ant < NUM_ANTS; ant++) {
        if (visited[ant].size() < 2) continue;

        // Update best (protect shared state)
        #pragma omp critical
        {
            if (best_tour.empty() || team_prizes[ant] > best_prize) {
                best_prize = team_prizes[ant];
                best_tour = visited[ant];
            }
        }

        double denom = team_costs[ant];
        if (denom <= 0) continue; // avoid div-by-zero
        double contribution = Q * team_prizes[ant] / denom;

        for (size_t i = 0; i + 1 < visited[ant].size(); ++i) {
            #pragma omp atomic
            pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
        }
    }

    return best_prize;
}



std::vector<int>
AntColonyPixelOverlap::ant_colony_optimization_top_pixel_overlap(
    double budget, int start_city, int end_city,
    const PixelOverlapSplitterFn& splitter, int m_routes) {
    if (!splitter || m_routes <= 0) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: top optimization requires a splitter and routes.");
    }
    if (!std::isfinite(budget) || budget < 0.0 ||
        start_city < 0 || end_city < 0 ||
        start_city >= NUM_CITIES || end_city >= NUM_CITIES ||
        start_city == end_city) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: invalid budget or start/end node.");
    }
    if (!member_pixels[start_city].empty() ||
        !member_pixels[end_city].empty()) {
        throw std::invalid_argument(
            "AntColonyPixelOverlap: dummy start/end pixel lists must be empty.");
    }
    if (distances[start_city][end_city] > budget + 1e-12) {
        return {};
    }

    std::vector<int> best_tour;
    double best_prize = -std::numeric_limits<double>::infinity();

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        simulate_ants_top(
            best_tour, best_prize, budget, start_city,
            end_city, splitter, m_routes);
    }
    if (best_tour.empty()) {
        return {start_city, end_city};
    }
    return best_tour;

}
