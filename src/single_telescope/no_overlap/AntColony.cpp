#include "AntColony.h"
#include "helplers.h"
#include "PathSplitter.h"
#include "ReadData.h"
#include "LocalSearch.h"

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

AntColony::AntColony(std::vector<std::vector<double>> costs, std::vector<double> prize_) {
    distances = costs;
    prizes = prize_;

    NUM_CITIES = costs.size();
    NUM_ANTS = 2 * NUM_CITIES;

    initialize_heuristic();
    pheromones.resize(NUM_CITIES, std::vector<double>(NUM_CITIES, 0.5)); // initialize_pheromones
}

void AntColony::initialize_heuristic() {
    heuristic.resize(NUM_CITIES, std::vector<double>(NUM_CITIES));
    for (int i = 0; i < NUM_CITIES; i++) {
        for (int j = 0; j < NUM_CITIES; j++) {
            heuristic[i][j] = std::pow(prizes[j] * 5, GAMMA) * 
                               std::pow(1.0 / distances[i][j], BETA);
        }
    }
}


// Compute next hop for an ant
std::pair<int, double> AntColony::choose_next_city(int current_city, int end_city,
                                                   const std::vector<bool>& visited_set,
                                                   double budget) {
    double sum = 0.0;
    std::vector<std::pair<int, double>> candidates; 

    // Single pass to compute sum and collect candidates
    for (int j = 0; j < NUM_CITIES; j++) {
        if (!visited_set[j] && distances[current_city][j] + distances[j][end_city] <= budget) {
            double prob = std::pow(pheromones[current_city][j], ALPHA) * heuristic[current_city][j];
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
std::vector<int> AntColony::ant_tour(double budget, int start_city, int end_city) {
    std::vector<int> visited;
    std::vector<bool> visited_set(NUM_CITIES, false);
    
    visited.push_back(start_city);
    visited_set[start_city] = true;
    visited_set[end_city] = true;

    int current_city = start_city;
    double remaining_budget = budget;

    while (remaining_budget > 0) {
        auto [next_city, cost] = choose_next_city(current_city, end_city, visited_set, remaining_budget);
        if (next_city == -1) break;

        visited.push_back(next_city);
        visited_set[next_city] = true;

        current_city = next_city;
        remaining_budget -= cost;
    }

    if (distances[current_city][end_city] <= remaining_budget) {
        visited.push_back(end_city);
        visited_set[end_city] = true;
    } else {
        std::cerr << "not enough budget to go end city\n";
        visited.push_back(end_city);  
    }

    return visited;
}


double AntColony::simulate_ants(std::vector<int>& best_tour, double &best_prize, 
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
        two_opt_st(visited[ant], distances);
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

        double tour_cost = 0.0;
        double tour_prize = prizes[visited[ant].back()];
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            tour_cost += distances[visited[ant][i]][visited[ant][i + 1]];
            tour_prize += prizes[visited[ant][i]];

        }

        #pragma omp critical
        if (tour_prize > best_prize) {
            best_prize = tour_prize;
            best_tour = visited[ant];
        }

        double contribution = Q * tour_prize / tour_cost;
        for (size_t i = 0; i < visited[ant].size() - 1; i++) {
            #pragma omp atomic
            pheromones[visited[ant][i]][visited[ant][i + 1]] += contribution;
        }
    }

    return best_prize;
}

std::vector<int> AntColony::ant_colony_optimization(double budget, int start_city, int end_city) {
    std::vector<int> best_tour;
    double best_prize = 0;

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        double prize = simulate_ants(best_tour, best_prize, budget, start_city, end_city);
    }
    return best_tour;

}


double AntColony::simulate_ants_top(std::vector<int>& best_tour, double &best_prize, 
                                double initial_budget, int start_city, int end_city, 
                                const SplitterFn& splitter, int m_routes) {

    std::vector<std::vector<int>> visited(NUM_ANTS);

    auto decode = [&](const std::vector<int>& path)
        -> std::pair<std::vector<std::vector<int>>, double>
    {

        Instance I(distances, prizes, start_city, end_city);
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

        auto [team_path, team_prize] = decode(visited[ant]);
        team_prizes[ant] = team_prize;
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
            if (team_prizes[ant] > best_prize) {
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



std::vector<int> AntColony::ant_colony_optimization_top(double budget, int start_city, int end_city,  
                                                        const SplitterFn& splitter, int m_routes) {
    std::vector<int> best_tour;
    double best_prize = 0;

    for (int iteration = 0; iteration < MAX_ITERATIONS; iteration++) {
        double prize = simulate_ants_top(best_tour, best_prize, budget, start_city, end_city, splitter, m_routes);
    }
    return best_tour;

}