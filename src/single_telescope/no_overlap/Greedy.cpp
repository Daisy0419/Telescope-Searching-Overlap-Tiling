#include "helplers.h"

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

void sortWithIndices(std::vector<double>& values, std::vector<int>& indices) {
    std::vector<std::pair<double, int>> valueIndexPairs;
    for (size_t i = 0; i < values.size(); ++i) {
        valueIndexPairs.emplace_back(values[i], i);
    }

    std::sort(valueIndexPairs.begin(), valueIndexPairs.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first;
              });

    for (size_t i = 0; i < valueIndexPairs.size(); ++i) {
        values[i] = valueIndexPairs[i].first;
        indices[i] = valueIndexPairs[i].second;
    }
}

// s- any path; greedy on prize
std::vector<int> prizeGreedyPathOneFixed(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   int init_tile, double budget) {

    std::vector<std::pair<double, int>> prize(probability.size());

    for (int i = 0; i < probability.size(); ++i) {
        if (i != init_tile) {
            prize[i] = {probability[i], i};
        } else {
            prize[i] = {0.0, i};  
        }
    }

    std::sort(prize.begin(), prize.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first;
              });

    if (init_tile == -1) {
        init_tile = prize.front().second;
    }

    std::vector<int> path;
    int current_tile = init_tile;
    path.push_back(current_tile);

    for (const auto& [tile_probability, tile] : prize) {

        if (tile == current_tile || tile_probability == 0.0 || costs[current_tile][tile] > budget) {
            continue;
        }

        budget -= costs[current_tile][tile];
        current_tile = tile;
        path.push_back(current_tile);

        if (budget <= 0) break;
    }

    return path;
}


// s- any path; greedy on cheapest node
std::vector<int> costGreedyPathOneFixed(const std::vector<std::vector<double>>& distances, int current_city, double budget) {
    std::vector<int> visited;
    visited.push_back(current_city);
    int num_cities = distances.size();

    while (budget > 0) {
        double min_dis = std::numeric_limits<double>::infinity();
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (std::find(visited.begin(), visited.end(), i) == visited.end() &&
                distances[current_city][i] < min_dis) {
                next_city = i;
                min_dis = distances[current_city][i];
            }
        }

        if(min_dis > budget) 
            break;

        visited.push_back(next_city);
        budget -= min_dis;
        current_city = next_city;
        // std::cout << "budget: " << budget << std::endl;
    }

    return visited;
}

// s- any path; greedy on prize/distance ratio
std::vector<int> prizeRatioGreedyPathOneFixed(const std::vector<std::vector<double>>& costs, 
                                   const std::vector<double>& probability, 
                                   int init_tile, double budget) {
    int n = probability.size();

    if (init_tile == -1) {
        init_tile = std::max_element(probability.begin(), probability.end()) - probability.begin();
    }

    std::vector<int> path;
    std::vector<bool> visited(n, false);

    int current_tile = init_tile;
    path.push_back(current_tile);
    visited[current_tile] = true;

    while (true) {
        double best_score = -1.0;
        int best_tile = -1;

        for (int j = 0; j < n; ++j) {
            if (visited[j] || j == current_tile) continue;

            double cost = costs[current_tile][j];
            if (cost <= budget) {
                double score = probability[j] / (cost + 1E-10);
                if (score > best_score) {
                    best_score = score;
                    best_tile = j;
                }
            }
        }

        if (best_tile == -1) break;

        budget -= costs[current_tile][best_tile];
        current_tile = best_tile;
        visited[current_tile] = true;
        path.push_back(current_tile);

        if (budget <= 0) break;
    }

    return path;
}

// s-t path; greedy on prize/distance ratio
std::vector<int> prizeRatioGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                              const std::vector<double>& probability, 
                                              double budget, int start_tile, int end_tile) {
    int n = probability.size();

    if (start_tile == -1) {
        start_tile = std::max_element(probability.begin(), probability.end()) - probability.begin();
    }

    std::vector<int> path;
    std::vector<bool> visited(n, false);

    int current_tile = start_tile;
    path.push_back(current_tile);
    visited[current_tile] = true;

    while (true) {
        double best_score = -1.0;
        int best_tile = -1;

        for (int j = 0; j < n; ++j) {
            if (visited[j] || j == current_tile || j == end_tile) continue;

            double cost = costs[current_tile][j] + costs[j][end_tile];
            if (cost <= budget) {
                double score = probability[j] / (costs[current_tile][j] + 1E-10);
                if (score > best_score) {
                    best_score = score;
                    best_tile = j;
                }
            }
        }

        if (best_tile == -1) break;

        double travel_cost = costs[current_tile][best_tile];
        if (travel_cost > budget) break;

        budget -= travel_cost;
        current_tile = best_tile;
        visited[current_tile] = true;
        path.push_back(current_tile);
    }

    if (costs[current_tile][end_tile] <= budget) {
        path.push_back(end_tile);
    }

    return path;
}

// s-t path; greedy on prize
std::vector<int> prizeGreedyPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                                         const std::vector<double>& probability, 
                                         double budget, int start_tile, int end_tile) {
    std::vector<std::pair<double, int>> prize(probability.size());
    for (int i = 0; i < probability.size(); ++i) {
        prize[i] = {probability[i], i};
    }

    std::sort(prize.begin(), prize.end(),
              [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
                  return a.first > b.first;
              });

    std::vector<bool> visited(probability.size(), false);
    std::vector<int> path{start_tile};
    visited[start_tile] = true;

    int current_tile = start_tile;

    for (const auto& [tile_probability, tile] : prize) {
        if (tile == current_tile || visited[tile] || tile_probability == 0.0) continue;

        double cost_to_tile = costs[current_tile][tile];
        double cost_to_end = costs[tile][end_tile];

        if (cost_to_tile + cost_to_end > budget) continue;

        budget -= cost_to_tile;
        current_tile = tile;
        path.push_back(current_tile);
        visited[tile] = true;
    }

    // Final check: is end_tile reachable?
    if (costs[current_tile][end_tile] <= budget) {
        path.push_back(end_tile);
    } else {
        std::cerr << "something went wrong in prizeGreedyPathTwoFixed\n";
        return {start_tile, end_tile};
    }

    return path;
}
