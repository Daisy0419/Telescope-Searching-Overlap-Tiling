#pragma once

#include "ReadData.h"

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <functional>

using PixelOverlapSplitterFn = std::function<
    std::pair<std::vector<std::vector<int>>, double>(
        const Instance&, double /*B*/, const std::vector<int>& /*perm*/, int /*m_routes*/
    )
>;

class AntColonyPixelOverlap {
private:
    int NUM_CITIES;
    int NUM_ANTS = 50;
    int MAX_ITERATIONS = 500;
    double ALPHA = 1.0; // Pheromone influence
    double BETA = 1.0;  // Heuristic influence
    double GAMMA = 1.0;  // Prize influence
    double EVAPORATION_RATE = 0.1;
    double Q = 20.0;   // Pheromone deposit factor

    std::vector<std::vector<double>> distances;
    std::vector<std::vector<double>> pheromones;
    std::vector<double> pixel_probs;
    std::vector<std::vector<int>> member_pixels;
    // Standalone tile-union prizes are retained only for the legacy
    // PathSplitter/Instance interface. The ant objective uses pixel_probs and
    // member_pixels directly.
    std::vector<double> standalone_prizes;
    std::vector<std::vector<double>> heuristic;

    void initialize_heuristic();

    // for s-t path
    std::pair<int, double> choose_next_city(int current_city, int end_city, 
                                            const std::vector<bool>& visited_set,
                                            const std::vector<int>& cover_counts,
                                            double budget);
    std::vector<int> ant_tour(double budget, int start_city, int end_city);
    double simulate_ants(std::vector<int> &best_tour, double &best_prize, 
                        double initial_budget, int start_city, int end_city);

    double simulate_ants_top(std::vector<int>& best_tour, double &best_prize, 
                                double initial_budget, int start_city, int end_city, 
                                const PixelOverlapSplitterFn& splitter,
                                int m_routes);

public:
    // Overlap-aware constructor. Dummy start/end nodes must have empty pixel
    // lists; edge costs already include dwell, settle, and terminal padding.
    AntColonyPixelOverlap(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels);

    // Backward-compatible non-overlap constructor: each node prize is modeled
    // as one private synthetic pixel.
    [[deprecated("Use the pixel_probs/member_pixels constructor for overlapping tiles")]]
    AntColonyPixelOverlap(std::vector<std::vector<double>> cost,
                          std::vector<double> prizes);

    std::vector<int> ant_colony_optimization_pixel_overlap(
        double budget, int start_city, int end_city);

    std::vector<int> ant_colony_optimization_top_pixel_overlap(
        double budget, int start_city, int end_city,
        const PixelOverlapSplitterFn& splitter, int m_routes);
};
