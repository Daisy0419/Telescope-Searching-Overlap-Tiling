#pragma once

#include "SubmodularSelection.h"

#include <cstdint>
#include <random>
#include <utility>
#include <vector>

// ant-colony heuristic for one no-revisit path
class AntColonyPixelOverlap {
public:
    AntColonyPixelOverlap(const std::vector<std::vector<double>>& costs,
                        const std::vector<double>& pixel_probabilities,
                        const std::vector<std::vector<int>>& member_pixels,
                        int ant_count = 0,
                        int iteration_count = 500,
                        std::uint64_t random_seed = 5489U,
                        int thread_count = 0);

    std::vector<int> ant_colony_optimization_pixel_overlap(double budget,
                                                        int start_node,
                                                        int end_node);

private:
    std::vector<std::vector<double>> costs_;
    CoverageOracle oracle_;
    std::vector<char> node_has_pixels_;
    std::vector<std::vector<double>> pheromones_;

    int node_count_ = 0;
    int ant_count_ = 0;
    int iteration_count_ = 0;
    int thread_count_ = 0;
    std::uint64_t random_seed_ = 0U;

    double pheromone_influence_ = 1.0;
    double distance_influence_ = 1.0;
    double reward_influence_ = 1.0;
    double evaporation_rate_ = 0.1;
    double deposit_scale_ = 20.0;

    std::pair<int, double> chooseNextNode(int current_node,
                                        int end_node,
                                        const std::vector<char>& used,
                                        const CoverageOracle& coverage,
                                        double remaining_budget,
                                        std::mt19937_64& engine) const;

    std::vector<int> constructAntPath(double budget,
                                        int start_node,
                                        int end_node,
                                        CoverageOracle& coverage,
                                        std::mt19937_64& engine) const;
};
