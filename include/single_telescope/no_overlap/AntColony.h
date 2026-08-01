#pragma once

#include "ReadData.h"

#include <iostream>
#include <vector>
#include <string>
#include <utility>
#include <functional>

using SplitterFn = std::function<
    std::pair<std::vector<std::vector<int>>, double>(
        const Instance&, double /*B*/, const std::vector<int>& /*perm*/, int /*m_routes*/
    )
>;

class AntColony {
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
    std::vector<double> prizes;
    std::vector<std::vector<double>> heuristic;

    void initialize_heuristic();

    // for s-t path
    std::pair<int, double> choose_next_city(int current_city, int end_city, 
                                            const std::vector<bool>& visited_set, double budget);
    std::vector<int> ant_tour(double budget, int start_city, int end_city);
    double simulate_ants(std::vector<int> &best_tour, double &best_prize, 
                        double initial_budget, int start_city, int end_city);

    double simulate_ants_top(std::vector<int>& best_tour, double &best_prize, 
                                double initial_budget, int start_city, int end_city, 
                                const SplitterFn& splitter, int m_routes);

public:
    AntColony(std::vector<std::vector<double>> cost, std::vector<double> prizes);
    std::vector<int> ant_colony_optimization(double budget, int start_city, int end_city);
    std::vector<int> ant_colony_optimization_top(double budget, int start_city, int end_city,  
                                                        const SplitterFn& splitter, int m_routes);
};
