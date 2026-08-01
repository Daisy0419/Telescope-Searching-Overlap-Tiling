#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <omp.h>
#include <iostream>

#include "SimulatedAnnealing.h"
#include "helplers.h"
#include "LocalSearch.h"
#include "PathOperations.h"

// std::mt19937 rng(std::random_device{}());

double path_cost(const std::vector<std::vector<double>>& costs, const std::vector<int>& path) {
    double cost_sum = 0.0;
    for (size_t i = 0; i < path.size() - 1; i++)
        cost_sum += costs[path[i]][path[i + 1]];
    return cost_sum;
}

double path_prize(const std::vector<double>& prizes, const std::vector<int>& path) {
    double prize_sum = 0.0;
    for (int tile : path)
        prize_sum += prizes[tile];
    return prize_sum;
}

std::vector<int> initial_solution(const std::vector<std::vector<double>> &costs, 
                                    int budget, int start_tile, int end_tile) {

    if (costs[start_tile][end_tile] > budget)
        return {};
    int N = costs.size();     

    std::vector<int> tiles;
    tiles.reserve(N - 2);   

    for (int i = 0; i < N; ++i)
        if (i != start_tile && i != end_tile)
            tiles.push_back(i);

    std::shuffle(tiles.begin(), tiles.end(), rng);
    
    std::vector<int> solution = {start_tile};
    double spent = 0.0;

    for (int tile : tiles) {
        double cost_to_tile = costs[solution.back()][tile];
        double cost_tile_to_end = costs[tile][end_tile];
        if (spent + cost_to_tile + cost_tile_to_end <= budget) {
            solution.push_back(tile);
            spent += cost_to_tile;
        }
    }
    solution.push_back(end_tile);
    // spent += costs[path[path.size() - 2]][end_tile]; 
    return solution;
}


inline std::vector<int> construct_path(const std::vector<std::vector<double>>& costs,
               const std::vector<double>& prizes,
               int start, int end,
               double budget, double greediness) {
    const double EPS = 1e-12;
    const int n = static_cast<int>(costs.size());

    if (costs[start][end] > budget + EPS) {
        return {start, end};
    }

    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double alpha = (greediness < 0.0 ? U01(rng) : std::clamp(greediness, 0.0, 1.0));

    // Route init: [s, t], cost = c(s,t)
    std::vector<int> path{start, end};
    double pathCost = costs[start][end];

    // Node usage to prevent duplicates; mark s,t as taken
    std::vector<char> used(n, 0);
    used[start] = 1; used[end] = 1;

    struct Cand { int v, pos; double dt, h; };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue;  
            auto [pos, dt] = cheapestInsertionOnePath(costs, path, v);
            if (pos == -1) continue; 
            if (pathCost + dt <= budget + EPS) {
                const double denom = std::max(dt, EPS);
                const double h = prizes[v] / denom; // S/delta_t
                C.push_back({v, pos, dt, h});
            }
        }
    };

    int steps = 0;
    while (steps++ < 1000) {
        std::vector<Cand> C; buildCandidates(C);
        if (C.empty()) break;

        // Value-based RCL: threshold = minH + alpha * (maxH - minH)
        double minH = C.front().h, maxH = C.front().h;
        for (const auto& c : C) {
            if (c.h < minH) minH = c.h;
            if (c.h > maxH) maxH = c.h;
        }
        const double thr = minH + alpha * (maxH - minH);

        std::vector<int> RCL; RCL.reserve(C.size());
        for (int i = 0; i < static_cast<int>(C.size()); ++i)
            if (C[i].h + EPS >= thr) RCL.push_back(i);
        if (RCL.empty()) break;

        std::uniform_int_distribution<int> pick(0, static_cast<int>(RCL.size()) - 1);
        const Cand& mv = C[RCL[pick(rng)]];

        // Apply insertion (keeps route within budget and doesn't touch endpoints)
        path.insert(path.begin() + mv.pos, mv.v);
        pathCost += mv.dt;
        used[mv.v] = 1;
        // assert(path.front() == start && path.back() == end);
        // assert(pathCost <= budgetT + 1e-9);
    }
    // assert(path.front() == start && path.back() == end);
    // assert(pathCost <= budgetT + 1e-9);

    return path;
}


std::vector<int> neighbor(const std::vector<int>& path,
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& prizes,
                          double budget,
                          std::vector<bool>& visited) {
    std::uniform_real_distribution<double> method_selector(0.0, 1.0);
    double method = method_selector(rng);
    // std::cout << "method: " << method << "\n";
    std::vector<int> new_path;
    if (method < 0.2) {
        new_path = std::move(swap_by_ratio(path, costs, prizes, budget, visited));
    }
    else if (method < 0.4) {
        new_path = std::move(change_tile(path, costs, prizes, budget, visited));
    }
    else if (method < 0.5) {
        new_path = std::move(shift_segment(path));
    }
    else if (method < 0.6) {
        new_path = std::move(shuffle_segment(path));
    }
    else {
        new_path = std::move(three_opt_one_step(path));
    }

    two_opt_st(new_path, costs);
    repairPath(costs, prizes, new_path, budget, visited);
    // std::cout << "new_path: " << path_prize(prizes, new_path) << "\n";
    return new_path;
}


std::vector<int> simulated_annealing(const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& prizes,
                                     double budget, int start_tile, int end_tile,
                                     std::vector<int> init_path,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate) {
    if (init_path.empty() || init_path[0] != start_tile || init_path.back() != end_tile) {
        init_path = initial_solution(costs, budget, start_tile, end_tile);
    }

    auto current_path = init_path;
    auto best_path = current_path;
    double best_prize = path_prize(prizes, best_path);

    std::vector<bool> visited(costs.size(), false);
    for (int tile : current_path) {
        visited[tile] = true;
    }

    double temp = initial_temp;
    int no_improvement_itr = 0;
    int MAX_NO_IMPROVEMENT = 300;

    for (int iter = 0; iter < iterations; iter++) {
        auto local_visited = visited;
        auto candidate_path = neighbor(current_path, costs, prizes, budget, local_visited);
        // std::cout << "new_path: " << path_prize(prizes, candidate_path) << "\n";
        if (candidate_path[0] != start_tile || candidate_path.back() != end_tile) {
            // std::cout << "wrong new_path! \n";
            continue;
        }

        double candidate_cost = path_cost(costs, candidate_path);
        if (candidate_cost > budget) {
            continue;
        }

        double candidate_prize = path_prize(prizes, candidate_path);
        double delta = candidate_prize - path_prize(prizes, current_path);

        if (delta >= 0 || exp(delta / temp) > std::uniform_real_distribution<>(0.0, 1.0)(rng)) {
            current_path = candidate_path;
            visited = local_visited;

            if (candidate_prize > best_prize) {
                best_path = candidate_path;
                best_prize = candidate_prize;
                // no_improvement_itr = 0;
                // std::cout << "current_best_prize: " << best_prize << "/n";
            } 
            // else{
            //     no_improvement_itr += 1;
            // }

            // if(no_improvement_itr > MAX_NO_IMPROVEMENT) {
            //     break;
            // }
        }

        temp *= cooling_rate;
        if (temp < 1e-6) break;
    }

    two_opt_st(best_path, costs);
    repairPath(costs, prizes, best_path, budget, visited);

    return best_path;
}

std::vector<int> simulated_annealing_parallel2(const std::vector<std::vector<double>>& costs,
                                     const std::vector<double>& prizes,
                                     double budget, int start_tile, int end_tile,
                                     std::vector<int> init_path,
                                     int iterations,
                                     double initial_temp,
                                     double cooling_rate,
                                     int num_threads) {
    if (init_path.empty() || init_path[0] != start_tile || init_path.back() != end_tile) {
        init_path = initial_solution(costs, budget, start_tile, end_tile);
    }

    auto current_path = init_path;
    auto best_path = current_path;
    double best_prize = pathPrize(prizes, best_path);

    std::vector<bool> visited(costs.size(), false);
    for (int tile : current_path) {
        visited[tile] = true;
    }

    double temp = initial_temp;
    int no_improvement_itr = 0;
    int MAX_NO_IMPROVEMENT = 300;
    std::vector<std::vector<int>> candidate_paths(num_threads);
    std::vector<std::vector<bool>> local_visited(num_threads);
    std::vector<double> candidate_costs(num_threads);
    std::vector<double> candidate_prizes(num_threads);



    for (int iter = 0; iter < iterations; iter++) {
        #pragma omp parallel for
        for (int thread = 0; thread < num_threads; thread++) {
            local_visited[thread] = visited;
            candidate_paths[thread] = neighbor(current_path, costs, prizes, budget, local_visited[thread]);

            if (candidate_paths[thread].front() != start_tile || candidate_paths[thread].back() != end_tile) {
                continue;
            }

            candidate_costs[thread] = pathCost(costs, candidate_paths[thread]);
            if (candidate_costs[thread] > budget) {
                candidate_prizes[thread] = -std::numeric_limits<double>::infinity();
            } else {
                candidate_prizes[thread] = path_prize(prizes, candidate_paths[thread]);
            }
        }

        auto best_prize_itr = std::max_element(candidate_prizes.begin(), candidate_prizes.end());
        auto best_idx =  best_prize_itr - candidate_prizes.begin();
        
        double delta =  *best_prize_itr - path_prize(prizes, current_path);
        

        if (delta >= 0 || exp(delta / temp) > std::uniform_real_distribution<>(0.0, 1.0)(rng)) {
            current_path = candidate_paths[best_idx];
            visited = local_visited[best_idx];

            if (*best_prize_itr  > best_prize) {
                best_path = current_path;
                best_prize = *best_prize_itr;
                no_improvement_itr = 0;
                // std::cout << "current_best_prize: " << best_prize << "/n";
            } else{
                no_improvement_itr += 1;
            }

            if(no_improvement_itr > MAX_NO_IMPROVEMENT) {
                break;
            }
        }

        temp *= cooling_rate;
        if (temp < 1e-6) break;
    }

    bool has_improve = true;
    while(has_improve) {
        has_improve = two_opt_st(best_path, costs);
    }
    repairPath(costs, prizes, best_path, budget, visited);
    return best_path;
}


std::vector<int> simulated_annealing_parallel(const std::vector<std::vector<double>>& costs,
                                              const std::vector<double>& prizes,
                                              double budget, int start_tile, int end_tile,
                                              std::vector<std::vector<int>> init_paths,
                                              int iterations, double initial_temp,
                                              double cooling_rate, int num_threads) {
    std::vector<std::vector<int>> paths(num_threads);
    std::vector<double> prizes_found(num_threads, 0.0);

    std::random_device rd;
    unsigned int base_seed = rd();

    #pragma omp parallel for num_threads(num_threads)
    for (int thread = 0; thread < num_threads; thread++) {
        std::mt19937 thread_rng(base_seed + thread);
        double thread_temp = initial_temp * (0.9 + 0.2 * (thread / double(num_threads)));
        // std::cout << "thread_temp: " << thread_temp << "\n";

        paths[thread] = simulated_annealing(costs, prizes, budget, start_tile, end_tile, 
                                        init_paths[thread], iterations, thread_temp, cooling_rate);
            
        prizes_found[thread] = path_prize(prizes, paths[thread]);
    }

    auto best_it = std::max_element(prizes_found.begin(), prizes_found.end());
    int best_idx = std::distance(prizes_found.begin(), best_it);

    return paths[best_idx];
}

std::vector<int> simulated_annealing_optimization(const std::vector<std::vector<double>>& costs,
                                                 const std::vector<double>& prizes,
                                                 double budget, int start_tile,
                                                 int end_tile, std::vector<int> init_path) {
    int iterations = 2000;
    double initial_temp = 100.0;
    double cooling_rate = 0.99;
    int num_threads = 40;

    std::vector<std::vector<int>> init_paths(num_threads);

    if (init_path.empty()) {
        for (int i = 0; i < num_threads; ++i) {
            init_paths[i] = initial_solution(costs, budget, start_tile, end_tile);
            // static thread_local std::mt19937_64 rng{std::random_device{}()};
            // std::uniform_real_distribution<double> U01(0.0, 1.0);
            // const double greediness = U01(rng);
            // init_paths[i] = construct_path(costs, prizes, start_tile, end_tile, budget, greediness);
        }
    } else {
        for (int i = 0; i < static_cast<int>(num_threads * 0.8); ++i) {
            init_paths[i] = init_path;
        }
        for (int i = static_cast<int>(num_threads * 0.8); i < num_threads; ++i) {
            init_paths[i] = initial_solution(costs, budget, start_tile, end_tile);
        }
    }
    
    return simulated_annealing_parallel(costs, prizes, budget, start_tile, end_tile, init_paths,
                                        iterations, initial_temp, cooling_rate, num_threads);

    // if (init_path.empty()) {
    //     static thread_local std::mt19937_64 rng{std::random_device{}()};
    //     std::uniform_real_distribution<double> U01(0.0, 1.0);
    //     const double greediness = U01(rng);
    //     init_path =construct_path(costs, prizes, start_tile, end_tile, budget, greediness);
    // }

    // return simulated_annealing_parallel2(costs, prizes, budget, start_tile, end_tile, init_path,
    //                                     iterations, initial_temp, cooling_rate, num_threads);
}