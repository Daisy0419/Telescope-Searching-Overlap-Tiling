// Genetic s-t search with exact union-of-HEALPixels scoring.
#include "helplers.h"
#include "GeneticPixelOverlap.h"
#include "LocalSearch.h"
#include "PathSplitter.h"
#include "ReadData.h"
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
#include <numeric>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <utility>   
#include <iterator> 
#include <stdexcept>

// Keep all genetic construction/crossover/mutation helpers private to this
// translation unit. The legacy Genetic.cpp defines helpers with several of
// the same names and may be linked into the same executable.
namespace {

inline std::vector<int> construct_path(const std::vector<std::vector<double>>& costs,
               const std::vector<double>& pixel_probs,
               const std::vector<std::vector<int>>& member_pixels,
               int start, int end,
               double budget, double greediness) {
    const double EPS = 1e-12;
    const int n = static_cast<int>(costs.size());

    if (costs[start][end] > budget + EPS) {
        return {};
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
    auto cover_counts = pixelCoverCounts(
        path, member_pixels, pixel_probs.size());

    struct Cand { int v, pos; double dt, h; };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue;  
            auto [pos, dt] = cheapestInsertionOnePath(costs, path, v);
            if (pos == -1) continue; 
            if (pathCost + dt <= budget + EPS) {
                const double denom = std::max(dt, EPS);
                const double gain = marginalPixelPrize(
                    v, pixel_probs, member_pixels, cover_counts);
                if (gain > EPS) {
                    C.push_back({v, pos, dt, gain / denom});
                }
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
        updatePixelCoverCounts(
            mv.v, 1, member_pixels, cover_counts);
        // assert(path.front() == start && path.back() == end);
        // assert(pathCost <= budgetT + 1e-9);
    }
    // assert(path.front() == start && path.back() == end);
    // assert(pathCost <= budgetT + 1e-9);

    return path;
}


// generate a random path from s-t within budget
std::vector<int> generate_random_st_path(const std::vector<std::vector<double>>& costs, int current_city, 
                                        int end_city, double budget) {
    if (costs[current_city][end_city] > budget + 1e-12) {
        return {};
    }
    std::vector<int> visited; 
    visited.push_back(current_city);

    int num_cities = costs.size();

    std::vector<int> unvisited(num_cities);
    std::iota(unvisited.begin(), unvisited.end(), 0);
    unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), current_city), unvisited.end());
    unvisited.erase(std::remove(unvisited.begin(), unvisited.end(), end_city), unvisited.end());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(unvisited.begin(), unvisited.end(), gen);

    for (int next_city : unvisited) {
        // Reserve the actual next_city -> end edge before accepting the
        // insertion. The former current_city -> end check could create an
        // over-budget individual.
        if (costs[current_city][next_city] +
            costs[next_city][end_city] <= budget + 1e-12) {
            visited.push_back(next_city);
            budget -= costs[current_city][next_city];
            current_city = next_city;
        }
    }
    visited.push_back(end_city);
    return visited;
}

// batch generate s-t random paths, each within budget
std::vector<std::vector<int>> random_st_paths(const std::vector<std::vector<double>>& costs, 
                                            int current_city, int end_city,
                                            double budget, int num_paths){
                                            
    std::vector<std::vector<int>> paths(num_paths);

    #pragma omp parallel for schedule(dynamic, 2)
    for (int i = 0; i < num_paths; ++i) {
        paths[i] = generate_random_st_path(costs, current_city, end_city, budget);
    }

    return paths;
}

// batch generate s-t random paths, each within budget
std::vector<std::vector<int>> random_greedy_st_paths(const std::vector<std::vector<double>>& costs,
                                            const std::vector<double>& pixel_probs,
                                            const std::vector<std::vector<int>>& member_pixels,
                                            int current_city, int end_city,
                                            double budget, int num_paths){
                                            
    std::vector<std::vector<int>> paths(num_paths);
    #pragma omp parallel for
    for (int i = 0; i < num_paths; ++i) {
        std::uniform_real_distribution<double> U01(0.0, 1.0);
        double greediness = U01(rng);
        paths[i] = construct_path(
            costs, pixel_probs, member_pixels,
            current_city, end_city, budget, greediness);
        // paths[i] = generate_random_st_path(costs, current_city, end_city, budget);
    }

    return paths;
}


// generate a path/segment of a path from s-t within budget using greedy-in-cost policy
std::vector<int> greedy_st_path(const std::vector<std::vector<double>>& costs, int start_city, 
                                    int end_city, double budget) {
    std::vector<int> visited;
    visited.push_back(start_city);
    int num_cities = costs.size();
    int current_city = start_city;

    while (budget >= costs[current_city][end_city]) {
        double min_dis = std::numeric_limits<double>::infinity();
        int next_city = -1;

        for (int i = 0; i < num_cities; ++i) {
            if (i != end_city && 
                std::find(visited.begin(), visited.end(), i) == visited.end() &&
                costs[current_city][i] < min_dis &&
                (budget > costs[current_city][i] + costs[i][end_city])) {
                next_city = i;
                min_dis = costs[current_city][i];
            }
        }

        if (next_city == -1) {
            break;
        }

        visited.push_back(next_city);
        budget -= min_dis;
        current_city = next_city;

    }
    visited.push_back(end_city);


    return visited;
}

// batch generate paths
std::vector<std::vector<int>> partial_greedy_paths(const std::vector<std::vector<double>>& costs, int current_city, 
                                                    int end_city, double budget, int num_paths) {

    std::vector<std::vector<int>> paths;
    for(int i = 0; i < num_paths; ++i) {
        std::vector<int> path =
            generate_random_st_path(costs, current_city, end_city, budget);
        if(path.size() < 4){
            paths.emplace_back(path);
            continue;
        }
        std::vector<int> greedy_segment = unique_random_ints(1, path.size()-2, 2);
        int segment1 = std::min(greedy_segment[0], greedy_segment[1]);
        int segment2 = std::max(greedy_segment[0], greedy_segment[1]);

        std::vector<int> head;

        head = std::vector<int>(path.begin(), path.begin() + segment1+1);
        std::vector<int> tail(path.begin()+ segment2, path.end());

        double cost = pathCost(costs, head) + pathCost(costs, tail);
        double remaining_budget = budget - cost;

        std::vector<int> middle = greedy_st_path(costs, path[segment1], path[segment2], remaining_budget);

        head.insert(head.end(), middle.begin()+1, middle.end()-1);
        head.insert(head.end(), tail.begin(), tail.end());

        paths.emplace_back(head);
    }
 
    return paths;
}


//randomly draw a segment from parent1, combine it with parent2
std::vector<int> ordered_crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    int size1 = parent1.size();
    std::vector<int> offspring;
    offspring.push_back(parent1[0]); // s

    std::unordered_set<int> used{parent1[0], parent1.back()}; // s, t
    if(size1 >=4){
        auto seg = unique_random_ints(1, size1 - 2, 2);
        int start = std::min(seg[0], seg[1]), end = std::max(seg[0], seg[1]);

        for (int i = start; i <= end; ++i) {
            offspring.push_back(parent1[i]);
            used.insert(parent1[i]);
        }
    } else {
        for (int i = 1; i < size1-1; ++i) {
            offspring.push_back(parent1[i]);
            used.insert(parent1[i]);
        }
    }

    for (int i = 1; i < (int)parent2.size() - 1; ++i) {
        int city = parent2[i];
        if (!used.count(city)) {
            offspring.push_back(city);
            used.insert(city);
        }
    }

    offspring.push_back(parent1.back()); // t
    return offspring;
}

// extern std::mt19937 rng;     
static inline std::pair<int,int> random_segment(int lo, int hi) {
    std::uniform_int_distribution<int> dist(lo, hi);
    int a = dist(rng), b = dist(rng);
    while (b == a) b = dist(rng);
    if (a > b) std::swap(a, b);
    return {a, b};
}


/*  Prize-/Cost-biased OX-ST crossover  */
std::vector<int> ordered_crossover(
        const std::vector<int>& pA,
        const std::vector<int>& pB,
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        int num_samples = 3) {
    const int n = static_cast<int>(pA.size());
    const int s = pA.front();
    const int t = pA.back();

    /* ── 1. find best random segment in pA by ( Σ prize / Σ edge-cost ) ── */
    int best_lo = 1, best_hi = n - 2;
    double best_density = -std::numeric_limits<double>::infinity();

    if (n >= 4) {                                       // need ≥2 interior nodes
        for (int k = 0; k < num_samples; ++k) {
            auto [lo, hi] = random_segment(1, n - 2);

            /* Exact union probability on the candidate segment. */
            const std::vector<int> segment(
                pA.begin() + lo, pA.begin() + hi + 1);
            const double seg_prize = pathPixelPrize(
                pixel_probs, member_pixels, segment);

            /* travel cost of the segment (edges inside lo..hi) */
            double seg_cost = 0.0;
            for (int i = lo; i < hi; ++i)
                seg_cost += costs[pA[i]][pA[i + 1]];

            double density;
            if (seg_cost < 1e-12)                     // zero-length safeguard
                density = std::numeric_limits<double>::infinity();
            else
                density = seg_prize / seg_cost;

            if (density > best_density) {
                best_density = density;
                best_lo = lo; best_hi = hi;
            }
        }
    }

    /* ── 2. assemble offspring ──────────────────────────────────────────── */
    std::vector<int> child;
    child.reserve(n);
    child.push_back(s);

    std::unordered_set<int> used{ s, t };

    for (int i = best_lo; i <= best_hi; ++i) {
        child.push_back(pA[i]);
        used.insert(pA[i]);
    }

    for (int v : pB) {
        if (!used.count(v) && v != t) {
            child.push_back(v);
            used.insert(v);
        }
    }

    child.push_back(t);                 // fixed endpoint

    return child;
}


// maximize prize/delta_cost ratio
int find_next_city(int current_city, int end_city,
                   const std::vector<std::vector<double>>& costs,
                   const std::vector<double>& pixel_probs,
                   const std::vector<std::vector<int>>& member_pixels,
                   const std::vector<int>& cover_counts,
                   const std::unordered_set<int>& unvisited,
                   double budget_remaining_to_insert) {

    int best_city  = -1;
    double best_ratio = -std::numeric_limits<double>::infinity();

    const double base_edge = costs[current_city][end_city];  

    for (int v : unvisited) {
        double delta = costs[current_city][v] 
                     + costs[v][end_city] 
                     - base_edge; 

        if (delta <= budget_remaining_to_insert) {
            const double gain = marginalPixelPrize(
                v, pixel_probs, member_pixels, cover_counts);
            if (gain <= 1e-12) {
                continue;
            }
            double ratio = gain / std::max(delta, 1e-12);  
            if (ratio > best_ratio) {
                best_ratio = ratio;
                best_city  = v;
            }
        }
    }
    return best_city; 
}


std::vector<int> merge_crossover (const std::vector<int>& parent1, const std::vector<int>& parent2,
                                        const std::vector<std::vector<double>>& costs, 
                                        const std::vector<double>& pixel_probs,
                                        const std::vector<std::vector<int>>& member_pixels,
                                        double budget) {
    int s = parent1[0];
    int t = parent1.back();

    std::vector<int> offspring{ s };
    std::unordered_set<int> unvisited(parent1.begin() + 1, parent1.end() - 1);
    unvisited.insert(parent2.begin() + 1, parent2.end() - 1);

    double remaining_budget = budget;
    int current_city = s;
    auto cover_counts = pixelCoverCounts(
        offspring, member_pixels, pixel_probs.size());

    while (!unvisited.empty()) {
        double cost_to_end = costs[current_city][t];

        int next_city = find_next_city(
            current_city, t, costs, pixel_probs, member_pixels,
            cover_counts, unvisited, remaining_budget - cost_to_end);

        if (next_city == -1) {
            break;
        }

        double cost_to_next = costs[current_city][next_city];
        offspring.push_back(next_city);
        remaining_budget -= cost_to_next;
        unvisited.erase(next_city);
        updatePixelCoverCounts(
            next_city, 1, member_pixels, cover_counts);
        current_city = next_city;
    }

    // Add t if budget allows
    double cost_to_t = costs[current_city][t];
    if (cost_to_t <= remaining_budget) {
        offspring.push_back(t);
    }

    return offspring;
}

std::vector<int> er_crossover(const std::vector<int>& pA,
                                 const std::vector<int>& pB)
{
    const int s = pA.front();
    const int t = pA.back();

    /* ---------- collect common nodes ----------------------------------- */
    std::unordered_set<int> inA(pA.begin() + 1, pA.end() - 1);
    std::unordered_set<int> common;
    for (int v : pB)
        if (v != s && v != t && inA.count(v)) common.insert(v);

    /* ---------- edge map: node → list of {nbr, sub-path} --------------- */
    using Segment = std::vector<int>;
    struct Edge { int nbr; Segment seg; };
    std::unordered_map<int, std::vector<Edge>> adj;

    auto add_edges = [&](const std::vector<int>& P)
    {
        int last = 0;                              // index of last “border” node
        for (size_t i = 1; i < P.size(); ++i) {
            int v = P[i];
            bool border = (v == t) || common.count(v);
            if (border) {
                int src = P[last], dst = v;
                Segment seg(P.begin() + last + 1, P.begin() + i + 1);  // excl src
                adj[src].push_back({ dst, seg });
                adj[dst];                     // make sure key exists
                last = static_cast<int>(i);
            }
        }
    };
    add_edges(pA);
    add_edges(pB);

    /* ---------- build offspring ---------------------------------------- */
    std::unordered_set<int> visited{ s };
    std::vector<int> child{ s };
    int current = s;

    while (true) {
        /* remove current from all adjacency lists to update degrees */
        for (auto& [node, vec] : adj)
            vec.erase(std::remove_if(vec.begin(), vec.end(),
                      [&](const Edge& e){ return e.nbr == current; }), vec.end());

        /* pick unvisited neighbour of minimum degree */
        int next = -1;
        const auto& nbrs = adj[current];
        if (!nbrs.empty()) {
            size_t best_deg = std::numeric_limits<size_t>::max();
            std::vector<int> candidates;
            for (size_t k = 0; k < nbrs.size(); ++k) {
                int v = nbrs[k].nbr;
                if (visited.count(v)) continue;
                size_t deg = adj[v].size();
                if (deg < best_deg) { best_deg = deg; candidates = { (int)k }; }
                else if (deg == best_deg) candidates.push_back((int)k);
            }
            if (!candidates.empty()) {
                std::uniform_int_distribution<size_t> pick(0, candidates.size()-1);
                next = nbrs[candidates[pick(rng)]].nbr;
            }
        }
        /* if none left, exit loop */
        if (next == -1) break;

        /* splice stored sub-path current → next */
        auto it = std::find_if(adj[current].begin(), adj[current].end(),
                               [&](const Edge& e){ return e.nbr == next; });
        const Segment& seg = it->seg;     // destination included
        for (int v : seg)
            if (!visited.count(v)) { child.push_back(v); visited.insert(v); }

        current = next;
    }

    child.push_back(t);                   // fixed endpoint
    return child;
}


std::vector<std::vector<int>> cross_over(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                        const std::vector<double>& pixel_probs,
                                        const std::vector<std::vector<int>>& member_pixels,
                                        double budget, int num_cross) {

    if (paths.size() < 2) {
        return std::vector<std::vector<int>>(1, paths[0]);  
    }

    std::vector<std::vector<int>> offsprings(num_cross);

    #pragma omp parallel for schedule(dynamic, 2)
    for (int i = 0; i < num_cross; ++i) {
        std::vector<int> parents = unique_random_ints(0, paths.size() - 1, 2);
        int parent1_idx = parents[0];
        int parent2_idx = parents[1];

        const auto& parent1 = paths[parent1_idx];
        const auto& parent2 = paths[parent2_idx];

        std::vector<int> offspring;
        std::uniform_real_distribution<double> method_selector(0.0, 1.0);
        double method = method_selector(rng);
        
        if(method <= 0.33)        
            offspring = std::move(ordered_crossover(
                parent1, parent2, costs, pixel_probs, member_pixels));
        else if(method > 0.33 && method <= 0.66)
            offspring = std::move(er_crossover(parent1, parent2));      
        else 
            offspring = std::move(merge_crossover(
                parent1, parent2, costs, pixel_probs,
                member_pixels, budget));

        repairPath(
            costs, pixel_probs, member_pixels, offspring, budget);
        two_opt_st(offspring, costs);
        // bool has_improve = two_opt_st(offspring, costs);
        // while(has_improve) {
        //     // has_improve = false;
        //     has_improve = two_opt_st(offspring, costs);
        // }
        
        // offsprings.push_back(offspring);
        offsprings[i] = std::move(offspring);
    }

    return offsprings;
}

//randomly swap num_swap_pairs of tiles
std::vector<int> swap_mutate_st(std::vector<int> path, int num_swap_pairs) {
    if (path.size() < 4 || num_swap_pairs <= 0) {
        return path;
    }
    for(int j = 0; j < num_swap_pairs; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-2, 2);
        std::iter_swap(path.begin()+swap_elements[0], path.begin()+swap_elements[1]);  
    } 
    return path;
}

std::vector<int> insert_mutate_st(std::vector<int> path, int num_insertion) {
    if (path.size() < 4 || num_insertion <= 0) {
        return path;
    }
    for(int j = 0; j < num_insertion; ++j) {
        std::vector<int> swap_elements = unique_random_ints(1, path.size()-2, 2);
        int city = path[swap_elements[0]];
        path.erase(path.begin() + swap_elements[0]);
        path.insert(path.begin() + swap_elements[1], city);
    } 

    return path;
}

std::vector<int> reverse_mutate_st(std::vector<int> path) {
    if (path.size() < 4) {
        return path;
    }
    std::vector<int> range = unique_random_ints(1, path.size()-2, 2);
    int head = std::min(range[0], range[1]);
    int tail = std::max(range[0], range[1]);
    std::reverse(path.begin() + head, path.begin() + tail + 1);

    return path;
}

std::vector<std::vector<int>> mutate_st(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int num_mutation) {
    // std::vector<std::vector<int>> mutate_paths;
    std::vector<std::vector<int>> mutate_paths(num_mutation);
    #pragma omp parallel for schedule(dynamic, 2)
    for(int i = 0; i < num_mutation; ++i) {
        int selected_path = random_int(0, paths.size()-1);
        std::vector<int> path;

        std::uniform_real_distribution<double> method_selector(0.0, 1.0);
        double method = method_selector(rng);
        
        if(method <= 0.33)
            path = std::move(swap_mutate_st(paths[selected_path], 5));
        else if(method > 0.33 && method <= 0.66)
            path = std::move(insert_mutate_st(paths[selected_path], 5));
        else
            path = std::move(reverse_mutate_st(paths[selected_path]));

        repairPath(
            costs, pixel_probs, member_pixels, path, budget);
        two_opt_st(path, costs);
        // mutate_paths.emplace_back(path);
        mutate_paths[i] = std::move(path);
    }
    return mutate_paths;
}



std::vector<std::vector<int>> mutate(const std::vector<std::vector<int>>& paths, const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int num_mutation) {
    // std::vector<std::vector<int>> mutate_paths;
    std::vector<std::vector<int>> mutate_paths(num_mutation);
    const int N = static_cast<int>(paths.size());
    std::vector<double> weights(N);
    for (int i = 0; i < N; ++i) weights[i] = N - i;
    #pragma omp parallel for 
    for(int i = 0; i < num_mutation; ++i) {
        std::discrete_distribution<int> pick(weights.begin(), weights.end());
        int selected_path = pick(rng);
        // int selected_path = random_int(0, paths.size()-1);
        std::vector<int> path;

        std::uniform_real_distribution<double> method_selector(0.0, 1.0);
        double method = method_selector(rng);

        if (method < 0.2) {
            path = std::move(swap_by_ratio(
                paths[selected_path], costs, pixel_probs,
                member_pixels, budget));
        }
        else if (method < 0.4) {
            path = std::move(change_tile(
                paths[selected_path], costs, pixel_probs,
                member_pixels, budget));
        }
        else if (method < 0.5) {
            path = std::move(shift_segment(paths[selected_path]));
        }
        else if (method < 0.6) {
            path = std::move(shuffle_segment(paths[selected_path]));
        }
        else {
            path = std::move(three_opt_one_step(paths[selected_path]));
        }

        // if(i < num_mutation/3)
        //     path = std::move(swap_mutate_st(paths[selected_path], 5));
        // else if(i >= num_mutation/3 && i < 2*num_mutation/3)
        //     path = std::move(insert_mutate_st(paths[selected_path], 5));
        // else
        //     path = std::move(reverse_mutate_st(paths[selected_path]));
        repairPath(
            costs, pixel_probs, member_pixels, path, budget);
        two_opt_st(path, costs);
        // bool has_improve = two_opt_st(path, costs);
        // while(has_improve) {
        //     // has_improve = false;
        //     has_improve = two_opt_st(path, costs);
        // }
    
        
        // mutate_paths.emplace_back(path);
        mutate_paths[i] = std::move(path);
    }
    return mutate_paths;
}



std::pair<double, double> fitness(const std::vector<int>& path, const std::vector<std::vector<double>>& costs, 
                                const std::vector<double>& pixel_probs,
                                const std::vector<std::vector<int>>& member_pixels) {
    return {
        pathPixelPrize(pixel_probs, member_pixels, path),
        pathCost(costs, path)
    };

}


std::vector<std::vector<int>> select_top_paths(const std::vector<std::vector<int>>& paths, 
                                               const std::vector<std::vector<double>>& costs, 
                                               const std::vector<double>& pixel_probs,
                                               const std::vector<std::vector<int>>& member_pixels,
                                               size_t top_n) {
    std::vector<std::tuple<std::vector<int>, double, double>> scored_paths;

    for (const auto& path : paths) {
        auto [total_prize, total_cost] = fitness(
            path, costs, pixel_probs, member_pixels);
        scored_paths.push_back({path, total_prize, total_cost});
    }

    std::sort(scored_paths.begin(), scored_paths.end(),
              [](const auto& a, const auto& b) {
                    if (std::get<1>(a) == std::get<1>(b)) {
                        return std::get<2>(a) < std::get<2>(b); 
                    }
                  return std::get<1>(a) > std::get<1>(b);
              });

    std::vector<std::vector<int>> selected_paths;
    for (size_t i = 0; i < std::min(top_n, scored_paths.size()); ++i) {
        selected_paths.push_back(std::get<0>(scored_paths[i]));
    }

    return selected_paths;
}


std::vector<int> get_best_path(
    const std::vector<std::vector<int>>& paths,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels) {

    double best_prize = -1.0;
    std::vector<int> best_path;
    for (const auto& path : paths) {
        const double current_prize = pathPixelPrize(
            pixel_probs, member_pixels, path);
        if (current_prize > best_prize) {
            best_prize = current_prize;
            best_path = path;
        }
    }
    return best_path;
}

std::vector<int> evolution_st(const std::vector<std::vector<double>>& costs,
                           const std::vector<double>& pixel_probs,
                           const std::vector<std::vector<int>>& member_pixels,
                           int start_city, int end_city, int num_path, double budget, int evolution_itr,
                           std::vector<std::vector<int>>& init_paths) {
    // std::vector<std::vector<int>> paths = std::move(partial_greedy_paths(costs, start_city, budget, num_path));
    int num_random_path = std::max(
        0, num_path - static_cast<int>(init_paths.size()));
    // std::vector<std::vector<int>> paths =  std::move(random_st_paths(costs, start_city, end_city, budget, num_random_path));
    std::vector<std::vector<int>> paths = random_greedy_st_paths(
        costs, pixel_probs, member_pixels,
        start_city, end_city, budget, num_random_path);
    paths.insert(paths.begin(), init_paths.begin(), init_paths.end());
    int num_mutation = 0.85 * num_path, num_cross = num_path - num_mutation;
    int top_n = num_path * 0.9, num_new_gen = num_path - top_n;


    std::vector<int> global_best_path = {};
    double global_best_prize = 0;
    int itr_no_improvement = 0;
    int MAX_NO_IMPROVEMENT = 500;
    for (int i = 0; i < evolution_itr; ++i) {
        //Mutate and cross-over paths
        std::vector<std::vector<int>> mutate_paths = mutate(
            paths, costs, pixel_probs, member_pixels,
            budget, num_mutation);
        std::vector<std::vector<int>> cross_over_paths = cross_over(
            paths, costs, pixel_probs, member_pixels,
            budget, num_cross);
        
        //Combine all paths
        std::vector<std::vector<int>> all_paths = paths;
        all_paths.insert(all_paths.end(), mutate_paths.begin(), mutate_paths.end());
        all_paths.insert(all_paths.end(), cross_over_paths.begin(), cross_over_paths.end());

        //top paths
        paths = std::move(select_top_paths(
            all_paths, costs, pixel_probs, member_pixels, top_n));
        std::vector<std::vector<int>> new_gen_paths = random_st_paths(
            costs, start_city, end_city, budget, num_new_gen);
        // std::vector<std::vector<int>> new_gen_paths = std::move(partial_greedy_paths(costs, start_city, budget, num_new_gen));
        paths.insert(paths.end(), new_gen_paths.begin(), new_gen_paths.end());
        double cur_best_prize = pathPixelPrize(
            pixel_probs, member_pixels, paths[0]);
        
        if(cur_best_prize > global_best_prize) {
            global_best_prize = cur_best_prize;
            global_best_path = paths[0];
            itr_no_improvement = 0;
        }else {
            itr_no_improvement += 1;
        }
        if(itr_no_improvement > MAX_NO_IMPROVEMENT) {
            break;
        }
    }

    std::vector<int> best_path = get_best_path(
        paths, pixel_probs, member_pixels);
    return best_path;
}


}  // namespace

std::vector<int> genetic_optimization_st_pixel_overlap(
                          const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& pixel_probs,
                          const std::vector<std::vector<int>>& member_pixels,
                          double budget, int start_city, int end_city,
                          std::vector<std::vector<int>> init_paths) {
    if (costs.size() < 2 || !std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "genetic optimization: invalid graph size or budget.");
    }
    for (const auto& row : costs) {
        if (row.size() != costs.size()) {
            throw std::invalid_argument(
                "genetic optimization: cost matrix must be square.");
        }
    }
    validatePixelCoverageData(costs.size(), pixel_probs, member_pixels);
    if (start_city < 0 || end_city < 0 ||
        start_city >= static_cast<int>(costs.size()) ||
        end_city >= static_cast<int>(costs.size()) ||
        start_city == end_city) {
        throw std::invalid_argument(
            "genetic optimization: invalid start/end node.");
    }
    if (!member_pixels[start_city].empty() ||
        !member_pixels[end_city].empty()) {
        throw std::invalid_argument(
            "genetic optimization: dummy start/end pixel lists must be empty.");
    }
    if (costs[start_city][end_city] > budget + 1e-12) {
        return {};
    }

    int population = 40, iterations = 1000;
    std::vector<std::vector<int>> valid_init_paths;
    valid_init_paths.reserve(std::min<std::size_t>(
        init_paths.size(), static_cast<std::size_t>(population)));
    for (auto& path : init_paths) {
        if (path.empty()) {
            continue;
        }
        if (path.size() < 2 || path.front() != start_city ||
            path.back() != end_city) {
            throw std::invalid_argument(
                "genetic optimization: an initial path has wrong endpoints.");
        }
        repairPath(
            costs, pixel_probs, member_pixels, path, budget);
        if (pathCost(costs, path) <= budget + 1e-12) {
            valid_init_paths.push_back(std::move(path));
        }
    }
    if (valid_init_paths.size() > static_cast<std::size_t>(population)) {
        valid_init_paths = select_top_paths(
            valid_init_paths, costs, pixel_probs, member_pixels,
            static_cast<std::size_t>(population));
    }

    std::vector<int> best_tour = evolution_st(
        costs, pixel_probs, member_pixels,
        start_city, end_city, population, budget, iterations,
        valid_init_paths);
    if (best_tour.size() < 2 || best_tour.front() != start_city ||
        best_tour.back() != end_city) {
        return {};
    }
    repairPath(
        costs, pixel_probs, member_pixels, best_tour, budget);
    while(two_opt_st(best_tour, costs)) {}
    repairPath(
        costs, pixel_probs, member_pixels, best_tour, budget);
    if (best_tour.front() != start_city || best_tour.back() != end_city ||
        pathCost(costs, best_tour) > budget + 1e-12) {
        return {};
    }
    
    return best_tour;

}
