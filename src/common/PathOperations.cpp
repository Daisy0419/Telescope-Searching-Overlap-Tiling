#include "PathOperations.h"

#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

/*****Functions for Single OP Path Operations ******/
double pathCost(const std::vector<std::vector<double>>& costs, const std::vector<int>& path) {
    double cost_sum = 0.0;
    for (size_t i = 1; i < path.size(); ++i)
        cost_sum += costs[path[i - 1]][path[i]];
    return cost_sum;
}

double pathPrize(const std::vector<double>& prizes, const std::vector<int>& path) {
    double prize_sum = 0.0;
    for (int tile : path)
        prize_sum += prizes[tile];
    return prize_sum;
}

void validatePixelCoverageData(
    std::size_t node_count,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels ) {
    if (pixel_probs.empty()) {
        throw std::invalid_argument(
            "pixel_probs must contain at least one compact HEALPixel.");
    }
    if (member_pixels.size() != node_count) {
        throw std::invalid_argument(
            "member_pixels must contain one list per graph node.");
    }

    for (double probability : pixel_probs) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "pixel probabilities must be finite and nonnegative.");
        }
    }

    for (const auto& pixels : member_pixels) {
        for (int pixel : pixels) {
            if (pixel < 0 ||
                static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
                throw std::out_of_range(
                    "member_pixels contains an invalid compact pixel id.");
            }
        }
    }
}

double tilePixelPrize(
    int tile,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels) {
    if (tile < 0 || static_cast<std::size_t>(tile) >= member_pixels.size()) {
        throw std::out_of_range("tilePixelPrize: invalid tile id.");
    }

    double prize = 0.0;
    std::unordered_set<int> seen;
    for (int pixel : member_pixels[tile]) {
        if (pixel < 0 ||
            static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
            throw std::out_of_range("tilePixelPrize: invalid pixel id.");
        }
        if (seen.insert(pixel).second) {
            prize += pixel_probs[pixel];
        }
    }
    return prize;
}

double pathPixelPrize(
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<int>& path) {
    std::vector<char> covered(pixel_probs.size(), 0);
    double prize = 0.0;

    for (int tile : path) {
        if (tile < 0 ||
            static_cast<std::size_t>(tile) >= member_pixels.size()) {
            throw std::out_of_range("pathPixelPrize: invalid tile id.");
        }
        for (int pixel : member_pixels[tile]) {
            if (pixel < 0 ||
                static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
                throw std::out_of_range("pathPixelPrize: invalid pixel id.");
            }
            if (!covered[pixel]) {
                covered[pixel] = 1;
                prize += pixel_probs[pixel];
            }
        }
    }
    return prize;
}

double pathsPixelPrize(
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<std::vector<int>>& paths ) {
    std::vector<char> covered(pixel_probs.size(), 0);
    double prize = 0.0;

    for (const auto& path : paths) {
        for (int tile : path) {
            if (tile < 0 ||
                static_cast<std::size_t>(tile) >= member_pixels.size()) {
                throw std::out_of_range("pathsPixelPrize: invalid tile id.");
            }
            for (int pixel : member_pixels[tile]) {
                if (pixel < 0 ||
                    static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
                    throw std::out_of_range(
                        "pathsPixelPrize: invalid pixel id.");
                }
                if (!covered[pixel]) {
                    covered[pixel] = 1;
                    prize += pixel_probs[pixel];
                }
            }
        }
    }
    return prize;
}

std::vector<int> pixelCoverCounts(
    const std::vector<int>& path,
    const std::vector<std::vector<int>>& member_pixels,
    std::size_t pixel_count ) {
    std::vector<int> counts(pixel_count, 0);
    std::vector<char> node_seen(member_pixels.size(), 0);

    for (int tile : path) {
        if (tile < 0 ||
            static_cast<std::size_t>(tile) >= member_pixels.size()) {
            throw std::out_of_range("pixelCoverCounts: invalid tile id.");
        }
        if (node_seen[tile]) {
            continue;
        }
        node_seen[tile] = 1;
        updatePixelCoverCounts(tile, 1, member_pixels, counts);
    }
    return counts;
}

void updatePixelCoverCounts(
    int tile,
    int delta,
    const std::vector<std::vector<int>>& member_pixels,
    std::vector<int>& cover_counts ) {
    if (tile < 0 || static_cast<std::size_t>(tile) >= member_pixels.size()) {
        throw std::out_of_range("updatePixelCoverCounts: invalid tile id.");
    }
    if (delta != -1 && delta != 1) {
        throw std::invalid_argument(
            "updatePixelCoverCounts: delta must be -1 or 1.");
    }

    std::unordered_set<int> seen;
    for (int pixel : member_pixels[tile]) {
        if (pixel < 0 ||
            static_cast<std::size_t>(pixel) >= cover_counts.size()) {
            throw std::out_of_range(
                "updatePixelCoverCounts: invalid pixel id.");
        }
        if (!seen.insert(pixel).second) {
            continue;
        }
        if (delta < 0 && cover_counts[pixel] == 0) {
            throw std::logic_error(
                "updatePixelCoverCounts: pixel count would become negative.");
        }
        cover_counts[pixel] += delta;
    }
}

double marginalPixelPrize(
    int tile,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<int>& cover_counts ) {
    if (cover_counts.size() != pixel_probs.size()) {
        throw std::invalid_argument(
            "marginalPixelPrize: count/probability size mismatch.");
    }
    if (tile < 0 || static_cast<std::size_t>(tile) >= member_pixels.size()) {
        throw std::out_of_range("marginalPixelPrize: invalid tile id.");
    }

    double gain = 0.0;
    std::unordered_set<int> seen;
    for (int pixel : member_pixels[tile]) {
        if (pixel < 0 ||
            static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
            throw std::out_of_range("marginalPixelPrize: invalid pixel id.");
        }
        if (seen.insert(pixel).second && cover_counts[pixel] == 0) {
            gain += pixel_probs[pixel];
        }
    }
    return gain;
}

double exclusivePixelPrize(
    int tile,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<int>& cover_counts ) {
    if (cover_counts.size() != pixel_probs.size()) {
        throw std::invalid_argument(
            "exclusivePixelPrize: count/probability size mismatch.");
    }
    if (tile < 0 || static_cast<std::size_t>(tile) >= member_pixels.size()) {
        throw std::out_of_range("exclusivePixelPrize: invalid tile id.");
    }

    double loss = 0.0;
    std::unordered_set<int> seen;
    for (int pixel : member_pixels[tile]) {
        if (pixel < 0 ||
            static_cast<std::size_t>(pixel) >= pixel_probs.size()) {
            throw std::out_of_range("exclusivePixelPrize: invalid pixel id.");
        }
        if (seen.insert(pixel).second && cover_counts[pixel] == 1) {
            loss += pixel_probs[pixel];
        }
    }
    return loss;
}

std::vector<double> contextualPixelValues(
    const std::vector<int>& path,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels ) {
    const auto counts = pixelCoverCounts(
        path, member_pixels, pixel_probs.size());
    std::vector<char> in_path(member_pixels.size(), 0);
    for (int tile : path) {
        if (tile >= 0 && static_cast<std::size_t>(tile) < in_path.size()) {
            in_path[tile] = 1;
        }
    }

    std::vector<double> values(member_pixels.size(), 0.0);
    for (std::size_t tile = 0; tile < member_pixels.size(); ++tile) {
        values[tile] = in_path[tile]
            ? exclusivePixelPrize(
                  static_cast<int>(tile), pixel_probs, member_pixels, counts)
            : marginalPixelPrize(
                  static_cast<int>(tile), pixel_probs, member_pixels, counts);
    }
    return values;
}


// calculate delta cost when removing a node
double removalDeltaCost(const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& path,
                         size_t node_index) {
    if (node_index == 0 || node_index >= path.size() - 1) {
        return 0.0; 
    }
    
    int u = path[node_index - 1];
    int v = path[node_index];
    int w = path[node_index + 1];
    
    return costs[u][v] + costs[v][w] - costs[u][w];
}

// calculate delta cost when inserting a node
double insertionDeltaCost(const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& path,
                         size_t node, size_t insert_index) {
    if (insert_index >= path.size() - 1) {
        return std::numeric_limits<double>::infinity();
    }
    
    int u = path[insert_index];
    int w = path[insert_index + 1];
    
    return costs[u][node] + costs[node][w] - costs[u][w];
}

// Find cheapest insertion 'pos' (between pos-1 and pos), pos ∈ [1, n-1] ensure s-t path
std::pair<int, double> cheapestInsertionOnePath(const std::vector<std::vector<double>>& costs,
                                        const std::vector<int>& path, int node) {
    const int n = (int)path.size();
    if (n < 2) return {-1, std::numeric_limits<double>::infinity()};
    int bestPos = -1; double bestDt = std::numeric_limits<double>::infinity();
    for (int pos = 1; pos <= n-1; ++pos) {   
        int a = path[pos-1], b = path[pos];
        double dt = costs[a][node] + costs[node][b] - costs[a][b];
        if (dt < bestDt) { bestDt = dt; bestPos = pos; }
    }
    return {bestPos, bestDt};
}


// Best insertion on a specific s–t route: insert at index 'pos' (between pos-1 and pos), pos ∈ [1, n-1]
std::pair<int,double> bestInsertionOnRoute(const std::vector<std::vector<double>>& costs,
                                                         const std::vector<int>& route,
                                                         int v) {
    const int n = (int)route.size();
    if (n < 2) return {-1, std::numeric_limits<double>::infinity()};
    int bestPos = -1; double bestDt = std::numeric_limits<double>::infinity();
    for (int pos = 1; pos <= n-1; ++pos) {
        int a = route[pos-1], b = route[pos];
        double dt = costs[a][v] + costs[v][b] - costs[a][b];
        if (dt < bestDt) { bestDt = dt; bestPos = pos; }
    }
    return {bestPos, bestDt};
}

void insert_at(std::vector<int>& path, int idx, int node) {
    path.insert(path.begin() + (idx + 1), node);
}

void remove_at(std::vector<int>& path, int idx) {
    path.erase(path.begin() + idx);
}


// Remove duplicate nodes in a path based, preserving the first occurance
void removeDuplicates(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& path) {
    if (path.size() < 2) return;

    std::unordered_set<int> unique_cities;
    std::vector<int> cleaned_path;
    int end = path.back();
    unique_cities.insert(end);
    for (int node : path) {
        if (unique_cities.insert(node).second) {
            cleaned_path.push_back(node);
        }
    }
    cleaned_path.push_back(end);
    path.swap(cleaned_path);
}


// Remove duplicate nodes in a path based on prize/delta_cost ratio
void removeDuplicatesByRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& path) {
    if (path.size() < 2) return;
    
    int start_node = path.front();
    int end_node = path.back();
    
    std::unordered_map<int, std::vector<size_t>> node_positions;
    
    // Find all positions of each node
    for (size_t i = 0; i < path.size(); ++i) {
        node_positions[path[i]].push_back(i);
    }
    
    std::vector<bool> to_remove(path.size(), false);
    
    // For each node that appears multiple times
    for (const auto& [node, positions] : node_positions) {
        if (positions.size() > 1) {
            double best_ratio = -1.0;
            size_t best_position = positions[0];
            
            // Find the position with the best prize/delta_cost ratio
            for (size_t pos : positions) {
                double delta_cost = removalDeltaCost(costs, path, pos);
                double ratio;
                
                if (delta_cost <= 0) {
                    ratio = std::numeric_limits<double>::max(); 
                } else {
                    ratio = prizes[node] / delta_cost;
                }
                
                if (ratio > best_ratio) {
                    best_ratio = ratio;
                    best_position = pos;
                }
            }
            
            // Mark all other positions for removal
            for (size_t pos : positions) {
                if (pos != best_position) {
                    to_remove[pos] = true;
                }
            }
        }
    }
    
    // Remove marked positions (from back to front to maintain indices)
    for (int i = path.size() - 1; i >= 0; --i) {
        if (to_remove[i]) {
            path.erase(path.begin() + i);
        }
    }
    
    // Ensure start and end cities are preserved
    if (path.empty() || path.front() != start_node) {
        path.insert(path.begin(), start_node);
    }
    if (path.size() < 2 || path.back() != end_node) {
        path.push_back(end_node);
    }
}


// Drop nodes with lowest prize until budget is met
void prunePathByPrize(const std::vector<std::vector<double>>& costs,
                         const std::vector<double>& prizes,
                         std::vector<int>& path,
                         double budget,
                         double& total_cost) {
    while (total_cost > budget && path.size() > 2) {
        auto min_itr = std::min_element(path.begin() + 1, path.end() - 1, [&](int a, int b) {
                                           return prizes[a] < prizes[b];
                                        });
        
        size_t index = min_itr - path.begin();
        
        int u = path[index - 1];
        int v = path[index];
        int w = path[index + 1];
        
        double cost_change = costs[u][v] + costs[v][w] - costs[u][w];
        total_cost -= cost_change;
        
        path.erase(min_itr);
    }
}


// Drop nodes with lowest prize/delta_cost ratio until budget is met
void prunePathByRatio(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      std::vector<int>& path,
                      double budget,
                      double& total_cost) {
    while (total_cost > budget && path.size() > 2) {
        double worst_ratio = std::numeric_limits<double>::max();
        size_t worst_index = 1;
        for (size_t i = 1; i + 1 < path.size(); ++i) {
            double dc = removalDeltaCost(costs, path, i);
            double ratio = (dc <= 0.0) ? std::numeric_limits<double>::max()
                                       : (prizes[path[i]] / dc);

            if (ratio < worst_ratio) { 
                worst_ratio = ratio; 
                worst_index = i; 
            }
        }
        total_cost -= removalDeltaCost(costs, path, worst_index);
        path.erase(path.begin() + worst_index);
    }
}


void prunePathBySegmentRatio(const std::vector<std::vector<double>>& costs,
                             const std::vector<double>& prizes,
                             std::vector<int>& path,
                             double budget,
                             double& total_cost) {
    if (path.size() < 2) return;

    auto route_cost = [&](const std::vector<int>& r)->double {
        double c = 0.0;
        for (size_t i = 0; i + 1 < r.size(); ++i) c += costs[r[i]][r[i+1]];
        return c;
    };

    // If caller didn't set or set a bogus total_cost, recompute once.
    if (!(total_cost > 0.0) || !std::isfinite(total_cost)) {
        total_cost = route_cost(path);
    }

    const double eps = 1e-12;

    // Keep endpoints; remove best prize/Δcost segment each round until within budget (with tol)
    while (total_cost > budget + eps && path.size() > 2) {
        const size_t m = path.size();
        const double over = total_cost - budget;

        // Prefix sums for O(1) segment prize and internal edge costs
        std::vector<double> prefP(m + 1, 0.0);
        for (size_t k = 0; k < m; ++k) prefP[k + 1] = prefP[k] + prizes[path[k]];

        std::vector<double> prefE(m, 0.0); // edges: k->k+1 accumulated at index k+1
        for (size_t k = 0; k + 1 < m; ++k) prefE[k + 1] = prefE[k] + costs[path[k]][path[k + 1]];

        auto segPrize = [&](size_t i, size_t j) -> double {
            return prefP[j + 1] - prefP[i]; // prizes of nodes i..j
        };
        auto internalSegEdgeCost = [&](size_t i, size_t j) -> double {
            // edges i->i+1, ..., j-1->j ; zero if i==j
            return prefE[j] - prefE[i];
        };
        auto segDeltaCost = [&](size_t i, size_t j) -> double {
            // Remove nodes i..j; reconnect a=path[i-1] to b=path[j+1]
            const int a = path[i - 1];
            const int b = path[j + 1];
            const double old_boundary = costs[a][path[i]] + costs[path[j]][b];
            const double old_internal = internalSegEdgeCost(i, j);
            const double new_edge = costs[a][b];
            return (old_boundary + old_internal) - new_edge; // savings (may be ≤ 0)
        };

        double best_ratio = std::numeric_limits<double>::infinity();
        size_t best_i = 1, best_j = 1;
        bool found = false;

        // i and j are inclusive indices of the removable segment; keep endpoints at 0 and m-1
        for (size_t i = 1; i + 1 < m; ++i) {
            const size_t max_j = m - 2;
            double best_local_ratio = std::numeric_limits<double>::infinity();
            size_t best_local_j = i;
            bool local_found = false;

            for (size_t j = i; j <= max_j; ++j) {
                double delta = segDeltaCost(i, j);
                if (delta <= eps) continue;  // no real savings (or worse)
                double ratio = segPrize(i, j) / delta;

                if (ratio < best_local_ratio) {
                    best_local_ratio = ratio;
                    best_local_j = j;
                    local_found = true;
                }
                // Stop early once this growing segment could cover the overage
                if (delta >= over - eps) break;
            }

            if (local_found && best_local_ratio < best_ratio) {
                best_ratio = best_local_ratio;
                best_i = i;
                best_j = best_local_j;
                found = true;
            }
        }

        if (!found || !std::isfinite(best_ratio)) {
            // Nothing yields positive savings → cannot meet budget by pruning
            break;
        }

        // Apply removal
        double delta_cost = segDeltaCost(best_i, best_j);
        if (delta_cost <= eps) break; // defensive

        total_cost -= delta_cost;
        path.erase(path.begin() + static_cast<std::ptrdiff_t>(best_i),
                   path.begin() + static_cast<std::ptrdiff_t>(best_j) + 1);
    }

    // Optional: clamp to zero noise
    if (std::abs(total_cost - route_cost(path)) > 1e-9) {
        total_cost = route_cost(path); // resync if needed
    }
}



// Insert nodes with highest possible value in cheapest position until budget exhausts
void insertHighValueCheapest(std::vector<int>& path,
                                const std::vector<std::vector<double>>& cost,
                                const std::vector<double>& prize,
                                double budget,
                                double& current_cost,
                                std::vector<bool>& visited) {
    const int n = static_cast<int>(cost.size());

    // Build candidate list ordered by decreasing prize
    std::vector<int> cand;
    cand.reserve(n);
    for (int v = 0; v < n; ++v)
        if (!visited[v] && v != path.front() && v != path.back())
            cand.push_back(v);

    std::sort(cand.begin(), cand.end(),
              [&](int a, int b){ return prize[a] > prize[b]; });

    // Try each node repeatedly until no more insertions possible
    bool progress = true;
    while (progress) {
        progress = false;

        for (int idx = 0; idx < static_cast<int>(cand.size()); ++idx) {
            int node = cand[idx];
            if (visited[node]) continue; 

            // Find cheapest insertion position
            double best_extra = std::numeric_limits<double>::infinity();
            std::size_t best_pos = 0;

            for (std::size_t j = 0; j < path.size() - 1; ++j) {
                double extra = cost[path[j]][node] +
                               cost[node][path[j+1]] -
                               cost[path[j]][path[j+1]];
                if (extra < best_extra) {
                    best_extra = extra;
                    best_pos = j + 1;
                }
            }

            if (current_cost + best_extra <= budget + 1e-12) {
                path.insert(path.begin() + best_pos, node);
                current_cost += best_extra;
                visited[node] = true;
                progress = true;
            }
        }
    }
}

// Insert nodes with highest possible value before end node until budget exhausts
void insertHighValueEnd(std::vector<int>& path,
                                  const std::vector<std::vector<double>>& costs,
                                  const std::vector<double>& prizes,
                                  double budget,
                                  double& current_cost,
                                  std::vector<bool>& visited) {
    if (path.size() < 2) return;

    const int n = static_cast<int>(costs.size());
    const double eps = 1e-12;

    // Build list of unvisited cities sorted by decreasing prize (static order).
    std::vector<int> candidates;
    candidates.reserve(n);
    for (int i = 0; i < n; ++i)
        if (!visited[i] && i != path.front() && i != path.back())
            candidates.push_back(i);

    std::sort(candidates.begin(), candidates.end(),
              [&](int a, int b) { return prizes[a] > prizes[b]; });

    // Revisit candidates until no further feasible insertion.
    bool progress = true;
    while (progress) {
        progress = false;

        for (int node : candidates) {
            if (visited[node]) continue;

            double best_extra = std::numeric_limits<double>::infinity();
            std::size_t best_pos = 0;

            // Recompute cheapest insertion w.r.t. the *current* path.
            for (std::size_t j = 0; j + 1 < path.size(); ++j) {
                double extra = costs[path[j]][node]
                             + costs[node][path[j+1]]
                             - costs[path[j]][path[j+1]];
                if (extra < best_extra) {
                    best_extra = extra;
                    best_pos = j + 1;
                }
            }

            // Insert if feasible (allowing zero/negative extra if graph permits).
            if (current_cost + best_extra <= budget + eps) {
                path.insert(path.begin() + static_cast<long>(best_pos), node);
                current_cost += best_extra;
                visited[node] = true;
                progress = true;
            }
        }
    }
}


// Insert nodes with global best prize/delta_cost ratio until budget exhausts
void insertBestRatio(std::vector<int>& path,
                     const std::vector<std::vector<double>>& cost,
                     const std::vector<double>& prize,
                     double budget,
                     double& current_cost,
                     std::vector<bool>& visited) {
    if (path.size() < 2) return;

    const int n = static_cast<int>(cost.size());
    const double eps = 1e-12;
    auto is_endpoint = [&](int v){ return v == path.front() || v == path.back(); };

    while (true) {
        int best_city = -1;
        std::size_t best_pos = 0;
        double best_extra = 0.0;

        // We maximize score = prize / extra (bigger is better).
        // For extra <= eps (free or cost-reducing), treat as +inf score.
        double best_score = -std::numeric_limits<double>::infinity();
        double best_prize = -std::numeric_limits<double>::infinity();

        for (int city = 0; city < n; ++city) {
            if (city < 0 || city >= n) continue;          // defensive
            if (visited[city] || is_endpoint(city)) continue;

            for (std::size_t j = 0; j + 1 < path.size(); ++j) {
                double extra = cost[path[j]][city]
                             + cost[city][path[j+1]]
                             - cost[path[j]][path[j+1]];

                if (std::isnan(extra)) continue;          //bad data
                if (current_cost + extra > budget + eps) continue;

                // Score: prefer cost-reducing/free insertions; otherwise prize/extra
                double score = (extra <= eps) ? std::numeric_limits<double>::infinity()
                                              : (prize[city] / extra);

                // Tie-break: higher score, then higher prize, then smaller city id
                bool better = false;
                if (score > best_score) {
                    better = true;
                } else if (std::abs(score - best_score) <= eps) {
                    if (prize[city] > best_prize + eps) {
                        better = true;
                    } else if (std::abs(prize[city] - best_prize) <= eps && city < best_city) {
                        better = true;
                    }
                }

                if (better) {
                    best_city  = city;
                    best_pos   = j + 1;
                    best_extra = extra;
                    best_score = score;
                    best_prize = prize[city];
                }
            }
        }

        if (best_city == -1) break;  // no feasible insertion

        path.insert(path.begin() + static_cast<long>(best_pos), best_city);
        current_cost += best_extra;
        visited[best_city] = true;
    }
}

void repairPath(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                std::vector<int>& path,
                double budget) {
    if (path.size() < 2) return;
    removeDuplicates(costs, prizes, path);
    double total_cost = pathCost(costs, path);
    prunePathByRatio(costs, prizes, path, budget, total_cost);

    std::vector<bool> visited(costs.size(), false);
    for (auto p : path) visited[p] = true;
    insertHighValueCheapest(path, costs, prizes, budget, total_cost, visited);
}

//overload
void repairPath(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                std::vector<int>& path,
                double budget, 
                std::vector<bool>& visited) {
    if (path.size() < 2) return;
    removeDuplicates(costs, prizes, path);
    double total_cost = pathCost(costs, path);
    prunePathByRatio(costs, prizes, path, budget, total_cost);

    std::fill(visited.begin(), visited.end(), false);
    for (auto p : path) visited[p] = true;
    insertHighValueCheapest(path, costs, prizes, budget, total_cost, visited);
}


void validatePixelPathInputs(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<int>& path,
    double budget ) {
    const std::size_t node_count = costs.size();
    if (node_count < 2 || !std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "overlap-aware path repair received invalid dimensions or budget.");
    }
    for (const auto& row : costs) {
        if (row.size() != node_count) {
            throw std::invalid_argument("cost matrix must be square.");
        }
    }
    validatePixelCoverageData(node_count, pixel_probs, member_pixels);
    if (path.size() < 2) {
        throw std::invalid_argument("an s-t path must contain both endpoints.");
    }
    for (int node : path) {
        if (node < 0 || static_cast<std::size_t>(node) >= node_count) {
            throw std::out_of_range("path contains an invalid node id.");
        }
    }
}

void removeDuplicateInteriorNodes(std::vector<int>& path) {
    const int start = path.front();
    const int end = path.back();
    std::vector<char> seen;
    int largest = std::max(start, end);
    for (int node : path) {
        largest = std::max(largest, node);
    }
    seen.assign(static_cast<std::size_t>(largest + 1), 0);
    seen[start] = 1;
    seen[end] = 1;

    std::vector<int> cleaned;
    cleaned.reserve(path.size());
    cleaned.push_back(start);
    for (std::size_t index = 1; index + 1 < path.size(); ++index) {
        const int node = path[index];
        if (!seen[node]) {
            seen[node] = 1;
            cleaned.push_back(node);
        }
    }
    cleaned.push_back(end);
    path.swap(cleaned);
}

void repairPixelPathImpl(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    std::vector<int>& path,
    double budget,
    std::vector<bool>* visited_out ) {
    validatePixelPathInputs(
        costs, pixel_probs, member_pixels, path, budget);
    removeDuplicateInteriorNodes(path);

    constexpr double tolerance = 1e-12;
    double total_cost = pathCost(costs, path);

    // A selected tile's removal loss is only the probability of pixels that
    // no other currently selected tile covers.
    while (total_cost > budget + tolerance && path.size() > 2) {
        const auto counts = pixelCoverCounts(
            path, member_pixels, pixel_probs.size());

        std::size_t remove_index = path.size();
        double best_ratio = std::numeric_limits<double>::infinity();
        double best_saving = -std::numeric_limits<double>::infinity();

        for (std::size_t index = 1; index + 1 < path.size(); ++index) {
            const double saving = removalDeltaCost(costs, path, index);
            if (saving <= tolerance) {
                continue;
            }
            const double loss = exclusivePixelPrize(
                path[index], pixel_probs, member_pixels, counts);
            const double ratio = loss / saving;
            if (ratio < best_ratio - tolerance ||
                (std::abs(ratio - best_ratio) <= tolerance &&
                 saving > best_saving + tolerance)) {
                best_ratio = ratio;
                best_saving = saving;
                remove_index = index;
            }
        }

        if (remove_index == path.size()) {
            break;
        }
        total_cost -= removalDeltaCost(costs, path, remove_index);
        path.erase(path.begin() + static_cast<std::ptrdiff_t>(remove_index));
    }

    std::vector<bool> visited(costs.size(), false);
    for (int node : path) {
        visited[node] = true;
    }

    if (total_cost <= budget + tolerance) {
        auto counts = pixelCoverCounts(
            path, member_pixels, pixel_probs.size());

        while (true) {
            int best_node = -1;
            std::size_t best_position = 0;
            double best_extra = 0.0;
            double best_gain = 0.0;
            double best_score = -std::numeric_limits<double>::infinity();

            for (int node = 0; node < static_cast<int>(costs.size()); ++node) {
                if (visited[node] || node == path.front() || node == path.back()) {
                    continue;
                }

                const double gain = marginalPixelPrize(
                    node, pixel_probs, member_pixels, counts);
                if (gain <= tolerance) {
                    continue;
                }

                const auto [position, extra] =
                    cheapestInsertionOnePath(costs, path, node);
                if (position < 0 ||
                    total_cost + extra > budget + tolerance) {
                    continue;
                }

                const double score = extra <= tolerance
                    ? std::numeric_limits<double>::infinity()
                    : gain / extra;
                const bool better =
                    score > best_score + tolerance ||
                    (std::abs(score - best_score) <= tolerance &&
                     (gain > best_gain + tolerance ||
                      (std::abs(gain - best_gain) <= tolerance &&
                       (extra < best_extra - tolerance || best_node < 0))));

                if (better) {
                    best_node = node;
                    best_position = static_cast<std::size_t>(position);
                    best_extra = extra;
                    best_gain = gain;
                    best_score = score;
                }
            }

            if (best_node < 0) {
                break;
            }
            path.insert(
                path.begin() + static_cast<std::ptrdiff_t>(best_position),
                best_node);
            total_cost += best_extra;
            visited[best_node] = true;
            updatePixelCoverCounts(
                best_node, 1, member_pixels, counts);
        }
    }

    if (visited_out != nullptr) {
        *visited_out = std::move(visited);
    }
}


void repairPath(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    std::vector<int>& path,
    double budget ) {
    repairPixelPathImpl(
        costs, pixel_probs, member_pixels, path, budget, nullptr);
}

void repairPath(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    std::vector<int>& path,
    double budget,
    std::vector<bool>& visited ) {
    repairPixelPathImpl(
        costs, pixel_probs, member_pixels, path, budget, &visited);
}


std::vector<int> repair_path(const std::vector<int>& path,
               const std::vector<std::vector<double>>& costs,
               const std::vector<double>& prizes,
               double budget, std::vector<bool>& visited) {
    //initial state
    std::vector<int> new_path = path;                 // mutable copy
    const int  s        = path.front();               // start
    const int  t        = path.back();                // end (fixed)
    auto edge = [&](int u,int v){ return costs[u][v]; };

    auto tour_cost = [&](const std::vector<int>& P){
        double c = 0.0;
        for (std::size_t i = 1; i < P.size(); ++i) c += edge(P[i-1], P[i]);
        return c;
    };

    double current_cost = tour_cost(new_path);
    // if (current_cost <= budget) return new_path; 

    // remove nodes with smallest prize until cost ≤ budget
    while (current_cost > budget && new_path.size() > 2) {
        std::size_t worst_pos   = std::numeric_limits<std::size_t>::max();
        double      min_prize   = std::numeric_limits<double>::max();

        for (std::size_t i = 1; i < new_path.size() - 1; ++i) {  
            int v = new_path[i];
            if (prizes[v] < min_prize) {
                min_prize = prizes[v];
                worst_pos = i;
            }
        }
        if (worst_pos == std::numeric_limits<std::size_t>::max()) break;

        int v_rem = new_path[worst_pos];
        new_path.erase(new_path.begin() + worst_pos);
        visited[v_rem] = false;

        current_cost = tour_cost(new_path); 
    }

    insertHighValueCheapest(new_path, costs, prizes, budget, current_cost, visited);

    return new_path; 
}


/*****Functions for TOP Paths Operations ******/

std::vector<double> multiPathsCost(const std::vector<std::vector<double>>& costs,
                                   const std::vector<std::vector<int>>& paths) {
    std::vector<double> paths_costs;
    paths_costs.reserve(paths.size());
    for (const auto& p : paths) 
        paths_costs.push_back(pathCost(costs, p));
    return paths_costs;
}

double multiPathsSumCost(const std::vector<std::vector<double>>& costs,
                       const std::vector<std::vector<int>>& paths) {
    double sum = 0.0;
    for (const auto& p : paths) sum += pathCost(costs, p);
    return sum;
}

double multiPathsMaxCost(const std::vector<std::vector<double>>& costs,
                       const std::vector<std::vector<int>>& paths) {
    double max_cost = 0.0;
    for (const auto& p : paths) max_cost = std::max(max_cost, pathCost(costs, p));
    return max_cost;
}

std::vector<double> multiPathsPrize(const std::vector<double>& prizes,
                       const std::vector<std::vector<int>>& paths) {
    std::vector<double> paths_prizes;
    paths_prizes.reserve(paths.size());
    for (const auto& p : paths) 
        paths_prizes.push_back(pathPrize(prizes, p));
    return paths_prizes;
}

double multiPathsSumPrize(const std::vector<double>& prizes,
                       const std::vector<std::vector<int>>& paths) {
    double sum = 0.0;
    for (const auto& p : paths) sum += pathPrize(prizes, p);
    return sum;
}

std::vector<bool> getVisited(const std::vector<std::vector<int>>& paths, int nNodes) {
    std::vector<bool> visited(nNodes, false);
    for (const auto& r : paths)
        for (int i = 1; i+1 < (int)r.size(); ++i)
            visited[r[i]] = true;
    return visited;
}

// Returns: {best_route, best_position, best_cost_increase} or {-1, -1, infinity} if no feasible insertion
std::tuple<int, int, double> cheapestInsertionAll(const std::vector<std::vector<double>>& costs,
                                                       const std::vector<std::vector<int>>& paths,
                                                       const std::vector<double>& pathCosts,
                                                       const std::vector<double>& prizes,
                                                       double budget, 
                                                       int node) {
    int best_route = -1;
    int best_pos = -1;
    double best_cost_increase = std::numeric_limits<double>::max();
    double best_heuristic = 0.0;
    
    for (int r = 0; r < paths.size(); r++) {
        const auto& path = paths[r];
        if (path.size() < 2) continue;
        
        double currentCost = pathCosts[r];
        
        for (int pos = 1; pos < path.size(); pos++) {
            int i = path[pos - 1];
            int j = path[pos];
            
            double costIncrease = costs[i][node] + costs[node][j] - costs[i][j];
            double newCost = currentCost + costIncrease;
            
            if (newCost <= budget && costIncrease > 0) {
                double heuristic = prizes[node] / costIncrease;
                
                if (costIncrease < best_cost_increase || 
                    (std::abs(costIncrease - best_cost_increase) < 1e-12 && heuristic > best_heuristic)) {
                    best_route = r;
                    best_pos = pos;
                    best_cost_increase = costIncrease;
                    best_heuristic = heuristic;
                }
            }
        }
    }
    
    return {best_route, best_pos, best_cost_increase};
}

bool is_valid_node(int v, int n) {
    return (v >= 0 && v < n);
}

// Remove duplicates across paths; keep the FIRST ownership
void removeInterPathDuplicatesKeepFirst(const std::vector<std::vector<double>>& costs,
                                        std::vector<std::vector<int>>& paths) {
    const int n = (int)costs.size();
    std::unordered_set<int> owned;
    for (auto& path : paths) {
        if (path.size() < 2) continue;
        const int s = path.front(), t = path.back();
        std::vector<int> out; out.reserve(path.size());
        out.push_back(s);
        for (std::size_t i = 1; i + 1 < path.size(); ++i) {
            int v = path[i];
            if (!is_valid_node(v, n)) continue;
            if (owned.insert(v).second) out.push_back(v);
        }
        out.push_back(t);
        path.swap(out);
    }
}

// Keep the single best (prize / delta_cost) occurrence across all paths per node
void removeInterPathDuplicatesByRatio(const std::vector<std::vector<double>>& costs,
                                      const std::vector<double>& prizes,
                                      std::vector<std::vector<int>>& paths) {

    double EPS = 1e-12;
    const int n = (int)costs.size();
    if ((int)prizes.size() != n) return;

    struct Occurrence { std::size_t p_idx, pos; double ratio, delta; };
    std::unordered_map<int, std::vector<Occurrence>> occs;

    // gather
    for (std::size_t p = 0; p < paths.size(); ++p) {
        const auto& path = paths[p];
        if (path.size() < 2) continue;
        for (std::size_t i = 1; i + 1 < path.size(); ++i) {
            int v = path[i];
            if (!is_valid_node(v, n)) continue;
            double dc = removalDeltaCost(costs, path, i);
            double r  = (dc <= EPS) ? std::numeric_limits<double>::infinity()
                                    : (prizes[v] / dc);
            occs[v].push_back({p, i, r, dc});
        }
    }

    // choose keep occurrence for each node
    std::vector<std::vector<char>> to_remove(paths.size());
    for (std::size_t p = 0; p < paths.size(); ++p)
        to_remove[p].assign(paths[p].size(), 0);

    for (auto& kv : occs) {
        const int node = kv.first;
        auto& list = kv.second;
        if (list.size() <= 1) continue;

        std::size_t keep = 0;
        for (std::size_t i = 1; i < list.size(); ++i) {
            const auto& a = list[keep];
            const auto& b = list[i];
            bool better = false;
            if (b.ratio > a.ratio) better = true;
            else if (std::abs(b.ratio - a.ratio) <= EPS) {
                // prefer keeping the one whose removal hurts more (smaller delta)
                if (b.delta < a.delta - EPS) better = true;
                else if (std::abs(b.delta - a.delta) <= EPS) {
                    if (b.p_idx < a.p_idx) better = true;
                    else if (b.p_idx == a.p_idx && b.pos < a.pos) better = true;
                }
            }
            if (better) keep = i;
        }
        for (std::size_t i = 0; i < list.size(); ++i) {
            if (i == keep) continue;
            const auto& oc = list[i];
            if (oc.pos == 0 || oc.pos + 1 >= paths[oc.p_idx].size()) continue;
            to_remove[oc.p_idx][oc.pos] = 1;
        }
    }

    // erase back-to-front per path
    for (std::size_t p = 0; p < paths.size(); ++p) {
        auto& path = paths[p];
        if (path.size() < 2) continue;
        for (std::ptrdiff_t i = (std::ptrdiff_t)path.size() - 2; i >= 1; --i) {
            if (to_remove[p][(std::size_t)i])
                path.erase(path.begin() + i);
        }
    }
}


// Greedy global insertion: pick (city, route, pos) with best prize/delta ratio each time.
void insertHighValueCheapestMulti(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& prizes,
    std::vector<std::vector<int>>& paths,
    double budget_per_route,
    int s, int t) {

    double EPS = 1e-12;

    const int n = (int)costs.size();
    auto is_endpoint = [&](int v){ return v == s || v == t; };

    // compute current route costs (by value, not ref to temporary!)
    std::vector<double> route_costs = multiPathsCost(costs, paths);

    // track used interior nodes globally
    std::vector<bool> used(n, false);
    for (const auto& r : paths)
        for (std::size_t i = 1; i + 1 < r.size(); ++i) used[r[i]] = true;

    double best_score, best_extra;
    int best_city, best_rid;
    std::size_t best_pos;

    while (true) {
        best_score = -std::numeric_limits<double>::infinity();
        best_extra = 0.0;
        best_city  = -1;
        best_rid   = -1;
        best_pos   = 0;

        for (int city = 0; city < n; ++city) {
            if (used[city] || is_endpoint(city)) continue;

            for (int rid = 0; rid < (int)paths.size(); ++rid) {
                const auto& r = paths[rid];
                if (r.size() < 2) continue;

                for (std::size_t j = 0; j + 1 < r.size(); ++j) {
                    double extra = costs[r[j]][city] + costs[city][r[j+1]] - costs[r[j]][r[j+1]];
                    if (std::isnan(extra)) continue;
                    if (route_costs[rid] + extra > budget_per_route + EPS) continue;

                    // prefer negative/zero extra as +inf score
                    double score = (extra <= EPS) ? std::numeric_limits<double>::infinity()
                                                  : (prizes[city] / extra);

                    bool better = (score > best_score) ||
                                  (std::abs(score - best_score) <= 1e-15 && (
                                      prizes[city] > (best_city >= 0 ? prizes[best_city] : -1e300) ||
                                      (std::abs(prizes[city] - (best_city >= 0 ? prizes[best_city] : 0.0)) <= 1e-15 && (
                                          extra < best_extra ||
                                          (std::abs(extra - best_extra) <= 1e-15 && (
                                              r.size() < (best_rid >= 0 ? paths[best_rid].size() : (std::size_t)-1) ||
                                              (r.size() == (best_rid >= 0 ? paths[best_rid].size() : r.size()) && city < best_city)
                                          ))
                                      ))
                                  ));

                    if (better) {
                        best_score = score;
                        best_city  = city;
                        best_rid   = rid;
                        best_pos   = j + 1;
                        best_extra = extra;
                    }
                }
            }
        }

        if (best_city == -1) break; // no feasible insertion left

        paths[best_rid].insert(paths[best_rid].begin() + (long)best_pos, best_city);
        route_costs[best_rid] += best_extra;
        used[best_city] = true;
    }
}

// full multi-route repair
void repairMultiPaths(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      std::vector<std::vector<int>>& paths,
                      double budget_per_route,
                      int s, int t) {

    // inter-path dedup
    removeInterPathDuplicatesKeepFirst(costs, paths);

    // per-route prune to budget
    std::vector<double> route_costs = multiPathsCost(costs, paths);
    for (std::size_t i = 0; i < paths.size(); ++i) {
        prunePathByRatio(costs, prizes, paths[i], budget_per_route, route_costs[i]);
        // prunePathBySegmentRatio(costs, prizes, paths[i], budget_per_route, route_costs[i]);
    }

    // global greedy insertion (unique interiors, s/t only at ends)
    insertHighValueCheapestMulti(costs, prizes, paths, budget_per_route, s, t);
}


// Greedy global insertion: best (city, route, pos) by prize/delta each time.
// Adapted for multiple starts: start_indices[rid] is the fixed start of route rid.
// We still assume a shared end t for all routes.
void insertHighValueCheapestMulti(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& prizes,
    std::vector<std::vector<int>>& paths,
    double budget_per_route,
    const std::vector<int>& start_indices,
    int t) {
    using std::size_t;
    const double EPS = 1e-12;

    const int n = (int)costs.size();
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    if ((int)start_indices.size() != (int)paths.size())
        throw std::runtime_error("start_indices size must match paths.size()");

    // Compute current route costs
    std::vector<double> route_costs = multiPathsCost(costs, paths);

    // Mark fixed endpoints: all provided starts, shared t, plus whatever endpoints already in paths.
    std::vector<bool> fixed_endpoint(n, false);
    for (int s : start_indices) {
        if (s < 0 || s >= n) throw std::runtime_error("bad start index");
        fixed_endpoint[s] = true;
    }
    if (t < 0 || t >= n) throw std::runtime_error("bad t");
    fixed_endpoint[t] = true;
    for (const auto& r : paths) {
        if (!r.empty()) {
            fixed_endpoint[r.front()] = true;
            fixed_endpoint[r.back()]  = true;
        }
    }

    // Track globally used interior nodes (kept unique across routes)
    std::vector<bool> used(n, false);
    for (const auto& r : paths)
        for (size_t i = 1; i + 1 < r.size(); ++i) used[r[i]] = true;

    while (true) {
        double best_score = -std::numeric_limits<double>::infinity();
        double best_extra = 0.0;
        int    best_city  = -1;
        int    best_rid   = -1;
        size_t best_pos   = 0;

        for (int city = 0; city < n; ++city) {
            if (used[city] || fixed_endpoint[city]) continue;

            for (int rid = 0; rid < (int)paths.size(); ++rid) {
                const auto& r = paths[rid];
                if (r.size() < 2) continue;

                for (size_t j = 0; j + 1 < r.size(); ++j) {
                    const int a = r[j], b = r[j+1];
                    double extra = costs[a][city] + costs[city][b] - costs[a][b];
                    if (std::isnan(extra)) continue;
                    if (route_costs[rid] + extra > budget_per_route + EPS) continue;

                    // Score = prize / extra; favor nonpositive extra
                    double score = (extra <= EPS) ? std::numeric_limits<double>::infinity()
                                                  : (prizes[city] / extra);

                    bool better =
                        (score > best_score) ||
                        (std::abs(score - best_score) <= 1e-15 && (
                            prizes[city] > (best_city >= 0 ? prizes[best_city] : -1e300) ||
                            (std::abs(prizes[city] - (best_city >= 0 ? prizes[best_city] : 0.0)) <= 1e-15 && (
                                extra < best_extra ||
                                (std::abs(extra - best_extra) <= 1e-15 && (
                                    r.size() < (best_rid >= 0 ? paths[best_rid].size()
                                                              : std::numeric_limits<size_t>::max()) ||
                                    (r.size() == (best_rid >= 0 ? paths[best_rid].size() : r.size()) &&
                                     city < best_city)
                                ))
                            ))
                        ));

                    if (better) {
                        best_score = score;
                        best_city  = city;
                        best_rid   = rid;
                        best_pos   = j + 1;
                        best_extra = extra;
                    }
                }
            }
        }

        if (best_city == -1) break; // no feasible insertion left

        paths[best_rid].insert(paths[best_rid].begin() + (long)best_pos, best_city);
        route_costs[best_rid] += best_extra;
        used[best_city] = true;
    }
}

// Full multi-route repair with multiple starts.
// Keeps routes within budget, de-dups across routes, then does global greedy insertions.
void repairMultiPaths(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      std::vector<std::vector<int>>& paths,
                      double budget_per_route,
                      const std::vector<int>& start_indices,
                      int t) {
    if ((int)start_indices.size() != (int)paths.size())
        throw std::runtime_error("start_indices size must match paths.size()");

    // Inter-path dedup (keeps first occurrence)
    removeInterPathDuplicatesKeepFirst(costs, paths);

    // Per-route prune to budget (assumes prune keeps endpoints intact)
    std::vector<double> route_costs = multiPathsCost(costs, paths);
    for (std::size_t i = 0; i < paths.size(); ++i) {
        prunePathByRatio(costs, prizes, paths[i], budget_per_route, route_costs[i]);
    }

    // Global greedy insertion (unique interiors; all starts and common t are fixed endpoints)
    insertHighValueCheapestMulti(costs, prizes, paths, budget_per_route, start_indices, t);
}



bool verify_routes(const std::vector<std::vector<double>>& costs,
                   const std::vector<double>& prizes,
                   int s, int t,
                   const std::vector<std::vector<int>>& routes,
                   double Budget,
                   int m_expected,  
                   double tol) {

    bool ok = true;
    int N = costs.size();

    if (m_expected >= 0 && (int)routes.size() != m_expected) {
        std::cout << "[verify] expected " << m_expected
                  << " routes, got " << routes.size() << "\n";
        ok = false;
    }

    std::unordered_set<int> seen_global;   // all intermediate nodes used

    auto route_cost = [&](const std::vector<int>& r)->double{
        double c=0.0;
        for (size_t i=0;i+1<r.size();++i) c += costs[r[i]][r[i+1]];
        return c;
    };

    for (size_t r_idx = 0; r_idx < routes.size(); ++r_idx) {
        const auto& path = routes[r_idx];

        // ---- structural checks ----
        if (path.empty()) {
            std::cout << "[verify] route " << r_idx << " is empty\n";
            ok = false; continue;
        }
        if (path.front() != s || path.back() != t) {
            std::cout << "[verify] route " << r_idx
                      << " does not start at s or end at t\n";
            ok = false;
        }

        // ---- budget check ----
        double cost = route_cost(path);
        if (cost > Budget + tol) {
            std::cout << "[verify] route " << r_idx << " cost "
                      << cost << " exceeds Budget " << Budget << "\n";
            ok = false;
        }

        // ---- duplicate checks ----
        std::unordered_set<int> seen_local;
        for (size_t pos = 0; pos < path.size(); ++pos) {
            int v = path[pos];
            if (v < 0 || v >= N) {
                std::cout << "[verify] node id " << v
                          << " out of range in route " << r_idx << "\n";
                ok = false;
            }
            if (v == s || v == t) continue;

            if (!seen_local.insert(v).second) {               // appears twice in same route
                std::cout << "[verify] node " << v
                          << " appears twice in route " << r_idx << "\n";
                ok = false;
            }
        }
        for (int v : seen_local) {
            if (!seen_global.insert(v).second) {              // appears in another route
                std::cout << "[verify] node " << v
                          << " appears in multiple routes\n";
                ok = false;
            }
        }
    }

    if (ok)
        std::cout << "[verify] all checks passed for "
                  << routes.size() << " route(s)\n";
    return ok;
}

