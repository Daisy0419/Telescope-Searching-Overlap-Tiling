#include "LocalSearch.h"
#include "PathOperations.h"

#include <vector>
#include <array>
#include <algorithm>
#include <random>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <random>
#include <numeric>
#include <limits>
#include <iostream>

/*****one-pass local improvement, 2-opt, 3-opt, or-opt*****/
// 2-opt, fixed start to fixed end

// bool two_opt_st(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
//     if (costs.empty() || costs.size() != costs[0].size()) {
//         throw std::invalid_argument("Cost matrix must be non-empty and square.");
//         return false;
//     }
//     if (path.size() > costs.size()) {
//         for(int p : path) {
//             std::cout << p << " ";
//         }
//         std::cout << "\n";
//         throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
//         return false;
//     }
//     if (path.size() < 4) {
//         return false;  
//     }
//     bool is_improved = false;
//     // Avoid modifying edges adjacent to start (s) or end (t)
//     for (size_t i = 1; i < path.size() - 2; ++i) {
//         for (size_t j = i + 1; j < path.size() - 1; ++j) {
//             int a = path[i - 1];
//             int b = path[i];
//             int c = path[j];
//             int d = path[j + 1];

//             if (i == 1 && j + 1 == path.size() - 1) continue; // affect both s and t
//             if (i == 1 && a == path[0]) continue;             // affects s
//             if (j + 1 == path.size() - 1 && d == path.back()) continue; // affects t

//             double original_cost = costs[a][b] + costs[c][d];
//             double uncrossed_cost = costs[a][c] + costs[b][d];

//             if (uncrossed_cost < original_cost - 1e-9) {
//                 std::reverse(path.begin() + i, path.begin() + j + 1);
//                 is_improved = true;
//             }
//         }
//     }
//     return is_improved;
// }

bool two_opt_st(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    if (costs.empty() || costs.size() != costs[0].size()) {
        throw std::invalid_argument("Cost matrix must be non-empty and square.");
        return false;
    }
    if (path.size() > costs.size()) {
        for(int p : path) {
            std::cout << p << " ";
        }
        std::cout << "\n";
        throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
        return false;
    }
    if (path.size() < 4) {
        return false;  
    }
    // Edge costs may be asymmetric after dwell time is moved onto the
    // destination edge. Include the direction change of every reversed
    // internal edge; the standard symmetric 2-opt boundary formula alone is
    // not valid for this graph.
    double best_delta = -1e-9;
    std::size_t best_i = 0;
    std::size_t best_j = 0;

    for (std::size_t i = 1; i + 1 < path.size(); ++i) {
        double internal_delta = 0.0;
        for (std::size_t j = i + 1; j + 1 < path.size(); ++j) {
            internal_delta +=
                costs[path[j]][path[j - 1]] -
                costs[path[j - 1]][path[j]];

            const double boundary_delta =
                costs[path[i - 1]][path[j]] +
                costs[path[i]][path[j + 1]] -
                costs[path[i - 1]][path[i]] -
                costs[path[j]][path[j + 1]];
            const double delta = boundary_delta + internal_delta;
            if (delta < best_delta) {
                best_delta = delta;
                best_i = i;
                best_j = j;
            }
        }
    }

    if (best_i == 0) {
        return false;
    }
    std::reverse(path.begin() + static_cast<std::ptrdiff_t>(best_i),
                 path.begin() + static_cast<std::ptrdiff_t>(best_j + 1));
    return true;
}

// 3-opt for an s-t path (2-opt case not included )
bool three_opt_st(std::vector<int>& path,
                const std::vector<std::vector<double>>& cost) {
    const int n = static_cast<int>(path.size());
    if (n < 6) return false;  

    auto seg_cost = [&](int u, int v) { return cost[path[u]][path[v]]; };

    double best_gain = 0.0;
    int bi = -1, bj = -1, bk = -1, best_case = -1; 

    for (int i = 1; i <= n - 5; ++i)
        for (int j = i + 1; j <= n - 4; ++j)
            for (int k = j + 1; k <= n - 3; ++k) {
                // edges removed: (i-1,i) , (j,j+1) , (k,k+1)
                const double rem =
                    seg_cost(i - 1, i) + seg_cost(j, j + 1) + seg_cost(k, k + 1);

                // P  : reverse both A=i..j and B=j+1..k
                const double addP =
                    seg_cost(i - 1, j) + seg_cost(i, k) + seg_cost(j + 1, k + 1);

                // X1 : B then A (no reversals)
                const double addX1 =
                    seg_cost(i - 1, j + 1) + seg_cost(k, i) + seg_cost(j, k + 1);

                // X2 : B then A reversed
                const double addX2 =
                    seg_cost(i - 1, j + 1) + seg_cost(k, j) + seg_cost(i, k + 1);

                // X3 : B reversed then A
                const double addX3 =
                    seg_cost(i - 1, k) + seg_cost(j + 1, i) + seg_cost(j, k + 1);

                const double gains[4] = { rem - addP, rem - addX1,
                                          rem - addX2, rem - addX3 };

                int best_here = std::max_element(gains, gains + 4) - gains;
                if (gains[best_here] > best_gain + 1e-12) {
                    best_gain = gains[best_here];
                    bi = i; bj = j; bk = k; best_case = best_here;
                }
            }

    if (best_case == -1) return false; 

    // ----- apply the chosen 3-opt move -------------------------------------
    switch (best_case) {
        case 0:  // P : double-reverse
            std::reverse(path.begin() + bi,     path.begin() + bj + 1);
            std::reverse(path.begin() + bj + 1, path.begin() + bk + 1);
            break;

        case 1:  // X1 : swap A,B (no reversals)
            std::rotate(path.begin() + bi, path.begin() + bj + 1,
                         path.begin() + bk + 1);
            break;

        case 2:  // X2 : swap, A reversed
        {
            std::reverse(path.begin() + bi, path.begin() + bk + 1);
            std::reverse(path.begin() + bi,
                         path.begin() + bi + (bk - bj));  // restore B orientation
            break;
        }

        case 3:  // X3 : swap, B reversed
            std::reverse(path.begin() + bj + 1, path.begin() + bk + 1); // B^R
            std::rotate(path.begin() + bi, path.begin() + bj + 1,
                         path.begin() + bk + 1);                        // B^R A
            break;
    }
    return true;
}


bool or_opt_improvement(std::vector<int>& path,
                             const std::vector<std::vector<double>>& dist,
                             int maxLen) {
    const int n = static_cast<int>(path.size());
    if (n < 4) return false;                 // need at least 2 interior vertices

    auto d = [&](int a, int b) { return dist[path[a]][path[b]]; };

    double best_gain = 0.0;
    int best_i = -1, best_len = -1, best_pos = -1;   // move [i … i+len-1] after pos

    for (int len = 1; len <= maxLen; ++len)
        for (int i = 1; i + len - 1 <= n - 2; ++i)            // keep endpoints fixed
        {
            int j = i + len - 1;                              // last vertex of the block
            double gain_remove =
                d(i-1, i) + d(j, j+1) - d(i-1, j+1);          // releasing the block

            for (int pos = 0; pos <= n - len - 2; ++pos)      // index in *shrunk* path
            {
                if (pos >= i-1 && pos <= j) continue;         // cannot insert inside

                /* map pos from shrunk path to index in the original path -------- */
                int k = (pos < i) ? pos : pos + len;          // ‘after k’ in old path

                double gain_insert =
                    d(k, k+1) -                                 // edge to be removed
                    (d(k, i) + d(j, k+1));                     // two new edges

                double total_gain = gain_remove + gain_insert;
                if (total_gain > best_gain + 1e-12) {
                    best_gain = total_gain;
                    best_i = i; best_len = len; best_pos = pos;
                }
            }
        }

    if (best_i == -1) return false;           // no improvement found

    /* --------- apply the move ---------------------------------------------- */
    int j = best_i + best_len - 1;
    std::vector<int> block(path.begin() + best_i, path.begin() + j + 1);
    path.erase(path.begin() + best_i, path.begin() + j + 1);

    int insert_pos = (best_i <= best_pos) ? best_pos + 1 - best_len
                                          : best_pos + 1;      // position *after*
    path.insert(path.begin() + insert_pos, block.begin(), block.end());

    return true;
}


/*****Intra-path one-step local search, 2-opt, 3-opt, or-opt*****/
// std::mt19937 rng(std::random_device{}());

std::vector<int> two_opt_one_step(const std::vector<int>& path) {
    if (path.size() <= 3) return path;  // No meaningful swap

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 
    
    size_t a = dist(rng);
    size_t b = dist(rng);
    if (a > b) std::swap(a, b);

    std::reverse(new_path.begin() + a, new_path.begin() + b + 1);
    return new_path;
}


std::vector<int> three_opt_one_step(const std::vector<int>& path) {
    const size_t n = path.size();
    if (n <= 5) return path;               // need at least two interior edges

    /* --- pick three distinct cut points a < b < c (endpoints stay fixed) --- */
    std::uniform_int_distribution<size_t> dist(1, n - 2);
    size_t a, b, c;
    do {
        a = dist(rng);
        b = dist(rng);
        c = dist(rng);
    } while (!(a < b && b < c));

    /* split into four segments S1 | A | B | S4 */
    std::vector<int> S1(path.begin(),         path.begin() + a);
    std::vector<int> A (path.begin() + a,     path.begin() + b);
    std::vector<int> B (path.begin() + b,     path.begin() + c);
    std::vector<int> S4(path.begin() + c,     path.end());

    /* choose one of the four genuine 3-opt reconnections */
    std::uniform_int_distribution<int> case_dist(0, 3);
    switch (case_dist(rng)) {
        case 0: {           // P : double-reverse, keep order A B
            std::reverse(A.begin(), A.end());
            std::reverse(B.begin(), B.end());
            break;
        }
        case 1: {           // X₁ : swap, no reversals
            std::swap(A, B);
            break;
        }
        case 2: {           // X₂ : swap, A reversed
            std::reverse(A.begin(), A.end());
            std::swap(A, B);
            break;
        }
        case 3: {           // X₃ : swap, B reversed
            std::reverse(B.begin(), B.end());
            std::swap(A, B);
            break;
        }
    }

    /* glue everything back together */
    std::vector<int> candidate;
    candidate.reserve(n);
    candidate.insert(candidate.end(), S1.begin(), S1.end());
    candidate.insert(candidate.end(), A.begin(),  A.end());
    candidate.insert(candidate.end(), B.begin(),  B.end());
    candidate.insert(candidate.end(), S4.begin(), S4.end());

    return candidate;
}


// One random “double-bridge” move
std::vector<int> double_bridge(const std::vector<int>& path) {
    const std::size_t n = path.size();
    if (n < 8) return path;                     // need 4 inner edges

    //choose four distinct cut indices in {1 .. n-2} 
    std::uniform_int_distribution<std::size_t> dist(1, n - 2);
    std::size_t a, b, c, d;
    do {
        a = dist(rng);  b = dist(rng);  c = dist(rng);  d = dist(rng);
        std::array<std::size_t,4> cuts = {a,b,c,d};
        std::sort(cuts.begin(), cuts.end());
        a = cuts[0];  b = cuts[1];  c = cuts[2];  d = cuts[3];
    } while (a+1 >= b || b+1 >= c || c+1 >= d); // every segment ≥1

    std::vector<int> new_path;
    new_path.reserve(n);

    new_path.insert(new_path.end(), path.begin(),     path.begin() + a); 
    new_path.insert(new_path.end(), path.begin() + b, path.begin() + c);
    new_path.insert(new_path.end(), path.begin() + a, path.begin() + b); 
    new_path.insert(new_path.end(), path.begin() + c, path.begin() + d); 
    new_path.insert(new_path.end(), path.begin() + d, path.end()); 

    return new_path;
}

//random swap two nodes
std::vector<int> swap_nodes(const std::vector<int>& path) {
    if (path.size() <= 2) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 

    size_t a = dist(rng);
    size_t b = dist(rng);
    std::swap(new_path[a], new_path[b]);

    return new_path;
}

// Relocate the single worst cost node to its cheapest position.
std::vector<int> relocate_node(const std::vector<int>& path,
                                const std::vector<std::vector<double>>& cost) {

    double eps = 1e-12;
                                        
    const std::size_t n = path.size();
    if (n <= 2) return path;

    // find worst cost node
    double best_excess = -std::numeric_limits<double>::infinity();
    std::size_t best_i = n;  

    for (std::size_t i = 1; i + 1 < n; ++i) {
        int u = path[i-1], v = path[i], w = path[i+1];
        double excess = cost[u][v] + cost[v][w] - cost[u][w];
        if (excess > best_excess + eps || (std::abs(excess - best_excess) <= eps && i < best_i)) {
            best_excess = excess;
            best_i = i;
        }
    }

    if (best_i >= n) return path; 

    // remove node
    std::vector<int> new_path = path;
    int v = new_path[best_i];
    int left_old  = new_path[best_i - 1];
    int right_old = new_path[best_i + 1];
    new_path.erase(new_path.begin() + static_cast<long>(best_i));

    // find cheapest re-insertion position (excluding the original slot)
    double best_extra = std::numeric_limits<double>::infinity();
    std::size_t best_pos = 1; // insert after index (i.e., between nodes [pos-1, pos])

    for (std::size_t pos = 1; pos < new_path.size(); ++pos) {
        int x = new_path[pos - 1];
        int y = new_path[pos];

        // // Avoid no-op reinsertion that recreates (left_old, v, right_old)
        // if (x == left_old && y == right_old) continue;

        double extra = cost[x][v] + cost[v][y] - cost[x][y];
        if (extra < best_extra - eps || (std::abs(extra - best_extra) <= eps && pos < best_pos)) {
            best_extra = extra;
            best_pos = pos;
        }
    }

    // apply only if it improves total cost: delta = insert_extra - removal_gain
    double removal_gain = best_excess; // by definition above
    double delta = best_extra - removal_gain;
    if (delta < eps) {
        new_path.insert(new_path.begin() + static_cast<long>(best_pos), v);
        return new_path;
    }
    return path; 
}

// Random select and shift a segment
std::vector<int> shift_segment(const std::vector<int>& path) {
    if (path.size() <= 3) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 

    size_t start = dist(rng);
    size_t end = dist(rng);
    if (start > end) std::swap(start, end);

    std::vector<int> segment(new_path.begin() + start, new_path.begin() + end + 1);
    new_path.erase(new_path.begin() + start, new_path.begin() + end + 1);

    // std::uniform_int_distribution<size_t> insert_dist(1, new_path.size()); 
    std::uniform_int_distribution<size_t> insert_dist(1, new_path.size() - 1);
    size_t insert_pos = insert_dist(rng);
    new_path.insert(new_path.begin() + insert_pos, segment.begin(), segment.end());

    return new_path;
}

// Cut at the two most expensive edges in the current s->t path, 
//and insert the segment in the cheapest position
std::vector<int> relocate_segment(const std::vector<int>& path,
                                        const std::vector<std::vector<double>>& cost,
                                        bool allow_reverse = true,
                                        double eps = 1e-12) {
    const std::size_t n = path.size();
    if (n <= 3) return path; 

    // Find the two most expensive *edges* in the current path.
    std::size_t e1 = 1, e2 = 2; 
    double max1 = -std::numeric_limits<double>::infinity();
    double max2 = -std::numeric_limits<double>::infinity();

    for (std::size_t i = 1; i <= n - 1; ++i) {
        double w = cost[path[i-1]][path[i]];
        if (w > max1) { max2 = max1; e2 = e1; max1 = w; e1 = i; }
        else if (w > max2) { max2 = w; e2 = i; }
    }
    if (e1 > e2) std::swap(e1, e2);   
    std::size_t s = e1;
    std::size_t t = e2 - 1;

    if (s < 1 || t > n - 2 || s > t) return path;

    // extract the segment and the remainder path
    std::vector<int> seg(path.begin() + static_cast<long>(s),
                         path.begin() + static_cast<long>(t) + 1);

    std::vector<int> new_path;
    new_path.reserve(n - (t - s + 1));
    new_path.insert(new_path.end(), path.begin(), path.begin() + static_cast<long>(s));
    new_path.insert(new_path.end(), path.begin() + static_cast<long>(t) + 1, path.end());

    const int L_old = path[s - 1];
    const int R_old = path[t + 1];

    // Precompute costs of the segment in both orientations.
    auto internal_cost = [&](const std::vector<int>& v)->double {
        double s = 0.0;
        for (std::size_t i = 1; i < v.size(); ++i) s += cost[v[i-1]][v[i]];
        return s;
    };
    const double internal = internal_cost(seg);

    std::vector<int> seg_rev = seg;
    std::reverse(seg_rev.begin(), seg_rev.end());
    const double internal_rev = internal_cost(seg_rev);

    // 4) Find the cheapest reinsertion position (and orientation).
    double best_extra = std::numeric_limits<double>::infinity();
    std::size_t best_pos = 1;
    bool best_rev = false;

    for (std::size_t pos = 1; pos < new_path.size(); ++pos) {
        int x = new_path[pos - 1];
        int y = new_path[pos];

        // // Skip exact no-op: same location with same orientation.
        // if (x == L_old && y == R_old) {
        //     // If reversing is allowed, we can still consider the reversed orientation here.
        //     if (allow_reverse) {
        //         double extra_rev = cost[x][seg_rev.front()] + internal_rev + cost[seg_rev.back()][y] - cost[x][y];
        //         if (extra_rev + eps < best_extra) {
        //             best_extra = extra_rev; best_pos = pos; best_rev = true;
        //         }
        //     }
        //     continue;
        // }

        // Same orientation
        double extra = cost[x][seg.front()] + internal + cost[seg.back()][y] - cost[x][y];
        bool use_rev = false;

        // Try reversed orientation if allowed
        if (allow_reverse) {
            double extra_r = cost[x][seg_rev.front()] + internal_rev + cost[seg_rev.back()][y] - cost[x][y];
            if (extra_r < extra) { extra = extra_r; use_rev = true; }
        }

        if (extra + eps < best_extra ||
           (std::abs(extra - best_extra) <= eps && pos < best_pos)) {
            best_extra = extra;
            best_pos = pos;
            best_rev = use_rev;
        }
    }

    // 5) Apply only if total cost improves
    const double removal_gain = cost[L_old][seg.front()] + internal + cost[seg.back()][R_old] - cost[L_old][R_old];
    const double delta = best_extra - removal_gain;

    if (delta < eps) {
        if (best_rev) seg.swap(seg_rev); 
        new_path.insert(new_path.begin() + static_cast<long>(best_pos), seg.begin(), seg.end());
        return new_path;
    }
    return path; 
}

std::vector<int> shuffle_segment(const std::vector<int>& path) {
    if (path.size() <= 3) return path;

    std::vector<int> new_path = path;
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2); 

    size_t start = dist(rng);
    size_t end = dist(rng);
    if (start > end) std::swap(start, end);

    std::shuffle(new_path.begin() + start, new_path.begin() + end + 1, rng);

    return new_path;
}


std::vector<int> change_tile(const std::vector<int>& path,
                             const std::vector<std::vector<double>>& costs,
                             const std::vector<double>& prizes,
                             double budget,
                             std::vector<bool>& visited) {
    if (path.size() <= 2) 
        return path;

    std::vector<int> new_path = path;
    double current_cost = pathCost(costs, path);

    // Remove a random tile 
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2);
    size_t remove_idx = dist(rng);
    int removed_tile = new_path[remove_idx];

    //no longer visited
    visited[removed_tile] = false;

    // Compute cost difference
    double cost_before = costs[new_path[remove_idx - 1]][removed_tile];
    double cost_after = costs[removed_tile][new_path[remove_idx + 1]];
    double new_segment_cost = costs[new_path[remove_idx - 1]][new_path[remove_idx + 1]];
    double freed_budget = cost_before + cost_after - new_segment_cost;

    new_path.erase(new_path.begin() + remove_idx);
    current_cost -= freed_budget;

    insertHighValueCheapest(new_path, costs, prizes, budget, current_cost, visited);

    return new_path;
}

std::vector<int> change_tile(const std::vector<int>& path,
                             const std::vector<std::vector<double>>& costs,
                             const std::vector<double>& prizes,
                             double budget) {
    if (path.size() <= 2) 
        return path;

    std::vector<bool> visited(costs.size(), false);
    for(auto p : path) visited[p] = true;

    std::vector<int> new_path = path;
    double current_cost = pathCost(costs, path);

    // Remove a random tile 
    std::uniform_int_distribution<size_t> dist(1, path.size() - 2);
    size_t remove_idx = dist(rng);
    int removed_tile = new_path[remove_idx];

    //no longer visited
    // visited[removed_tile] = false;

    // Compute cost difference
    double cost_before = costs[new_path[remove_idx - 1]][removed_tile];
    double cost_after = costs[removed_tile][new_path[remove_idx + 1]];
    double new_segment_cost = costs[new_path[remove_idx - 1]][new_path[remove_idx + 1]];
    double freed_budget = cost_before + cost_after - new_segment_cost;

    new_path.erase(new_path.begin() + remove_idx);
    current_cost -= freed_budget;

    insertHighValueCheapest(new_path, costs, prizes, budget, current_cost, visited);

    return new_path;
}

std::vector<int> swap_by_ratio(std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes,
              double budget, std::vector<bool>& visited) {
    const std::size_t n = path.size();
    if (n <= 3) return path;   

    auto cost_edge = [&](int u,int v){ return costs[u][v]; };

    /* current tour cost */
    double current_cost = 0.0;
    for (std::size_t i = 1; i < n; ++i)
        current_cost += cost_edge(path[i-1], path[i]);

    //1. find worst prize/delta‑cost node
    int worst_idx = -1;
    double worst_ratio = std::numeric_limits<double>::infinity();

    for (std::size_t i = 1; i < n-1; ++i) {
        int u = path[i-1], v = path[i], w = path[i+1];
        double delta = cost_edge(u,v) + cost_edge(v,w) - cost_edge(u,w);
        double ratio = prizes[v] / std::max(delta, 1e-12);
        if (ratio < worst_ratio) {
            worst_ratio = ratio;
            worst_idx   = static_cast<int>(i);
        }
    }
    if (worst_idx == -1) return path;    

    //remove it and update cost & visited
    int removed = path[worst_idx];
    double gain =
        cost_edge(path[worst_idx-1], removed) +
        cost_edge(removed, path[worst_idx+1]) -
        cost_edge(path[worst_idx-1], path[worst_idx+1]);

    current_cost -= gain;
    path.erase(path.begin()+worst_idx);
    visited[removed] = false;

    //evaluate each unvisited node for best ratio insert
    int best_city = -1;
    double best_ratio = -1.0;
    std::size_t best_pos = 0;
    double best_extra = 0.0;

    for (int city = 0; city < static_cast<int>(costs.size()); ++city)
        if (!visited[city] && city != path.front() && city != path.back())
        {
            //find cheapest edge to insert
            double cheapest_extra = std::numeric_limits<double>::infinity();
            std::size_t cheapest_pos = 0;

            for (std::size_t j = 0; j < path.size()-1; ++j) {
                double extra = cost_edge(path[j], city) +
                               cost_edge(city,    path[j+1]) -
                               cost_edge(path[j], path[j+1]);
                if (extra < cheapest_extra) {
                    cheapest_extra = extra;
                    cheapest_pos   = j+1;
                }
            }
            /* append just before t */
            double extra_tail = cost_edge(path.back(), city);
            if (extra_tail < cheapest_extra) {
                cheapest_extra = extra_tail;
                cheapest_pos   = path.size();              // insert before t
            }

            /* does it fit?  evaluate ratio */
            if (current_cost + cheapest_extra <= budget + 1e-12) {
                double ratio = prizes[city] / std::max(cheapest_extra, 1e-12);
                if (ratio > best_ratio) {
                    best_ratio = ratio;
                    best_city  = city;
                    best_pos   = cheapest_pos;
                    best_extra = cheapest_extra;
                }
            }
        }

    //insert the chosen city if any fits
    if (best_city != -1) {
        path.insert(path.begin()+best_pos, best_city);
        visited[best_city] = true;
        // current_cost += best_extra;
    }
    return path; 
}


std::vector<int> swap_by_ratio(std::vector<int> path,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes, double budget) {
    const std::size_t n = path.size();
    if (n <= 3) return path;   

    std::vector<bool> visited(costs.size(), false);
    for(auto p : path) visited[p] = true;

    auto cost_edge = [&](int u,int v){ return costs[u][v]; };

    // current tour cost
    double current_cost = 0.0;
    for (std::size_t i = 1; i < n; ++i)
        current_cost += cost_edge(path[i-1], path[i]);

    //find worst prize/delta‑cost node
    int worst_idx = -1;
    double worst_ratio = std::numeric_limits<double>::infinity();

    for (std::size_t i = 1; i < n-1; ++i) {
        int u = path[i-1], v = path[i], w = path[i+1];
        double delta = cost_edge(u,v) + cost_edge(v,w) - cost_edge(u,w);
        double ratio = prizes[v] / std::max(delta, 1e-12);
        if (ratio < worst_ratio) {
            worst_ratio = ratio;
            worst_idx   = static_cast<int>(i);
        }
    }
    if (worst_idx == -1) return path; 

    //remove it and update cost & visited
    int removed = path[worst_idx];
    double gain =
        cost_edge(path[worst_idx-1], removed) +
        cost_edge(removed, path[worst_idx+1]) -
        cost_edge(path[worst_idx-1], path[worst_idx+1]);

    current_cost -= gain;
    path.erase(path.begin()+worst_idx);
    visited[removed] = false;

    //evaluate each unvisited node for best ratio insert
    int best_city = -1;
    double best_ratio = -1.0;
    std::size_t best_pos = 0;
    double best_extra = 0.0;

    for (int city = 0; city < static_cast<int>(costs.size()); ++city)
        if (!visited[city] && city != path.front() && city != path.back())
        {
            //find cheapest insertion edge
            double cheapest_extra = std::numeric_limits<double>::infinity();
            std::size_t cheapest_pos = 0;

            for (std::size_t j = 0; j < path.size()-1; ++j) {
                double extra = cost_edge(path[j], city) +
                               cost_edge(city,    path[j+1]) -
                               cost_edge(path[j], path[j+1]);
                if (extra < cheapest_extra) {
                    cheapest_extra = extra;
                    cheapest_pos   = j+1;
                }
            }
            //append just before t
            double extra_tail = cost_edge(path.back(), city);
            if (extra_tail < cheapest_extra) {
                cheapest_extra = extra_tail;
                cheapest_pos   = path.size();
            }

            //evaluate ratio 
            if (current_cost + cheapest_extra <= budget + 1e-12) {
                double ratio = prizes[city] / std::max(cheapest_extra, 1e-12);
                if (ratio > best_ratio) {
                    best_ratio = ratio;
                    best_city  = city;
                    best_pos   = cheapest_pos;
                    best_extra = cheapest_extra;
                }
            }
        }

    //insert the chosen city if any fits
    if (best_city != -1) {
        path.insert(path.begin()+best_pos, best_city);
        visited[best_city] = true;
        // current_cost += best_extra;
    }

    return path; 
}

//********************** */
std::vector<int> change_tile(
    const std::vector<int>& path,
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    std::vector<bool>& visited ) {
    if (path.size() <= 2) {
        visited.assign(costs.size(), false);
        for (int node : path) visited[node] = true;
        return path;
    }

    std::vector<int> candidate = path;
    std::uniform_int_distribution<std::size_t> pick(
        1, candidate.size() - 2);
    candidate.erase(candidate.begin() +
                    static_cast<std::ptrdiff_t>(pick(rng)));
    repairPath(
        costs, pixel_probs, member_pixels, candidate, budget, visited);
    return candidate;
}

std::vector<int> change_tile(
    const std::vector<int>& path,
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget) {
    std::vector<bool> visited;
    return change_tile(
        path, costs, pixel_probs, member_pixels, budget, visited);
}

std::vector<int> swap_by_ratio(
    std::vector<int> path,
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget,
    std::vector<bool>& visited ) {
    if (path.size() <= 3) {
        repairPath(
            costs, pixel_probs, member_pixels, path, budget, visited);
        return path;
    }

    const auto counts = pixelCoverCounts(
        path, member_pixels, pixel_probs.size());
    constexpr double tolerance = 1e-12;
    std::size_t worst_index = path.size();
    double worst_ratio = std::numeric_limits<double>::infinity();
    double best_saving = -std::numeric_limits<double>::infinity();

    for (std::size_t index = 1; index + 1 < path.size(); ++index) {
        const double saving = removalDeltaCost(costs, path, index);
        if (saving <= tolerance) {
            continue;
        }
        const double loss = exclusivePixelPrize(
            path[index], pixel_probs, member_pixels, counts);
        const double ratio = loss / saving;
        if (ratio < worst_ratio - tolerance ||
            (std::abs(ratio - worst_ratio) <= tolerance &&
             saving > best_saving + tolerance)) {
            worst_ratio = ratio;
            best_saving = saving;
            worst_index = index;
        }
    }

    if (worst_index < path.size()) {
        path.erase(path.begin() +
                   static_cast<std::ptrdiff_t>(worst_index));
    }
    repairPath(
        costs, pixel_probs, member_pixels, path, budget, visited);
    return path;
}

std::vector<int> swap_by_ratio(
    std::vector<int> path,
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    double budget ) {
    std::vector<bool> visited;
    return swap_by_ratio(
        std::move(path), costs, pixel_probs, member_pixels, budget, visited);
}

/************************* */



// // *********interpath local serch ******

// 2-Opt per route
bool twoOptAllRoutes(const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<int>>& paths,
                    std::vector<double>& pathCosts,
                    double budgetT) {

    bool any = false;
    for (int r = 0; r < (int)paths.size(); ++r) {
        auto& route = paths[r];
        int n = (int)route.size();
        if (n <= 4) continue;

        bool improved = true;
        while (improved) {
            improved = false;
            double bestDelta = -1e-12; int bi=-1, bk=-1;
            const double oldT = pathCosts[r];

            for (int i = 1; i <= n-3; ++i) {
                for (int k = i+1; k <= n-2; ++k) {
                    std::reverse(route.begin()+i, route.begin()+k+1);
                    double newT = pathCost(costs, route);
                    std::reverse(route.begin()+i, route.begin()+k+1);
                    double d = newT - oldT;
                    if (d < bestDelta && newT <= budgetT + 1e-12) {
                        bestDelta = d; bi=i; bk=k;
                    }
                }
            }
            if (bi != -1) {
                std::reverse(route.begin()+bi, route.begin()+bk+1);
                pathCosts[r] = pathCost(costs, route);
                improved = true; any = true;
            }
        }
    }
    return any;
}

// Exchange nodes between tours (cost-reducing only)
bool swapNodesBetweenPaths(const std::vector<std::vector<double>>& costs,
                           std::vector<std::vector<int>>& paths,
                           std::vector<double>& pathCosts,
                           double budget) { 
    const int R = (int)paths.size();
    bool improved = false;

    double bestDeltaSum = 0.0; // positive means total reduction
    int br1=-1, br2=-1, bi=-1, bj=-1;
    double bestDelta1=0.0, bestDelta2=0.0; // store deltas to update pathCosts correctly

    for (int r1 = 0; r1 < R; ++r1) {
        if ((int)paths[r1].size() <= 2) continue;
        for (int r2 = r1 + 1; r2 < R; ++r2) {
            if ((int)paths[r2].size() <= 2) continue;

            for (int i = 1; i < (int)paths[r1].size() - 1; ++i) {
                for (int j = 1; j < (int)paths[r2].size() - 1; ++j) {
                    int a1 = paths[r1][i-1], x1 = paths[r1][i], b1 = paths[r1][i+1];
                    int a2 = paths[r2][j-1], x2 = paths[r2][j], b2 = paths[r2][j+1];

                    // delta for route r1 if we replace x1 by x2 at position i
                    double old1 = costs[a1][x1] + costs[x1][b1];
                    double new1 = costs[a1][x2] + costs[x2][b1];
                    double d1 = new1 - old1;

                    // delta for route r2 if we replace x2 by x1 at position j
                    double old2 = costs[a2][x2] + costs[x2][b2];
                    double new2 = costs[a2][x1] + costs[x1][b2];
                    double d2 = new2 - old2;

                    double newCost1 = pathCosts[r1] + d1;
                    double newCost2 = pathCosts[r2] + d2;

                    // cost-reducing swap only
                    double totalReduction = -(d1 + d2);
                    if (totalReduction > bestDeltaSum &&
                        newCost1 <= budget + 1e-12 &&
                        newCost2 <= budget + 1e-12) {
                        bestDeltaSum = totalReduction;
                        br1 = r1; br2 = r2; bi = i; bj = j;
                        bestDelta1 = d1; bestDelta2 = d2;
                    }
                }
            }
        }
    }

    if (bestDeltaSum > 0.0) {
        std::swap(paths[br1][bi], paths[br2][bj]);
        // Use the deltas found during the search (do NOT recompute after swap)
        pathCosts[br1] += bestDelta1;
        pathCosts[br2] += bestDelta2;
        improved = true;
    }
    return improved;
}



bool replaceNodesInPaths(const std::vector<std::vector<double>>& costs,
                                   std::vector<std::vector<int>>& paths,
                                   std::vector<double>& pathCosts,
                                   const std::vector<double>& prizes,
                                   double budget) {
    int n = costs.size();
    bool improved = false;
    bool localImproved = true;
    
    std::vector<bool> visited = getVisited(paths, n);
    
    double bestScoreIncrease = 0.0;
    int bestNode = -1;
    int bestRemoveRoute = -1, bestRemovePos = -1;
    int bestInsertRoute = -1, bestInsertPos = -1;

    for (int node = 0; node < n; ++node) {
        if(visited[node]) continue;
        // Try replacing existing nodes with lower scores
        for (int r = 0; r < paths.size(); r++) {
            if (paths[r].size() <= 2) continue;
            
            for (int pos = 1; pos < paths[r].size() - 1; pos++) {
                int currentNode = paths[r][pos];
                
                if (prizes[currentNode] < prizes[node]) {
                    // Calculate removal cost
                    double removalCost = costs[paths[r][pos-1]][currentNode] + 
                                        costs[currentNode][paths[r][pos+1]] - 
                                        costs[paths[r][pos-1]][paths[r][pos+1]];
                    
                    // Try inserting the new node in the same position
                    double insertCost = costs[paths[r][pos-1]][node] + 
                                        costs[node][paths[r][pos+1]] - 
                                        costs[paths[r][pos-1]][paths[r][pos+1]];
                    
                    double netCostChange = insertCost - removalCost;
                    double newPathCost = pathCosts[r] + netCostChange;
                    
                    if (newPathCost <= budget) {
                        double scoreIncrease = prizes[node] - prizes[currentNode];
                        if (scoreIncrease > bestScoreIncrease) {
                            bestScoreIncrease = scoreIncrease;
                            bestNode = node;
                            bestRemoveRoute = r;
                            bestRemovePos = pos;
                            bestInsertRoute = r;
                            bestInsertPos = pos; 
                        }
                    }
                }
            }
        }
    }
    
    if (bestScoreIncrease > 0) {
        // Simple replacement in same position
        int oldNode = paths[bestRemoveRoute][bestRemovePos];
        double oldCost = costs[paths[bestRemoveRoute][bestRemovePos-1]][oldNode] + 
                        costs[oldNode][paths[bestRemoveRoute][bestRemovePos+1]] - 
                        costs[paths[bestRemoveRoute][bestRemovePos-1]][paths[bestRemoveRoute][bestRemovePos+1]];
        
        paths[bestRemoveRoute][bestRemovePos] = bestNode;
        
        double newCost = costs[paths[bestRemoveRoute][bestRemovePos-1]][bestNode] + 
                        costs[bestNode][paths[bestRemoveRoute][bestRemovePos+1]] - 
                        costs[paths[bestRemoveRoute][bestRemovePos-1]][paths[bestRemoveRoute][bestRemovePos+1]];
        
        pathCosts[bestRemoveRoute] += (newCost - oldCost);
        
        // Update visited state
        visited[oldNode] = false;
        visited[bestNode] = true;
        improved = true;
    }
    return improved;
}


// Insert neighborhood - increase score by inserting feasible locations
bool insertNodesToPaths(const std::vector<std::vector<double>>& costs,
                         std::vector<std::vector<int>>& paths,
                         std::vector<double>& pathCosts,
                         const std::vector<double>& prizes,
                         double budget) {
    int n = costs.size();
    bool improved = false;
    
    // Build initial visited array and unvisited list once
    std::vector<bool> visited = getVisited(paths, n);

    double bestHeuristic = 0.0;
    int bestNode = -1;
    int bestRoute = -1, bestPos = -1;
    double bestInsertCost = 0.0;
    
    // For every non-included location
    for (int node = 0; node < n; node++) {
        if(visited[node]) continue;
        auto [insertRoute, insertPos, insertCost] = 
            cheapestInsertionAll(costs, paths, pathCosts, prizes, budget, node);
        
        if (insertRoute != -1) {
            double heuristic = prizes[node] / std::max(insertCost, 1e-12);  // Prize to cost ratio
            
            if (heuristic > bestHeuristic) {
                bestHeuristic = heuristic;
                bestNode = node;
                bestRoute = insertRoute;
                bestPos = insertPos;
                bestInsertCost = insertCost;
            }
        }
    }
    
    if (bestHeuristic > 0) {
        // Insert the node and update path cost
        paths[bestRoute].insert(paths[bestRoute].begin() + bestPos, bestNode);
        pathCosts[bestRoute] += bestInsertCost;
        
        improved = true;
    }
    
    return improved;
}



//*****Inter-path search

// Relocate a contiguous block of k internal nodes from one route to another.
// Accepts iff total prize strictly increases (uses multiPathsSumPrize).
bool relocate_k(std::vector<std::vector<int>>& paths,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                int k) {
    if (paths.empty() || k <= 0) return false;

    const double basePrize = multiPathsSumPrize(prizes, paths);

    for (size_t from = 0; from < paths.size(); ++from) {
        // Need at least [start, k nodes, t]
        if (paths[from].size() <= static_cast<size_t>(k) + 2) continue;

        for (size_t to = 0; to < paths.size(); ++to) {
            if (from == to) continue;
            if (paths[to].size() < 2) continue; // need [start, t] to insert between

            // Valid cut positions are internal and must leave at least k nodes before the terminal:
            // cut ∈ [1, paths[from].size() - k - 1]
            const size_t lo = 1;
            const size_t hi = paths[from].size() - k - 1;
            if (hi < lo) continue;

            std::uniform_int_distribution<size_t> pick(lo, hi);
            size_t cut = pick(rng);

            auto first = paths[from].begin() + static_cast<std::ptrdiff_t>(cut);
            auto last  = first + k;
            std::vector<int> block(first, last);

            // Remove block from 'from'
            paths[from].erase(first, last);

            // Find cheapest insertion position in 'to' (between any consecutive pair)
            double best_extra = std::numeric_limits<double>::infinity();
            size_t best_pos = 1; // insert before index best_pos
            for (size_t j = 0; j + 1 < paths[to].size(); ++j) {
                int a = paths[to][j];
                int b = paths[to][j + 1];
                double extra = costs[a][block.front()] +
                               costs[block.back()][b] -
                               costs[a][b];
                if (extra < best_extra) { best_extra = extra; best_pos = j + 1; }
            }

            // Insert block into 'to'
            paths[to].insert(paths[to].begin() + static_cast<std::ptrdiff_t>(best_pos),
                             block.begin(), block.end());

            // Accept iff prize improves
            const double newPrize = multiPathsSumPrize(prizes, paths);
            if (newPrize > basePrize + 1e-12) {
                return true;
            }

            // Roll back and try next
            paths[to].erase(paths[to].begin() + static_cast<std::ptrdiff_t>(best_pos),
                            paths[to].begin() + static_cast<std::ptrdiff_t>(best_pos + k));
            paths[from].insert(paths[from].begin() + static_cast<std::ptrdiff_t>(cut),
                               block.begin(), block.end());
        }
    }
    return false;
}


// tile swap between two paths (internal nodes only)
// Accept if total prize increases, or (optionally) if prize ties but total travel decreases.
bool swap_1_1(std::vector<std::vector<int>>& paths,
              const std::vector<std::vector<double>>& costs,
              const std::vector<double>& prizes,
              bool accept_if_cost_reduced) {
    if (paths.size() < 2) return false;

    const double basePrize = multiPathsSumPrize(prizes, paths);
    const double baseCost  = multiPathsMaxCost(costs, paths); // used only for tie-break

    for (size_t p = 0; p < paths.size(); ++p) {
        if (paths[p].size() < 3) continue; // need [start, v, t]
        for (size_t q = p + 1; q < paths.size(); ++q) {
            if (paths[q].size() < 3) continue;

            for (size_t i = 1; i + 1 < paths[p].size(); ++i) {
                for (size_t j = 1; j + 1 < paths[q].size(); ++j) {
                    std::swap(paths[p][i], paths[q][j]);

                    const double newPrize = multiPathsSumPrize(prizes, paths);
                    if (newPrize > basePrize + 1e-12) {
                        return true;
                    }
                    if (accept_if_cost_reduced && std::abs(newPrize - basePrize) <= 1e-12) {
                        const double newCost = multiPathsMaxCost(costs, paths);
                        if (newCost < baseCost - 1e-12) {
                            return true;
                        }
                    }

                    // rollback
                    std::swap(paths[p][i], paths[q][j]);
                }
            }
        }
    }
    return false;
}


/* ---------------- Cross-exchange / 2-opt* (tail-safe, O(1) Δ) ------------- */
bool cross_exchange(std::vector<std::vector<int>>& paths,
                    const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes) {
    if (paths.size() < 2) return false;

    constexpr double EPS = 1e-12;

    for (size_t p = 0; p < paths.size(); ++p) {
        auto &P = paths[p];
        if (P.size() < 4) continue;              // need [s, x, y, t]
        for (size_t q = p + 1; q < paths.size(); ++q) {
            auto &Q = paths[q];
            if (Q.size() < 4) continue;

            for (size_t i = 1; i + 1 < P.size(); ++i) {
                int aP = P[i - 1];  
                int uP = P[i];

                for (size_t j = 1; j + 1 < Q.size(); ++j) {
                    int aQ = Q[j - 1]; 
                    int uQ = Q[j];

                    double delta = (costs[aP][uQ] - costs[aP][uP]) + (costs[aQ][uP] - costs[aQ][uQ]);

                    if (delta < -EPS) {
                        // Perform the swap: exchange interior tails
                        std::vector<int> tailP(P.begin() + static_cast<std::ptrdiff_t>(i), P.end() - 1);
                        std::vector<int> tailQ(Q.begin() + static_cast<std::ptrdiff_t>(j), Q.end() - 1);

                        P.erase(P.begin() + static_cast<std::ptrdiff_t>(i), P.end() - 1);
                        Q.erase(Q.begin() + static_cast<std::ptrdiff_t>(j), Q.end() - 1);

                        P.insert(P.end() - 1, tailQ.begin(), tailQ.end());
                        Q.insert(Q.end() - 1, tailP.begin(), tailP.end());
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

// Random-pair merge -> greedy split under per-route budget B.
// bool merge_split_greedy_budget(std::vector<std::vector<int>>& paths,
//                                 const std::vector<std::vector<double>>& costs,
//                                 const std::vector<double>& prizes,
//                                 double budget) {
//     if (paths.size() <= 2) return false;

//     // Pick a random pair (a, b), a != b
//     int m = (int)paths.size();
//     std::uniform_int_distribution<size_t> pickA(0, m-1);
//     size_t a = pickA(rng);
//     std::uniform_int_distribution<size_t> pickB(0, m-1);
//     size_t b = pickB(rng);
//     if (b == a) b = (b + 1) % m;

//     const auto& PA = paths[a];
//     const auto& PB = paths[b];
//     if (PA.size() <= 2 || PB.size() <= 2) return false; 

//     // Baselines
//     const double basePrize = multiPathsSumPrize(prizes, paths);

//     // Build merged internal sequence: internals of A, then internals of B
//     std::vector<int> merged;
//     merged.reserve((PA.size() > 1 ? PA.size()-2 : 0) + (PB.size() > 1 ? PB.size()-2 : 0));
//     for (size_t i = 0; i + 1 < PA.size()-1; ++i) {
//         int v = PA[i]; 
//         merged.push_back(v);
//     }
//     for (size_t i = 1; i + 1 < PB.size(); ++i) {
//         int v = PB[i];
//         merged.push_back(v);
//     }

//     while(two_opt_st(merged, costs)){};

//     int t = PA.back();
//     double cur_cost = 0;
//     int i = 1;
//     for(; i < merged.size(); ++i) {
//         if(cur_cost + costs[merged[i-1]][merged[i]] + costs[merged[i]][t] > budget)
//             break;
//     }
//     std::vector<int> Anew(merged.begin(), merged.begin() + i);
//     Anew.push_back(t);
//     std::vector<int> Bnew{PB.front()};
//     Bnew.insert(Bnew.end(), merged.begin()+i, merged.end());


//     paths[a] = std::move(Anew);
//     paths[b] = std::move(Bnew);

//     double newPrize = multiPathsSumPrize(prizes, paths);
//     if (newPrize > basePrize + 1e-12) {
//         return true;
//     }
//     return false;
// }

bool merge_split_greedy_budget(std::vector<std::vector<int>>& paths,
                                const std::vector<std::vector<double>>& costs,
                                const std::vector<double>& prizes,
                                double budget) {
    if (paths.size() <= 2) return false;

    int m = (int)paths.size();
    std::uniform_int_distribution<size_t> pickA(0, m-1);
    size_t a = pickA(rng);
    std::uniform_int_distribution<size_t> pickB(0, m-1);
    size_t b = pickB(rng);
    if (b == a) b = (b + 1) % m;

    const auto& PA = paths[a];
    const auto& PB = paths[b];
    if (PA.size() <= 2 || PB.size() <= 2) return false;

    const double basePrize = multiPathsSumPrize(prizes, paths);

    // Merge interior nodes only (skip first and last = depots)
    std::vector<int> merged;
    merged.reserve((PA.size() - 2) + (PB.size() - 2));
    for (size_t i = 1; i + 1 < PA.size(); ++i)
        merged.push_back(PA[i]);
    for (size_t i = 1; i + 1 < PB.size(); ++i)
        merged.push_back(PB[i]);

    if (merged.empty()) return false;

    while (two_opt_st(merged, costs)) {};

    // Split: greedily assign to A under budget, rest to B
    int startA = PA.front(), endA = PA.back();
    int startB = PB.front(), endB = PB.back();

    std::vector<int> Anew{startA};
    double cur_cost = 0;
    size_t split = 0;
    for (size_t i = 0; i < merged.size(); ++i) {
        int prev = Anew.back();
        double add_cost = costs[prev][merged[i]] + costs[merged[i]][endA]
                        - costs[prev][endA];
        if (cur_cost + add_cost > budget) { split = i; break; }
        cur_cost += costs[prev][merged[i]];
        Anew.push_back(merged[i]);
        split = i + 1;
    }
    Anew.push_back(endA);

    std::vector<int> Bnew{startB};
    Bnew.insert(Bnew.end(), merged.begin() + split, merged.end());
    Bnew.push_back(endB);

    // Save originals for rollback
    auto oldA = paths[a];
    auto oldB = paths[b];

    paths[a] = std::move(Anew);
    paths[b] = std::move(Bnew);

    double newPrize = multiPathsSumPrize(prizes, paths);
    if (newPrize > basePrize + 1e-12) {
        return true;
    }

    // Rollback
    paths[a] = std::move(oldA);
    paths[b] = std::move(oldB);
    return false;
}

bool ejection_chain(std::vector<std::vector<int>>& paths,
                    const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes,
                    int depth) {
    if (paths.empty() || depth <= 0) return false;

    // Choose only among routes that have an internal node to remove.
    std::vector<size_t> candidates;
    candidates.reserve(paths.size());
    for (size_t r = 0; r < paths.size(); ++r) {
        if (paths[r].size() > 3) candidates.push_back(r); // need [s, v, ..., t]
    }
    if (candidates.empty()) return false;

    const auto original = paths;
    const double basePrize = multiPathsSumPrize(prizes, paths);

    // Pick a random route and a random internal position
    std::uniform_int_distribution<size_t> pick_route(0, candidates.size() - 1);
    size_t p = candidates[pick_route(rng)];

    std::uniform_int_distribution<size_t> pick_pos(1, paths[p].size() - 2);
    size_t pos = pick_pos(rng);

    int moving = paths[p][pos];
    paths[p].erase(paths[p].begin() + static_cast<std::ptrdiff_t>(pos));

    for (int step = 0; step < depth; ++step) {
        // (1) Best insertion position over all paths/gaps by extra travel
        size_t best_i = SIZE_MAX, best_pos = 1;
        double best_extra = std::numeric_limits<double>::infinity();

        for (size_t i = 0; i < paths.size(); ++i) {
            const auto& R = paths[i];
            if (R.size() < 2) continue;            // need at least [s, t]
            for (size_t j = 0; j + 1 < R.size(); ++j) {
                double extra = costs[R[j]][moving] + costs[moving][R[j+1]] - costs[R[j]][R[j+1]];
                if (extra < best_extra) { best_extra = extra; best_i = i; best_pos = j + 1; }
            }
        }
        if (best_i == SIZE_MAX) { paths = original; return false; }

        paths[best_i].insert(paths[best_i].begin() + static_cast<std::ptrdiff_t>(best_pos), moving);

        // Stop if we reached max depth or no internal nodes to kick
        if (step == depth - 1 || paths[best_i].size() <= 3) break;

        // (2) Choose a "kick" node using removalDeltaCost:
        //     prefer large time saved per unit prize -> (saved_time / max(ε, prize))
        int kick_idx = -1;
        double best_score = -1.0;
        for (size_t j = 1; j + 1 < paths[best_i].size(); ++j) {
            int u = paths[best_i][j];
            double saved = removalDeltaCost(costs, paths[best_i], j); // costs[a][u] + costs[u][b] - costs[a][b]
            double denom = std::max(1e-12, prizes[u]);
            double score = saved / denom;  // higher is better: big time save, small prize loss
            if (score > best_score) { best_score = score; kick_idx = static_cast<int>(j); }
        }
        if (kick_idx < 0) break; // nothing to eject (shouldn't happen if size>3)

        moving = paths[best_i][kick_idx];
        paths[best_i].erase(paths[best_i].begin() + kick_idx);
    }

    // Accept iff total prize increases
    const double newPrize = multiPathsSumPrize(prizes, paths);
    if (newPrize > basePrize + 1e-12) {
        return true;
    }
    paths = original;
    return false;
}



