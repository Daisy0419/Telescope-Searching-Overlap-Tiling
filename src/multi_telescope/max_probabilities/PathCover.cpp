#include "GCP.h"
#include "PathOperations.h"
#include "LocalSearch.h"
#include "PathCover.h"

#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <unordered_set>
#include <future>


// Components of G(J, λ): keep edges (u,v) with w(u,v) ≤ λ (among customers only)
std::vector<std::vector<int>> lambda_components(const std::vector<std::vector<double>>& costs,
                                                const std::vector<int>& customers,
                                                double lambda) {
    int m = (int)customers.size();
    if (m == 0) return {};

    std::vector<std::vector<int>> adj(m);
    for (int i = 0; i < m; ++i) {
        int u = customers[i];
        for (int j = i + 1; j < m; ++j) {
            int v = customers[j];
            if (costs[u][v] <= lambda) { adj[i].push_back(j); adj[j].push_back(i); }
        }
    }
    std::vector<int> vis(m, 0);
    std::vector<std::vector<int>> comps;
    for (int i = 0; i < m; ++i) if (!vis[i]) {
        std::vector<int> comp;
        std::queue<int> q; q.push(i); vis[i] = 1;
        while (!q.empty()) {
            int a = q.front(); q.pop();
            comp.push_back(customers[a]);
            for (int b : adj[a]) if (!vis[b]) { vis[b] = 1; q.push(b); }
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}


// // Greedily split a customer-only cycle into contiguous segments.
// // Each segment is attached to the closest available start and then to a common end t,
// // enforcing: cost(start->first) + segment_edge_sum + cost(last->t) <= B.
// // Uses at most 'remaining_num_path' routes and requires that to match the # of unvisited starts.
// // Returns true on success and fills routes_out with [start, ..., t] per segment; false otherwise.
// bool split_cycle(const std::vector<std::vector<double>>& costs, 
//                         const std::vector<int>& cyc,
//                         const std::vector<int>& starts,
//                         std::vector<int>& visited_starts,
//                         int& remaining_num_path, 
//                         double B,
//                         int t_node,
//                         std::vector<std::vector<int>>& routes_out) {

//     // ---- validation ----
//     if (B <= 0.0) {                    // <-- fix braces
//         std::cout << "Budget B must be positive.";
//         return false;
//     }

//     if (cyc.size() == 0) {             // nothing to cover
//         routes_out.clear();
//         return true;
//     }

//     // const int m = static_cast<int>(cyc.size()) - 1; // number of edges

//     // ---- available starts = starts \ visited_starts ----
//     std::unordered_set<int> used(visited_starts.begin(), visited_starts.end());
//     std::vector<int> avail_starts;
//     avail_starts.reserve(starts.size());
//     for (int s : starts) if (!used.count(s)) avail_starts.push_back(s);

//     if (remaining_num_path < 0) return false;
//     if (static_cast<int>(avail_starts.size()) != remaining_num_path) {
//         // Must match by spec
//         return false;
//     }

//     auto take_closest_start = [&](int v_first, int& idx, double& c_sf) -> bool {
//         idx = -1; c_sf = std::numeric_limits<double>::infinity();
//         for (int i = 0; i < static_cast<int>(avail_starts.size()); ++i) {
//             const int s = avail_starts[i];
//             const double c = costs[s][v_first];
//             if (c < c_sf) { c_sf = c; idx = i; }
//         }
//         return idx >= 0 && std::isfinite(c_sf);
//     };

//     routes_out.clear();
//     routes_out.reserve(remaining_num_path);

//     int start_idx = 0; // index of first vertex of the next (uncovered) segment in 'cyc'
//     while (start_idx < m) {
//         if (remaining_num_path == 0 || avail_starts.empty()) {
//             routes_out.clear();
//             return false; // no starts left, yet cycle not fully covered
//         }

//         const int v_first = cyc[start_idx];

//         // Choose closest available start for this segment
//         int best_start_pos = -1;
//         double c_start_to_first = 0.0;
//         if (!take_closest_start(v_first, best_start_pos, c_start_to_first)) {
//             routes_out.clear();
//             return false;
//         }

//         // If even a single-vertex segment can't fit (start->v_first + v_first->t), fail
//         const double c_first_to_t = costs[v_first][t_node];
//         if (c_start_to_first + c_first_to_t > B) {
//             routes_out.clear();
//             return false;
//         }

//         // Greedily extend segment while keeping:
//         // c_start_to_first + sum(segment edges) + cost(current_end -> t) <= B
//         int end_idx = start_idx;            // will include cyc[end_idx]
//         double acc_edges = 0.0;
//         double c_end_to_t = c_first_to_t;   // cost(last -> t) for current end

//         while (end_idx + 1 <= m) {
//             const double next_edge = costs[cyc[end_idx]][cyc[end_idx + 1]];
//             const double new_c_end_to_t = costs[cyc[end_idx + 1]][t_node];
//             if (c_start_to_first + (acc_edges + next_edge) + new_c_end_to_t <= B) {
//                 acc_edges += next_edge;
//                 end_idx += 1;
//                 c_end_to_t = new_c_end_to_t;
//             } else {
//                 break; // can't extend further within budget
//             }
//         }

//         std::vector<int> route;
//         route.reserve(2 + (end_idx - start_idx + 1));
//         const int chosen_start = avail_starts[best_start_pos];  // <-- capture
//         route.push_back(chosen_start);                          // chosen start
//         for (int k = start_idx; k <= end_idx; ++k) route.push_back(cyc[k]);
//         route.push_back(t_node);                                // common end
//         routes_out.push_back(std::move(route));

//         // consume this start
//         visited_starts.push_back(chosen_start);                 // <-- NEW: record as visited
//         avail_starts.erase(avail_starts.begin() + best_start_pos);
//         --remaining_num_path;

//         // next uncovered vertex (AFTER end_idx)
//         start_idx = end_idx + 1;
//     }

//     return true; // covered entire cycle
// }
bool split_cycle(const std::vector<std::vector<double>>& costs, 
                 const std::vector<int>& cyc,
                 const std::vector<int>& starts,
                 std::vector<int>& visited_starts,
                 int& remaining_num_path, 
                 double B,
                 int t_node,
                 std::vector<std::vector<int>>& routes_out) {

    // ---- validation ----
    if (B <= 0.0) {                    // <-- fix braces
        std::cout << "Budget B must be positive.";
        return false;
    }
    const int N = (int)costs.size();
    if (cyc.size() == 0) {             // nothing to cover
        routes_out.clear();
        return true;
    }
    // if ((int)costs.size() == 0 || (int)costs[0].size() != N) return false;
    // if (t_node < 0 || t_node >= N) return false;

    // Determine whether cyc is closed (last==first). If closed, ignore the
    // duplicated last vertex when iterating; we still DO NOT add the wrap edge.
    const bool closed = (cyc.size() >= 2 && cyc.front() == cyc.back());
    const int last_vertex_index = closed ? (int)cyc.size() - 2 : (int)cyc.size() - 1;
    if (last_vertex_index < 0) { routes_out.clear(); return true; }

    // ---- available starts = starts \ visited_starts ----
    std::unordered_set<int> used(visited_starts.begin(), visited_starts.end());
    std::vector<int> avail_starts;
    avail_starts.reserve(starts.size());
    for (int s : starts) if (!used.count(s)) avail_starts.push_back(s);

    if (remaining_num_path < 0) return false;
    if ((int)avail_starts.size() != remaining_num_path) return false;

    auto take_closest_start = [&](int v_first, int& idx, double& c_sf) -> bool {
        idx = -1; c_sf = std::numeric_limits<double>::infinity();
        for (int i = 0; i < (int)avail_starts.size(); ++i) {
            int s = avail_starts[i];
            double c = costs[s][v_first];
            if (c < c_sf) { c_sf = c; idx = i; }
        }
        return idx >= 0 && std::isfinite(c_sf);
    };

    routes_out.clear();
    routes_out.reserve(remaining_num_path);

    int start_idx = 0; // index of first vertex of the next (uncovered) segment in 'cyc'
    while (start_idx <= last_vertex_index) {
        if (remaining_num_path == 0 || avail_starts.empty()) {
            routes_out.clear();
            return false; // no starts left, yet still uncovered vertices
        }

        const int v_first = cyc[start_idx];

        // Choose closest available start for this segment
        int best_start_pos = -1;
        double c_start_to_first = 0.0;
        if (!take_closest_start(v_first, best_start_pos, c_start_to_first)) {
            routes_out.clear();
            return false;
        }

        // If even single-vertex segment can't fit (start->first + first->t), fail
        const double c_first_to_t = costs[v_first][t_node];
        if (c_start_to_first + c_first_to_t > B) {
            routes_out.clear();
            return false;
        }

        // Greedily extend segment while keeping:
        // c_start_to_first + sum(segment edges) + cost(current_end -> t) <= B
        int end_idx = start_idx;         // include cyc[end_idx]
        double acc_edges = 0.0;

        while (end_idx < last_vertex_index) {
            double next_edge = costs[cyc[end_idx]][cyc[end_idx + 1]];
            double new_c_end_to_t = costs[cyc[end_idx + 1]][t_node];
            if (c_start_to_first + (acc_edges + next_edge) + new_c_end_to_t <= B) {
                acc_edges += next_edge;
                end_idx += 1;
            } else {
                break;
            }
        }

        // Build route: [start, cyc[start_idx..end_idx], t]
        std::vector<int> route;
        route.reserve(2 + (end_idx - start_idx + 1));
        const int chosen_start = avail_starts[best_start_pos];
        route.push_back(chosen_start);
        for (int k = start_idx; k <= end_idx; ++k) route.push_back(cyc[k]);
        route.push_back(t_node);
        routes_out.push_back(std::move(route));

        // consume this start
        visited_starts.push_back(chosen_start);
        avail_starts.erase(avail_starts.begin() + best_start_pos);
        --remaining_num_path;

        // next uncovered vertex (AFTER end_idx)
        start_idx = end_idx + 1;
    }

    return true; // covered entire (open or closed) cycle without wrapping
}



void buildComponentGraph(const std::vector<std::vector<double>>& costs, 
                    const std::vector<int>& component,
                    Graph& g, Graph::EdgeMap<double>& weight) {

    int N = (int)component.size();
    std::vector<Node> nodes(N);
    for (int i = 0; i < N; ++i)
        nodes[i] = g.addNode();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Edge e = g.addEdge(nodes[i], nodes[j]);
            weight[e] = costs[component[i]][component[j]];
        }
    }
}


// Traditional Christofides (tour) without a fixed start
std::vector<int> christofidesTour(const Graph& g, const std::vector<std::vector<double>>& costs,
                                  const std::vector<int>& component,
                                  std::vector<Edge>& mst_edges) {
    // Degrees in the MST
    int N_nodes = component.size();
    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    // Odd-degree vertices of the MST
    std::vector<int> odd_vertices;
    for (int i = 0; i < N_nodes; ++i) if (degree[i] % 2 != 0) odd_vertices.push_back(i);

    // Matching graph on odd vertices (min-weight via negative weights)
    Graph odd_g;
    std::vector<Node> odd_nodes;
    int m = (int)odd_vertices.size();
    odd_nodes.reserve(m);
    for (int i = 0; i < m; ++i) odd_nodes.push_back(odd_g.addNode());

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -costs[component[odd_vertices[i]]][component[odd_vertices[j]]];
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Multigraph = MST + matching edges
    std::vector<std::vector<int>> multigraph(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e));
        int v = g.id(g.v(e));
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    }
    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
        if (matching.matching(e)) {
            int u_idx = odd_g.id(odd_g.u(e));
            int v_idx = odd_g.id(odd_g.v(e));
            int u = odd_vertices[u_idx];
            int v = odd_vertices[v_idx];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    // Eulerian circuit start: first node that has edges
    int start = 0;
    while (start < N_nodes && multigraph[start].empty()) ++start;

    // Hierholzer
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;
    std::vector<int> stack;
    if (start < N_nodes) stack.push_back(start);

    while (!stack.empty()) {
        int v = stack.back();
        if (!temp_graph[v].empty()) {
            int u = temp_graph[v].back();
            temp_graph[v].pop_back();
            auto it = std::find(temp_graph[u].begin(), temp_graph[u].end(), v);
            if (it != temp_graph[u].end()) temp_graph[u].erase(it);
            stack.push_back(u);
        } else {
            circuit.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(circuit.begin(), circuit.end());

    // Shortcut to Hamiltonian cycle order (tour)
    std::vector<bool> visited(N_nodes, false);
    std::vector<int> tour;
    tour.reserve(N_nodes);
    for (int v : circuit) {
        if (v >= 0 && v < N_nodes && !visited[v]) {
            tour.push_back(v);
            visited[v] = true;
        }
    }

    return tour;
}

// void binarySearchBestPath(const std::vector<std::vector<double>>& costs, 
//                     const std::vector<double>& prizes,
//                     const std::vector<int>& indices,
//                     const std::vector<int>& starts,
//                     int t,
//                     int N_min, int N_max, double budget,
//                     std::vector<std::vector<int>>& best_mst_path) {
    
//     double best_prize = 0.0;
//     N_max += 1;
//     while (N_max - N_min > 1) {
//         int N = (N_max + N_min) / 2;
//         std::vector<std::vector<int>> paths;

//         std::vector<int> selected_nodes(indices.begin(), indices.begin()+N);
//         std::vector<std::vector<int>> components = std::move(lambda_components(costs, selected_nodes, budget));
//         std::vector<int> visited_starts;
//         int remaining_num_path = starts.size();
//         bool feasibale = true;
//         paths.clear();
//         for(auto component : components) {
//             Graph g;
//             Graph::EdgeMap<double> weight(g);
//             buildComponentGraph(costs, component, g, weight);

//             std::vector<Edge> mst_edges;
//             buildMST(g, weight, mst_edges);
//             double total_mst_weight = computeMSTWeight(g, weight, mst_edges);

//             std::vector<int> mst_path = std::move(christofidesTour(g, costs, component, N, mst_edges));
//             std::vector<int> mst_path_convert;
//             for (int node : mst_path) {
//                 mst_path_convert.push_back(indices[node]);
//             }        

//             bool has_improve = two_opt_st(mst_path_convert, costs);

//             while(has_improve) {
//                 has_improve = false;
//                 has_improve = two_opt_st(mst_path_convert, costs);
//             }

//             std::vector<std::vector<int>> subpaths;
//             feasiable &= split_cycle(costs, mst_path_convert, starts, visited_starts, remaining_num_path, budget, t, subpaths);
//             if(feasiable) {
//                 paths.insert(paths.end(), subpaths.begin(), subpaths.end());
//             }
//             else break;

//         }

//         if (!feasiable) {
//             N_max = N;
//         } else {
//             repairMultiPaths(costs, prizes, paths, budget, starts,t);

//             double new_prize = multiPathsSumPrize(prizes, paths);
//             if (new_prize > best_prize) {
//                 best_prize = new_prize;
//                 best_mst_paths = paths;
//                 // std::cout << "best_prize = " << best_prize << std::endl;
//             }

//             N_min = N;
//         }
//     }
// }


double min_cost(const std::vector<std::vector<double>>& costs) {
    double m = std::numeric_limits<double>::infinity();
    for (const auto& row : costs) {
        if (row.empty()) continue;
        double rmin = *std::min_element(row.begin(), row.end());
        if (rmin < m) m = rmin;
    }
    return m;
}


double binarySearchBestPathCover(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    const std::vector<int>& indices,
                    const std::vector<int>& starts,
                    int t,
                    double budget,
                    std::vector<std::vector<int>>& best_mst_paths) {
    double best_prize = 0.0;
    double minCost = min_cost(costs);
    
    int N_min = 0, N_max = indices.size();
    // std::ceil(budget * (double)starts.size()/minCost + (double)starts.size());
    // N_max += 1;
    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        // std::cout << "N : " << N <<std::endl;
        std::vector<std::vector<int>> paths;

        std::vector<int> selected_nodes(indices.begin(), indices.begin()+N);
        auto components = lambda_components(costs, selected_nodes, budget);
        // std::cout << "Num Components : " << components.size() <<std::endl;
        std::vector<int> visited_starts;
        int remaining_num_path = starts.size();

        bool feasible = true;

        paths.clear();
        for (const auto& component : components) {
            if (!feasible) break;

            Graph g;
            Graph::EdgeMap<double> weight(g);
            buildComponentGraph(costs, component, g, weight);

            std::vector<Edge> mst_edges;
            buildMST(g, weight, mst_edges);
            double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
            // std::cout << "mst size : " << mst_edges.size() << std::endl;

            std::vector<int> mst_path = christofidesTour(g, costs, component, mst_edges);
            // std::cout << "cycle size : " << mst_path.size() << std::endl;

            std::vector<int> mst_path_convert;
            for (int node : mst_path) {
                mst_path_convert.push_back(component[node]);
            }

            while (two_opt_st(mst_path_convert, costs)) {}

            std::vector<std::vector<int>> subpaths;
            feasible &= split_cycle(costs, mst_path_convert, starts, visited_starts, remaining_num_path, budget, t, subpaths);
            if (feasible) {
                paths.insert(paths.end(), subpaths.begin(), subpaths.end());
            }
            else break;
        }

        if (!feasible) {
            N_max = N;
        } else {
            // Pad unused starts as [start, t] routes 
            if (remaining_num_path > 0) {
                std::unordered_set<int> used(visited_starts.begin(), visited_starts.end());
                for (int s : starts) {
                    if (used.count(s)) continue; 
                    double c_st = costs[s][t];
                    if (!(c_st <= budget)) {  
                        feasible = false;
                        break;
                    }
                    paths.push_back({s, t});
                    visited_starts.push_back(s);
                    --remaining_num_path;
                }
            }
            if (!feasible) {
                N_max = N;
                continue;
            }

            // repairMultiPaths(costs, prizes, paths, budget, starts, t);

            double new_prize = multiPathsSumPrize(prizes, paths);
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_paths = paths;
            }

            N_min = N;
        }
    }
    return best_prize;
}



std::pair<std::vector<std::vector<int>>, double> bestPrizePathCover(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m){

    std::vector<int> indices;
    std::vector<int> indices_all(prizes.size());
    std::iota (indices_all.begin(), indices_all.end(), 0);
    std::sort(indices_all.begin(), indices_all.end(),
            [&](const int& a, const int& b){ return prizes[a] > prizes[b]; });



    for (int i : indices_all) {
        if (i == t) continue; 
        if (std::count(start_indices.begin(), start_indices.end(), i) > 0) continue;
        indices.push_back(i);
    }

    std::vector<std::vector<int>> best_paths;
    double best_prize = binarySearchBestPathCover(costs, prizes, indices, start_indices, t, budget, best_paths);
    for(auto& path : best_paths) {
        while(two_opt_st(path, costs)) {}
    }
    repairMultiPaths(costs, prizes, best_paths, budget, start_indices,t);
    best_prize = multiPathsSumPrize( prizes, best_paths);
    return {best_paths, best_prize};

}


inline void selectNodesPrizeRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices,
                            const std::vector<int>& start_indices,
                            int t) {
    
    const int n = static_cast<int>(costs.size());
    if (n != static_cast<int>(prizes.size())) {
        throw std::invalid_argument("costs/prizes size mismatch");
    }

    std::vector<bool> is_start(n, false);
    for (int s_idx : start_indices) {
        if (s_idx < 0 || s_idx >= n) continue; 
        is_start[s_idx] = true;
    }

    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double,int>> heap;
    const double eps = 1e-12;

    // mark starts selected 
    for (int s_idx : start_indices) {
        selected[s_idx] = true;
    }
    selected[t] = true;

    // init minDist from starts; push only non-start, non-t nodes
    for (int v = 0; v < n; ++v) {
        if (is_start[v] || v == t) continue;
        for (int s_idx : start_indices)
            minDist[v] = std::min(minDist[v], costs[s_idx][v]);
        double score = prizes[v] / (minDist[v] + eps);
        heap.push({score, v});
    }

    // grow set by best prize/distance; push chosen nodes into indices
    while (!heap.empty()) {
        auto [best_score, u] = heap.top(); heap.pop();
        if (selected[u]) continue;  
        selected[u] = true;
        indices.push_back(u);

        for (int v = 0; v < n; ++v) if (!selected[v]) {
            double nd = costs[u][v];
            if (nd < minDist[v]) {
                minDist[v] = nd;
                double new_score = prizes[v] / (minDist[v] + eps);
                heap.push({new_score, v});
            }
        }
    }
}


std::pair<std::vector<std::vector<int>>, double> bestPrizeRatioPathCover(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m){

    std::vector<int> indices;
    selectNodesPrizeRatio(costs, prizes, indices, start_indices,t);
    std::vector<std::vector<int>> best_paths;
    double best_prize = binarySearchBestPathCover(costs, prizes, indices, start_indices, t, budget, best_paths);
    for(auto& path : best_paths) {
        while(two_opt_st(path, costs)) {}
    }
    repairMultiPaths(costs, prizes, best_paths, budget, start_indices,t);
    best_prize = multiPathsSumPrize( prizes, best_paths);
    return {best_paths, best_prize};

}

std::pair<std::vector<std::vector<int>>, double> PathCoverWraper(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m) {

    auto f1 = std::async(std::launch::async, bestPrizePathCover,  costs, prizes, budget, start_indices, t, m);
    auto f2 = std::async(std::launch::async, bestPrizeRatioPathCover, costs, prizes, budget, start_indices, t, m);
    

    std::pair<std::vector<std::vector<int>>, double> p1 = f1.get();
    // return p1;
    std::pair<std::vector<std::vector<int>>, double> p2 = f2.get();

    if (p1.second < p2.second) {
        return p2;
    } else {
        return p1;
    }
}








