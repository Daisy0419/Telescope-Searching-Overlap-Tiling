#include "GCP.h"
#include "PathOperations.h"
#include "LocalSearch.h"
#include "greedy_mink.h"

#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>
#include <tuple>
#include <iostream>
#include <unordered_set>


// Components of G(J, λ): keep edges (u,v) with w(u,v) <= λ
// Only considers customer nodes (excludes start nodes).
std::vector<std::vector<int>> connected_components(
        const std::vector<std::vector<double>>& slew_costs,
        const std::vector<double>& dwell_costs,
        const std::unordered_set<int>& start_set,
        double lambda) {

    int N = (int)slew_costs.size();
    if (N == 0) return {};

    std::vector<int> vis(N, 0);
    // Mark starts as already visited so they never enter any component
    for (int s : start_set) vis[s] = 1;

    std::vector<std::vector<int>> comps;
    for (int i = 0; i < N; ++i) {
        if (vis[i]) continue;
        std::vector<int> comp;
        std::queue<int> q;
        q.push(i);
        vis[i] = 1;
        while (!q.empty()) {
            int u = q.front(); q.pop();
            comp.push_back(u);
            for (int v = 0; v < N; ++v) {
                if (vis[v] || v == u) continue;
                // Edge weight for connectivity test: slew + destination dwell
                if (slew_costs[u][v] + dwell_costs[v] <= lambda) {
                    vis[v] = 1;
                    q.push(v);
                }
            }
        }
        comps.push_back(std::move(comp));
    }
    return comps;
}


// Greedy split of a tour into open-ended rooted paths.
// Route cost = slew costs on edges + dwell costs on visited customer nodes.
// No terminal node: each route is [start, customer, customer, ...].
bool greedy_split(const std::vector<std::vector<double>>& slew_costs,
                  const std::vector<double>& dwell_costs,
                  const std::vector<int>& cyc,
                  const std::vector<int>& starts,
                  std::vector<int>& visited_starts,
                  int& remaining_num_path,
                  double B,
                  std::vector<std::vector<int>>& routes_out) {

    if (B <= 0.0) {
        std::cerr << "Budget B must be positive.\n";
        return false;
    }
    const int N = (int)slew_costs.size();
    if (cyc.empty()) {
        routes_out.clear();
        return true;
    }
    if (N == 0 || (int)slew_costs[0].size() != N) return false;

    const bool closed = (cyc.size() >= 2 && cyc.front() == cyc.back());
    const int last_vertex_index = closed ? (int)cyc.size() - 2 : (int)cyc.size() - 1;
    if (last_vertex_index < 0) { routes_out.clear(); return true; }

    // available starts = starts \ visited_starts
    std::unordered_set<int> used(visited_starts.begin(), visited_starts.end());
    std::vector<int> avail_starts;
    avail_starts.reserve(starts.size());
    for (int s : starts)
        if (!used.count(s)) avail_starts.push_back(s);

    if (avail_starts.empty()) {
        std::cerr << "Not enough starts to cover the graph!\n";
        return false;
    }

    auto take_closest_start = [&](int v_first, int& idx, double& c_sf) -> bool {
        idx = -1;
        c_sf = std::numeric_limits<double>::infinity();
        for (int i = 0; i < (int)avail_starts.size(); ++i) {
            double c = slew_costs[avail_starts[i]][v_first];
            if (c < c_sf) { c_sf = c; idx = i; }
        }
        return idx >= 0 && std::isfinite(c_sf);
    };

    routes_out.clear();

    int start_idx = 0;
    while (start_idx <= last_vertex_index) {
        if (avail_starts.empty() || remaining_num_path <= 0) {
            routes_out.clear();
            return false;
        }

        const int v_first = cyc[start_idx];

        int best_start_pos = -1;
        double c_start_to_first = 0.0;
        if (!take_closest_start(v_first, best_start_pos, c_start_to_first)) {
            routes_out.clear();
            return false;
        }

        // Cost of the minimal single-node route: slew(start->first) + dwell(first)
        double route_cost = c_start_to_first + dwell_costs[v_first];
        if (route_cost > B) {
            routes_out.clear();
            return false;
        }

        // Greedily extend while budget allows
        int end_idx = start_idx;

        while (end_idx < last_vertex_index) {
            int cur_node = cyc[end_idx];
            int next_node = cyc[end_idx + 1];
            double extend_cost = slew_costs[cur_node][next_node] + dwell_costs[next_node];
            if (route_cost + extend_cost <= B) {
                route_cost += extend_cost;
                ++end_idx;
            } else {
                break;
            }
        }

        // Build route: [start, cyc[start_idx..end_idx]]  (no terminal)
        const int chosen_start = avail_starts[best_start_pos];
        std::vector<int> route;
        route.reserve(1 + (end_idx - start_idx + 1));
        route.push_back(chosen_start);
        for (int k = start_idx; k <= end_idx; ++k)
            route.push_back(cyc[k]);
        routes_out.push_back(std::move(route));

        visited_starts.push_back(chosen_start);
        avail_starts.erase(avail_starts.begin() + best_start_pos);
        --remaining_num_path;

        start_idx = end_idx + 1;
    }

    return true;
}


void buildComponentGraphRooted(const std::vector<std::vector<double>>& slew_costs,
                         const std::vector<int>& component,
                         Graph& g, Graph::EdgeMap<double>& weight) {
    int N = (int)component.size();
    std::vector<Node> nodes(N);
    for (int i = 0; i < N; ++i)
        nodes[i] = g.addNode();
    for (int i = 0; i < N; ++i) {
        for (int j = i + 1; j < N; ++j) {
            Edge e = g.addEdge(nodes[i], nodes[j]);
            weight[e] = slew_costs[component[i]][component[j]];
        }
    }
}


std::vector<int> buildChristofidesTour(const Graph& g,
                                       const std::vector<std::vector<double>>& slew_costs,
                                       const std::vector<int>& component,
                                       std::vector<Edge>& mst_edges) {
    int N_nodes = (int)component.size();

    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    std::vector<int> odd_vertices;
    for (int i = 0; i < N_nodes; ++i)
        if (degree[i] % 2 != 0) odd_vertices.push_back(i);

    Graph odd_g;
    int m = (int)odd_vertices.size();
    std::vector<Node> odd_nodes;
    odd_nodes.reserve(m);
    for (int i = 0; i < m; ++i)
        odd_nodes.push_back(odd_g.addNode());

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -slew_costs[component[odd_vertices[i]]][component[odd_vertices[j]]];
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    std::vector<std::vector<int>> multigraph(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e));
        int v = g.id(g.v(e));
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    }
    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e) {
        if (matching.matching(e)) {
            int u = odd_vertices[odd_g.id(odd_g.u(e))];
            int v = odd_vertices[odd_g.id(odd_g.v(e))];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    int start = 0;
    while (start < N_nodes && multigraph[start].empty()) ++start;

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


int RootedOpenMinKCover(const std::vector<std::vector<double>>& slew_costs,
                        const std::vector<double>& dwell_costs,
                        const std::vector<int>& starts,
                        double budget,
                        std::vector<std::vector<int>>& best_paths) {

    // Build start set so components exclude start nodes
    std::unordered_set<int> start_set(starts.begin(), starts.end());

    auto components = connected_components(slew_costs, dwell_costs, start_set, budget);

    std::vector<int> visited_starts;
    int remaining_num_path = (int)starts.size();
    bool feasible = true;

    std::vector<std::vector<int>> paths;

    for (const auto& component : components) {
        if (!feasible) break;


        if (component.size() == 1) {
            int v = component[0];

            std::unordered_set<int> used(visited_starts.begin(), visited_starts.end());

            int best_s = -1;
            double best_cost = std::numeric_limits<double>::infinity();

            for (int s : starts) {
                if (used.count(s)) continue;

                double c = slew_costs[s][v] + dwell_costs[v];
                if (c < best_cost) {
                    best_cost = c;
                    best_s = s;
                }
            }

            if (best_s < 0 || best_cost > budget) {
                feasible = false;
                break;
            }

            paths.push_back({best_s, v});
            visited_starts.push_back(best_s);
            --remaining_num_path;
            continue;
        }

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildComponentGraphRooted(slew_costs, component, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        std::vector<int> tour = buildChristofidesTour(g, slew_costs, component, mst_edges);

        // Map local indices to global node IDs
        std::vector<int> tour_global;
        tour_global.reserve(tour.size());
        for (int node : tour)
            tour_global.push_back(component[node]);

        while (two_opt_st(tour_global, slew_costs)) {}

        std::vector<std::vector<int>> subpaths;
        feasible = greedy_split(slew_costs, dwell_costs, tour_global, starts,
                                visited_starts, remaining_num_path, budget, subpaths);
        if (feasible) {
            paths.insert(paths.end(), subpaths.begin(), subpaths.end());
        }
    }

    if (!feasible) {
        std::cerr << "No feasible solution!\n";
        return 0;
    }

    best_paths = std::move(paths);
    return (int)best_paths.size();
}


