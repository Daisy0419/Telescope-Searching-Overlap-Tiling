#include "GCP.h"
#include "PathOperations.h"
#include "LocalSearch.h"
#include "greedy_mink.h"

#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <algorithm>
#include <iostream>


// Components of G(J, λ): keep edges (u,v) with w(u,v) <= λ
// All nodes are customers (no starts to exclude).
std::vector<std::vector<int>> connected_components_unrooted(
        const std::vector<std::vector<double>>& slew_costs,
        const std::vector<double>& dwell_costs,
        double lambda) {

    int N = (int)slew_costs.size();
    if (N == 0) return {};

    std::vector<int> vis(N, 0);
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


// Greedy split of a tour into unrooted open paths.
// Route cost = sum of slew costs on edges + sum of dwell costs on visited nodes.
// Each route is just [customer, customer, ...] with no designated start or end.
bool greedy_split_unrooted(const std::vector<std::vector<double>>& slew_costs,
                           const std::vector<double>& dwell_costs,
                           const std::vector<int>& cyc,
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

    routes_out.clear();

    int start_idx = 0;
    while (start_idx <= last_vertex_index) {
        if (remaining_num_path <= 0) {
            routes_out.clear();
            return false;
        }

        const int v_first = cyc[start_idx];

        // Minimal single-node route: just dwell cost of the first node
        double route_cost = dwell_costs[v_first];
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

        // Build route: [cyc[start_idx..end_idx]]
        std::vector<int> route;
        route.reserve(end_idx - start_idx + 1);
        for (int k = start_idx; k <= end_idx; ++k)
            route.push_back(cyc[k]);
        routes_out.push_back(std::move(route));

        --remaining_num_path;
        start_idx = end_idx + 1;
    }

    return true;
}


void buildComponentGraphUnrooted(const std::vector<std::vector<double>>& slew_costs,
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


std::vector<int> buildChristofidesTourUnrooted(const Graph& g,
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


int UnrootedMinKCover(const std::vector<std::vector<double>>& slew_costs,
                      const std::vector<double>& dwell_costs,
                      double budget,
                      std::vector<std::vector<int>>& best_paths) {


    int N = (int)dwell_costs.size();

    // Global feasibility check: every tile must be coverable alone.
    for (int i = 0; i < N; ++i) {
        if (dwell_costs[i] > budget + 1e-9) {
            std::cerr << "No feasible solution: tile " << i
                      << " has dwell time " << dwell_costs[i]
                      << " > budget " << budget << "\n";
            best_paths.clear();
            return 0;
        }
    }

    auto components = connected_components_unrooted(slew_costs, dwell_costs, budget);

    int max_routes = dwell_costs.size();
    int remaining_num_path = max_routes;
    bool feasible = true;

    std::vector<std::vector<int>> paths;

    for (const auto& component : components) {
        if (!feasible) break;

        // Single-node component: just check dwell fits in budget
        if (component.size() == 1) {
            if (remaining_num_path <= 0) { feasible = false; break; }
            if (dwell_costs[component[0]] > budget) { feasible = false; break; }
            paths.push_back({component[0]});
            --remaining_num_path;
            continue;
        }

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildComponentGraphUnrooted(slew_costs, component, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        std::vector<int> tour = buildChristofidesTourUnrooted(g, slew_costs, component, mst_edges);

        // Map local indices to global node IDs
        std::vector<int> tour_global;
        tour_global.reserve(tour.size());
        for (int node : tour)
            tour_global.push_back(component[node]);

        while (two_opt_st(tour_global, slew_costs)) {}

        std::vector<std::vector<int>> subpaths;
        feasible = greedy_split_unrooted(slew_costs, dwell_costs, tour_global,
                                         remaining_num_path, budget, subpaths);
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