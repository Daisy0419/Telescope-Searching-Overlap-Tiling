#include "GCP.h"
#include "ReadData.h"
#include "helplers.h"
#include "LocalSearch.h"
#include "PathOperations.h"

#include <vector>
#include <list>
#include <queue>
#include <limits>
#include <algorithm>
#include <utility>
#include <tuple>
#include <functional>
#include <unordered_set>
#include <numeric>
#include <iostream>
#include <stdexcept>
#include <stack>
#include <fstream>
#include <sstream>
#include <future>
#include <cmath>
#include <chrono>


//sort by prize
void sortNodes(const std::vector<double>& prizes, std::vector<double>& prizes_sorted, 
               std::vector<int>& indices, int init_pos) {
    int n = prizes.size();
    if (init_pos < 0 || init_pos >= n) {
        prizes_sorted.clear();
        indices.clear();
        return;
    }
    std::vector<int> indices_temp(n);
    std::iota(indices_temp.begin(), indices_temp.end(), 0); // 0-based

    std::sort(indices_temp.begin(), indices_temp.end(), 
              [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

    prizes_sorted.resize(n);
    indices.resize(n);
    prizes_sorted[0] = prizes[init_pos];
    indices[0] = init_pos;
    size_t j = 1;
    for (size_t i = 0; i < n; ++i) {
        if (indices_temp[i] == init_pos)
            continue;
        prizes_sorted[j] = prizes[indices_temp[i]];
        indices[j] = indices_temp[i];
        ++j;
    }
}


// sort nodes in a given range by prize
void sortNodesInRange(const std::vector<double>& prizes, std::vector<int>& indices, std::pair<int, int> range) {
    int n = prizes.size();
    if (range.second <= 0 || range.first < 0 || range.first >= range.second) {
        std::cerr << "Invalid inputs to sortNodesInRange\n";
        return;
    }
    if (range.second > n) range.second = n;

    //partial partition
    std::nth_element(indices.begin() + range.first, indices.begin() + range.second, indices.end(),
                     [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

    //sort within [range.first, range.second)
    std::sort(indices.begin() + range.first, indices.begin() + range.second,
              [&prizes](int a, int b) { return prizes[a] > prizes[b]; });

}

void buildGraphUpdate(const std::vector<std::vector<double>>& costs, 
            const std::vector<int>& indices, int N_nodes,
            Graph& g, Graph::EdgeMap<double>& weight) {
    std::vector<Node> nodes(N_nodes);
    for (int i = 0; i < N_nodes; ++i)
        nodes[i] = g.addNode();
    for (int i = 0; i < N_nodes; ++i) {
        for (int j = i + 1; j < N_nodes; ++j) {
            Edge e = g.addEdge(nodes[i], nodes[j]);
            weight[e] = costs[indices[i]][indices[j]];
        }
    }
}

void buildMST(const Graph& g, const Graph::EdgeMap<double>& weight, std::vector<Edge>& mst_edges) {
    mst_edges.clear();
    lemon::kruskal(g, weight, std::back_inserter(mst_edges));
}

double computeMSTWeight(const Graph& g, const Graph::EdgeMap<double>& weight, 
                        const std::vector<Edge>& mst_edges) {
    double total_weight = 0.0;
    for (const Edge& e : mst_edges) {
        total_weight += weight[e];
    }
    return total_weight;
}


//  Find odd-degree vertices (J.A. Hoogeveen)
std::vector<int> findOddDegreeVertices(const std::list<std::pair<int, int>>& mst) {
    std::unordered_map<int, int> degree;
    for (const auto& [u, v] : mst) {
        degree[u]++;
        degree[v]++;
    }
    std::vector<int> odd;
    for (const auto& [node, deg] : degree) {
        if (deg % 2 == 1)
            odd.push_back(node);
    }
    return odd;
}


// Find vertices of "wrong degree" for one fixed endpoint s (J.A. Hoogeveen)
std::vector<int> findWrongDegreeVertices(const std::list<std::pair<int, int>>& mst, int s, int numNodes) {
    std::unordered_map<int, int> degree;
    for (const auto& [u, v] : mst) {
        degree[u]++;
        degree[v]++;
    }
    std::vector<int> S;
    for (int v = 0; v < numNodes; ++v) {
        if (v == s) {
            if (degree[v] % 2 == 0) { // s has even degree
                S.push_back(v);
            }
        } else {
            if (degree[v] % 2 == 1) { // v != s has odd degree
                S.push_back(v);
            }
        }
    }
    return S;
}

// Greedy matching (leaves one vertex exposed if |S| is odd)
std::list<std::pair<int, int>> greedyMatching(const std::vector<int>& nodes,
                                                 const std::vector<std::vector<double>>& costMatrix) {
    std::list<std::pair<int, int>> matching;
    std::unordered_set<int> used;

    for (size_t i = 0; i < nodes.size(); ++i) {
        int u = nodes[i];
        if (used.count(u)) continue;

        double minCost = std::numeric_limits<double>::infinity();
        int bestV = -1;

        for (size_t j = i + 1; j < nodes.size(); ++j) {
            int v = nodes[j];
            if (used.count(v)) continue;

            if (costMatrix[u][v] < minCost) {
                minCost = costMatrix[u][v];
                bestV = v;
            }
        }

        if (bestV != -1) {
            matching.emplace_back(u, bestV);
            used.insert(u);
            used.insert(bestV);
        }
    }
    return matching;
}

// Traditional Christofides (tour) without a fixed start
std::vector<int> christofidesTour(const Graph& g, const std::vector<std::vector<double>>& costs,
                                  const std::vector<int>& indices, int N_nodes,
                                  std::vector<Edge>& mst_edges) {
    // Degrees in the MST
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
            odd_weight[e] = -costs[indices[odd_vertices[i]]][indices[odd_vertices[j]]];
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
    // if (!tour.empty()) tour.push_back(tour.front());

    return tour;
}


// Traditional Christofides (tour)
std::vector<int> christofidesPath(const Graph& g, const std::vector<std::vector<double>>& costs,
                                          const std::vector<int>& indices, int N_nodes,
                                          std::vector<Edge>& mst_edges, int fixed_start_idx) {
    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    // Use standard Christofides: vertices of odd degree in the MST
    std::vector<int> odd_vertices;
    for (int i = 0; i < N_nodes; ++i) {
        if (degree[i] % 2 != 0) odd_vertices.push_back(i);
    }

    // Construct matching graph on odd vertices (no dummy node)
    Graph odd_g;
    std::vector<Node> odd_nodes;
    int m = (int)odd_vertices.size();
    for (int i = 0; i < m; ++i) odd_nodes.push_back(odd_g.addNode());

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            // MaxWeightedPerfectMatching on negative weights == min-weight matching
            odd_weight[e] = -costs[indices[odd_vertices[i]]][indices[odd_vertices[j]]];
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Build multigraph (MST + matching edges)
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

    // Eulerian circuit (all degrees are now even)
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;

    // Start at fixed_start_idx if it has edges; else first node with edges
    int start = fixed_start_idx;
    if (start < 0 || start >= N_nodes || temp_graph[start].empty()) {
        start = 0;
        for (int i = 0; i < N_nodes; ++i) if (!temp_graph[i].empty()) { start = i; break; }
    }
    std::vector<int> stack = { start };

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

    // Shortcut to Hamiltonian cycle order (tour); you can append start again if you want explicit closure
    std::vector<bool> visited(N_nodes, false);
    std::vector<int> path;
    for (int v : circuit) {
        if (!visited[v]) {
            path.push_back(v);
            visited[v] = true;
        }
    }
    // If you prefer an explicit cycle, uncomment:
    // if (!path.empty()) path.push_back(path.front());

    return path;
}


// Adapted christofides algorithm for a path with one fixed end (J.A. Hoogeveen)
std::vector<int> christofidesPathOneFixed(const Graph& g, const std::vector<std::vector<double>>& costs,
                                             const std::vector<int>& indices, int N_nodes,
                                             std::vector<Edge>& mst_edges, int fixed_start_idx) {
    std::vector<int> degree(N_nodes, 0);
    for (Edge e : mst_edges) {
        degree[g.id(g.u(e))]++;
        degree[g.id(g.v(e))]++;
    }

    // Determine vertices of wrong degree
    std::vector<int> wrong_vertices;
    for (int i = 0; i < N_nodes; ++i) {
        bool is_fixed = (i == fixed_start_idx);
        if ((is_fixed && degree[i] % 2 == 0) || (!is_fixed && degree[i] % 2 != 0)) {
            wrong_vertices.push_back(i);
        }
    }

    // Construct matching graph with dummy node if needed
    Graph odd_g;
    std::vector<Node> odd_nodes;
    int m = wrong_vertices.size();
    bool needs_dummy = (m % 2 != 0);

    for (int i = 0; i < m; ++i)
        odd_nodes.push_back(odd_g.addNode());

    Node dummy_node;
    if (needs_dummy) {
        dummy_node = odd_g.addNode();
        odd_nodes.push_back(dummy_node);
    }

    Graph::EdgeMap<double> odd_weight(odd_g);
    for (int i = 0; i < m; ++i) {
        for (int j = i + 1; j < m; ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_weight[e] = -costs[indices[wrong_vertices[i]]][indices[wrong_vertices[j]]];
        }
        if (needs_dummy) {
            Edge e = odd_g.addEdge(odd_nodes[i], dummy_node);
            odd_weight[e] = -1e6;  // large negative weight
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_weight);
    matching.run();

    // Build multigraph (MST + matching edges)
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
            // std::cout << "u_idx: " << u_idx << ", v_idx: " << v_idx << "\n";
            if (needs_dummy && (u_idx == m || v_idx == m)) continue; // ignore dummy edges
            int u = wrong_vertices[u_idx];
            int v = wrong_vertices[v_idx];
            multigraph[u].push_back(v);
            multigraph[v].push_back(u);
        }
    }

    // std::cout << "multigraph side: " << multigraph.size() << std::endl;
    // Find Eulerian path starting from fixed_start_idx
    std::vector<int> circuit;
    std::vector<std::vector<int>> temp_graph = multigraph;
    std::vector<int> stack = {fixed_start_idx};

    while (!stack.empty()) {
        int v = stack.back();
        if (!temp_graph[v].empty()) {
            int u = temp_graph[v].back();
            // if (u >= N_nodes || v >= N_nodes) {
            //     std::cerr << "ERROR: out-of-bounds node: u=" << u << ", v=" << v << std::endl;
            // }
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
    // std::cout << "circuit.size()" << circuit.size() << std::endl;

    // Shortcut
    std::vector<bool> visited(N_nodes, false);
    std::vector<int> path;
    for (int v : circuit) {
        if (!visited[v]) {
            path.push_back(v);
            visited[v] = true;
        }
    }
    // std::cout << "path.size()" << path.size() << std::endl;

    return path;
}

void binarySearchBestPathOneFixed(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    const std::vector<int>& indices,
                    std::vector<int>& best_mst_path,
                    int N_min, int N_max, double budget) {
    
    double best_prize = 0.0;
    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        std::vector<int> mst_path = std::move(christofidesPathOneFixed(g, costs, indices, N, mst_edges, 0));
        std::vector<int> mst_path_convert;
        for (int node : mst_path) {
            mst_path_convert.push_back(indices[node]);
        }   
        // auto start = std::chrono::high_resolution_clock::now();
        while(two_opt_st(mst_path_convert, costs)){};
        // auto end = std::chrono::high_resolution_clock::now();
        // auto elapsed_seconds = end - start;
        // std::cout << "one pass 2-opt running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
        // three_opt_st(mst_path_convert, costs);

        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path_convert.size(); ++i) {
            pathCost += costs[mst_path_convert[i-1]][mst_path_convert[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            double new_prize = 0.0;
            for (int node : mst_path_convert) {
                new_prize += prizes[node];
            }
            //update best path
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path_convert;
            }

            N_min = N;
        }
    }
}


//s-any path - Best Prize Node selection with one fixed end Christofedis
std::vector<int> GCP_BestPrizeOneFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {
    if (init_pos < 0 || init_pos >= (int)prizes.size() || budget <= 0) {
        return {init_pos};
    }

    std::vector<int> indices(prizes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::swap(indices[0], indices[init_pos]);

    int N = 2;
    while(true) {
        sortNodesInRange(prizes, indices, {N/2, N});

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
            break;
        }
    }
    
    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


// Adapted christofides algorithm for a path with two fixed ends (J.A. Hoogeveen)
std::vector<int> christofidesPathTwoFixed(const Graph& g,
                         const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& indices,
                         int N_nodes,
                         const std::vector<Edge>& mst_edges,
                         int s_idx,
                         int t_idx) {

    //build adjacency lists
    std::vector<std::vector<int>> adj(N_nodes);
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e)), v = g.id(g.v(e));
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    std::vector<int> parent(N_nodes, -1);
    std::queue<int> Q;
    Q.push(s_idx);
    parent[s_idx] = s_idx;
    while (!Q.empty() && parent[t_idx] == -1) {
        int u = Q.front(); Q.pop();
        for (int v : adj[u]) {
            if (parent[v] == -1) {
                parent[v] = u;
                Q.push(v);
            }
        }
    }

    // find heaviest edge on the s-t path in MST
    double max_cost = -1;
    std::pair<int, int> e_to_remove;
    for (int v = t_idx; v != s_idx; v = parent[v]) {
        int u = parent[v];
        double w = costs[indices[u]][indices[v]];
        if (w > max_cost) {
            max_cost = w;
            e_to_remove = {u, v};
        }
    }

    //Build T_bar = (T − {e_max}) ∪ {s,t}
    std::vector<std::pair<int, int>> tree_edges;
    for (Edge e : mst_edges) {
        int u = g.id(g.u(e)), v = g.id(g.v(e));
        if (!((u == e_to_remove.first && v == e_to_remove.second) ||
              (u == e_to_remove.second && v == e_to_remove.first)))
        {
            tree_edges.emplace_back(u, v);
        }
    }
    tree_edges.emplace_back(s_idx, t_idx);  // dummy edge

    //find odd-degree vertices
    std::vector<int> degree(N_nodes, 0);
    for (auto [u, v] : tree_edges) {
        degree[u]++;
        degree[v]++;
    }

    std::vector<int> odd;
    for (int i = 0; i < N_nodes; ++i)
        if (degree[i] % 2 != 0)
            odd.push_back(i);

    //compute minimum-weight perfect matching on odd vertices
    Graph odd_g;
    std::vector<Node> odd_nodes;
    for (int v : odd)
        odd_nodes.push_back(odd_g.addNode());

    Graph::EdgeMap<double> odd_w(odd_g);
    for (size_t i = 0; i < odd.size(); ++i) {
        for (size_t j = i + 1; j < odd.size(); ++j) {
            Edge e = odd_g.addEdge(odd_nodes[i], odd_nodes[j]);
            odd_w[e] = -costs[indices[odd[i]]][indices[odd[j]]];
        }
    }

    lemon::MaxWeightedPerfectMatching<Graph, Graph::EdgeMap<double>> matching(odd_g, odd_w);
    matching.run();

    //create multigraph with T_bar ∪ matching edges
    std::vector<std::vector<int>> multigraph(N_nodes);
    auto add = [&](int u, int v) {
        multigraph[u].push_back(v);
        multigraph[v].push_back(u);
    };
    for (auto [u, v] : tree_edges) add(u, v);
    for (Graph::EdgeIt e(odd_g); e != lemon::INVALID; ++e)
        if (matching.matching(e)) {
            int u = odd[odd_g.id(odd_g.u(e))];
            int v = odd[odd_g.id(odd_g.v(e))];
            add(u, v);
        }

    //compute eulerian cycle
    std::vector<int> cyc, stack = {s_idx};
    auto temp = multigraph;
    while (!stack.empty()) {
        int v = stack.back();
        if (!temp[v].empty()) {
            int u = temp[v].back();
            temp[v].pop_back();
            auto it = std::find(temp[u].begin(), temp[u].end(), v);
            if (it != temp[u].end()) temp[u].erase(it);
            stack.push_back(u);
        } else {
            cyc.push_back(v);
            stack.pop_back();
        }
    }
    std::reverse(cyc.begin(), cyc.end());
    cyc.pop_back(); // remove duplicate start

    //break Eulerian cycle at dummy edge (s,t)
    size_t cut = cyc.size();
    for (size_t i = 0; i < cyc.size(); ++i) {
        size_t j = (i + 1) % cyc.size();
        int u = cyc[i], v = cyc[j];
        if ((u == s_idx && v == t_idx) || (u == t_idx && v == s_idx)) {
            cut = i;
            break;
        }
    }
    if (cut == cyc.size()) throw std::logic_error("Dummy edge not found!");

    std::rotate(cyc.begin(), cyc.begin() + cut, cyc.end());

    //shortcut duplicate nodes (keep first occurrence)
    // for(int p : cyc) {
    //     std::cout << p << " ";
    // }
    // std::cout << "\n";      
    std::vector<bool> seen(N_nodes, false);
    std::vector<int> path_rev;
    for (int i = 0; i < (int)cyc.size(); ++i) {
        int v = cyc[i];
        if (!seen[v]) {
            path_rev.push_back(v);
            seen[v] = true;
        }
    }
    std::reverse(path_rev.begin() + 1, path_rev.end());


    if (path_rev.front() != s_idx || path_rev.back() != t_idx){
        std::reverse(path_rev.begin(), path_rev.end());
        if (path_rev.front() != s_idx || path_rev.back() != t_idx)
            throw std::logic_error("Christofides 2-fixed: endpoints wrong.");
    }

    // for(int p : path_rev) {
    //     std::cout << p << " ";
    // }
    // std::cout << "\n";      
    return path_rev;
}


// Greedy extend the path by adding a node in the end before t
void greedyExtendPath(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes,
                      std::vector<int>& path, double& total_cost, double budget) {
    int end_city = path.back();
    std::vector<bool> visited(costs.size(), false);
    for (int city : path) visited[city] = true;

    while (true) {
        double max_prize = -1.0;
        int best_insert_pos = -1;
        int next_city = -1;

        // Always consider inserting before the final node
        int before_end = path[path.size() - 2];

        for (int i = 0; i < (int)costs.size(); ++i) {
            if (visited[i]) continue;

            double extra_cost = costs[before_end][i] + costs[i][end_city] - costs[before_end][end_city];
            if (total_cost + extra_cost <= budget && prizes[i] > max_prize) {
                max_prize = prizes[i];
                next_city = i;
            }
        }

        if (next_city == -1) break;

        // Insert next_city before end_city
        path.insert(path.end() - 1, next_city);
        visited[next_city] = true;

        // Recalculate last segment before end
        int before_end_new = path[path.size() - 2];
        int before_before = path[path.size() - 3];
        total_cost += costs[before_before][before_end_new] + costs[before_end_new][end_city]
                      - costs[before_before][end_city];
    }
}



void binarySearchBestPathTwoFixed(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    const std::vector<int>& indices,
                    std::vector<int>& best_mst_path,
                    int N_min, int N_max, double budget) {
    
    double best_prize = 0.0;
    N_max += 1;
    while (N_max - N_min > 1) {
        int N = (N_max + N_min) / 2;
        // std::cout << "N_max: " << N_max << " N_min: " << N_min << " N: " << N << std::endl;;
        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        // std::cout << "MST weight = " << total_mst_weight << std::endl;
        if(total_mst_weight > budget) {
            N_max = N;
            continue;
        }
        if (N == 2) {
            double path_cost = costs[indices[0]][indices[1]];
            if (path_cost <= budget) {
                best_mst_path = {indices[0], indices[1]};
            } else {
                best_mst_path = {indices[0]};
            }
            return;
        }
        std::vector<int> mst_path = std::move(christofidesPathTwoFixed(g, costs, indices, N, mst_edges, 0, 1));
        std::vector<int> mst_path_convert;
        for (int node : mst_path) {
            mst_path_convert.push_back(indices[node]);
        }  
        // for(int p : mst_path_convert) {
        //     std::cout << p << " ";
        // }
        // std::cout << "\n";      
        // //2-opt to remove cross
        while(two_opt_st(mst_path_convert, costs)){};
        // // std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        // bool has_improve = two_opt_st(mst_path_convert, costs);
        // // std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        // // std::chrono::duration<double> elapsed_seconds = end - start;
        // // std::cout << "one pass 2-opt running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
        
        // // int count = 1;
        // while(has_improve) {
        //     has_improve = false;
        //     has_improve = two_opt_st(mst_path_convert, costs);
        //     // ++ count;
        // }
        // std::cout << "count : " << count << std::endl;
        // three_opt_st(mst_path_convert, costs);


        double pathCost = 0.0;
        for (size_t i = 1; i < mst_path_convert.size(); ++i) {
            pathCost += costs[mst_path_convert[i-1]][mst_path_convert[i]];
        }

        if (pathCost > budget) {
            N_max = N;
        } else {
            // greedyExtendPath(costs, prizes, mst_path_convert, pathCost, budget);
            // for(int p : mst_path_convert) {
            //     std::cout << p << " ";
            // }
            // std::cout << "\n";
            std::vector<bool> visited(costs.size(), false);
            for(auto p : mst_path_convert) visited[p] = true;
            insertHighValueCheapest(mst_path_convert, costs, prizes, budget, pathCost, visited);
            // insertBestRatio(mst_path_convert, costs, prizes, budget, pathCost, visited);
            // for(int p : mst_path_convert) {
            //     std::cout << p << " ";
            // }
            // std::cout << "\n";

            double new_prize = 0.0;
            for (int node : mst_path_convert) {
                new_prize += prizes[node];
            }
            // std::cout << "new_prize: " << new_prize << std::endl;
            //update best path
            if (new_prize > best_prize) {
                best_prize = new_prize;
                best_mst_path = mst_path_convert;
                // std::cout << "best_prize = " << best_prize << std::endl;
            }

            N_min = N;
        }
    }
}


// //s-t path - Best Prize Node selection with one fixed end Christofedis
// std::vector<int> GCP_BestPrizeTwoFix(const std::vector<std::vector<double>>& costs, 
//                           const std::vector<double>& prizes, 
//                           double budget, int s_idx, int t_idx) {
                     
//     if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
//         return {s_idx};
//     }

//     std::vector<int> indices(prizes.size());
//     std::iota(indices.begin(), indices.end(), 0);
//     std::swap(indices[0], indices[s_idx]);
//     std::swap(indices[1], indices[t_idx]);

    
//     int N_pre = 2;
//     int N = std::min(int(prizes.size()), int(2*N_pre));
//     while(true) {
//         if(N < 3) break;
//         sortNodesInRange(prizes, indices, {N_pre, N});

//         Graph g;
//         Graph::EdgeMap<double> weight(g);
//         buildGraphUpdate(costs, indices, N, g, weight);

//         std::vector<Edge> mst_edges;
//         buildMST(g, weight, mst_edges);

//         double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
//         if(total_mst_weight > budget) {
//             break;
//         }

//         if (N >= (int)prizes.size())
//             break;
        
//         N_pre = N, N *= 2;
//         if (N >= (int)prizes.size()) { 
//             N = prizes.size();
//         }

//     }

//     int N_min = N/3, N_max = N;
//     std::vector<int> best_mst_path;
//     binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

//     while(best_mst_path.empty()) {
//         N_max = N_min, N_min /= 2;
//         binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
//     }
//     return best_mst_path;
// }


std::vector<int> GCP_BestPrizeTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
                     
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    if (costs[s_idx][t_idx] > budget) {
        return {};
    }

    std::vector<int> indices(prizes.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::swap(indices[0], indices[s_idx]);
    std::swap(indices[1], indices[t_idx]);

    int N_pre = 2;
    int N = std::min(int(prizes.size()), int(2*N_pre));
    while(true) {
        if(N < 3) break;
        sortNodesInRange(prizes, indices, {N_pre, N});

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N;
        N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }
    }

    // N_min must be at least 2 (s and t are always included)
    int N_min = std::max(2, N/3);
    int N_max = N;

    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes, indices, best_mst_path, N_min, N_max, budget);

    // Fall back by shrinking the window, but never below the base case of 2 nodes
    while(best_mst_path.empty()) {
        if (N_max <= 2) {
            // Direct s->t is the only option; we already verified it fits in budget
            return {s_idx, t_idx};
        }
        N_max = N_min;
        N_min = std::max(2, N_min / 2);
        binarySearchBestPathTwoFixed(costs, prizes, indices, best_mst_path, N_min, N_max, budget);

        // Prevent infinite loop if N_min == N_max (both are 2)
        if (N_min == N_max) {
            // Still empty means even 2-node path wasn't returned by search; force it
            return {s_idx, t_idx};
        }
    }

    return best_mst_path;
}


// Select nodes using nonlinear prize boosting: score = asin(prize^4) / (minDist + epsilon)
// This emphasizes higher-prize nodes more aggressively and damps extreme values using asin().
//^4 is a tuned value, may not be optimal
void selectNodesByPrizeBoosting(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices, 
                            std::vector<bool>& selected,
                            std::vector<double>& minDist,
                            std::priority_queue<std::pair<double, int>>& heap,
                            int num) {

    int n = prizes.size();
    const double epsilon = 1;

    while (num > 0 && !heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);
        num --;
        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = std::asin(std::pow(prizes[v], 4))/ (minDist[v] + epsilon);
                heap.push({new_score, v});
            }
        }
    }

}

//s-any path - Prize Boosting Node selection with one fixed end Christofedis
std::vector<int> GCP_PrizeBoostOneFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, 
                          int init_pos) {

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[init_pos] = true;
    indices.push_back(init_pos);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != init_pos) {
            minDist[v] = costs[init_pos][v];
            double score = std::asin(std::pow(prizes[v], 4)) / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 1, N = 2;
    while(true) {
        int num_to_select = N - N_pre;
        selectNodesByPrizeBoosting(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if(N == (int)prizes.size())
            break;

        N_pre = N;
        N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }
    }

    int N_min = N/3, N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathOneFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
                
}

//s-t path
std::vector<int> GCP_PrizeBoostTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[s_idx] = true;
    indices.push_back(s_idx);
    selected[t_idx] = true;
    indices.push_back(t_idx);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != s_idx && v != t_idx ) {
            minDist[v] = std::min(costs[s_idx][v], costs[t_idx][v]);
            double score = std::asin(std::pow(prizes[v], 4)) / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 2;
    int N = std::min(int(prizes.size()), int(2*N_pre));
    while(true) {
        if(N < 3) break;
        int num_to_select = N - N_pre;
        selectNodesByPrizeBoosting(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N, N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }

    }
    
    int N_min = std::max(2, N/3);
    int N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}



// Select nodes greedily based on prize-to-distance ratio: score = prize / (minDist + epsilon)
// Which favors nodes that are close to already selected nodes and have high prize values.
void selectNodesByPrizeDistanceRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices, 
                            std::vector<bool>& selected,
                            std::vector<double>& minDist,
                            std::priority_queue<std::pair<double, int>>& heap,
                            int num) {

    int n = prizes.size();
    const double epsilon = 1;

    while (num > 0 && !heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);
        num --;
        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = prizes[v]/ (minDist[v] + epsilon);
                heap.push({new_score, v});
            }
        } 
    }

}

//s-t path
std::vector<int> GCP_PrizeRatioTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx) {
    if (s_idx < 0 || s_idx >= (int)prizes.size() || budget <= 0) {
        return {s_idx};
    }

    std::vector<double> prizes_selected; 
    std::vector<std::vector<double>> costs_selected; 
    std::vector<int> indices;
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double epsilon = 1;

    selected[s_idx] = true;
    indices.push_back(s_idx);
    selected[t_idx] = true;
    indices.push_back(t_idx);

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (v != s_idx && v != t_idx ) {
            minDist[v] = std::min(costs[s_idx][v], costs[t_idx][v]);
            double score = prizes[v] / (minDist[v] + epsilon);
            heap.push({score, v});
        }
    }

    int N_pre = 2;
    int N = std::min(int(prizes.size()), int(2*N_pre));
    while(true) {
        if(N < 3) break;
        int num_to_select = N - N_pre;
        selectNodesByPrizeDistanceRatio(costs, prizes, indices, selected, minDist, heap, num_to_select);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);

        double total_mst_weight = computeMSTWeight(g, weight, mst_edges);
        if(total_mst_weight > budget) {
            break;
        }

        if (N >= (int)prizes.size())
            break;
        
        N_pre = N, N *= 2;
        if (N >= (int)prizes.size()) { 
            N = prizes.size();
        }

    }
    
    int N_min = std::max(2, N/3);
    int N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);

    while(best_mst_path.empty()) {
        N_max = N_min, N_min /= 2;
        binarySearchBestPathTwoFixed(costs, prizes,indices,best_mst_path, N_min, N_max, budget);
    }
    return best_mst_path;
}


double sumPrizes(const std::vector<double>& prizes, const std::vector<int>& path) {
    double total_prizes = 0.0;
    for (int node : path) {
        total_prizes += prizes[node];
    }
    return total_prizes;
}

std::vector<int> GCP(const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, 
                                    double budget, int s_idx, int t_idx) {

    // std::cerr << "running gcp";                                        
    auto f1 = std::async(std::launch::async, GCP_BestPrizeTwoFix,  costs, prizes, budget, s_idx, t_idx);
    auto f2 = std::async(std::launch::async, GCP_PrizeBoostTwoFix, costs, prizes, budget, s_idx, t_idx);
    auto f3 = std::async(std::launch::async, GCP_PrizeRatioTwoFix, costs, prizes, budget, s_idx, t_idx);

    std::vector<int> p1 = f1.get();
    std::vector<int> p2 = f2.get();
    std::vector<int> p3 = f3.get();

    for(auto p : {p1, p2, p3}) {
        bool is_improve = true;
        while(is_improve) {
            is_improve = false;
            is_improve |= two_opt_st(p, costs);
        }
        repairPath(costs, prizes, p, budget);
    }


    double s1 = sumPrizes(prizes, p1);
    double s2 = sumPrizes(prizes, p2);
    double s3 = sumPrizes(prizes, p3);

    if (s2 > s1 && s2 > s3) {
        return p2;
    } else if (s3 > s1 && s3 > s2) {
        return p3;
    } else {
        return p1;
    }
}

