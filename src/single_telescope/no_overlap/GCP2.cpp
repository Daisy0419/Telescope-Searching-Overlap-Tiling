#include "GCP.h"
#include "ReadData.h"
#include "helplers.h"
#include "LocalSearch.h"
#include "PathOperations.h"

#include <vector>
#include <queue>
#include <tuple>
#include <limits>
#include <algorithm>
#include <cmath>
#include <unordered_set>


// Budget-based k-NN (direct edges)
std::vector<std::vector<int>> knn_budget_edges(const std::vector<std::vector<double>>& costs, double budget) {
    const int n = (int)costs.size();

    std::vector<std::vector<int>> knn(n);
    for (int i = 0; i < n; ++i) {
        std::vector<int> idx; idx.reserve(n-1);
        for (int j = 0; j < n; ++j) if (j != i) idx.push_back(j);

        std::sort(idx.begin(), idx.end(),
                  [&](int a, int b){ return costs[i][a] < costs[i][b]; });

        double sum = 0.0;
        for (int j : idx) {
            double w = costs[i][j];
            if (sum + w > budget) break;
            sum += w;
            knn[i].push_back(j);
        }
    }
    return knn;
}

std::vector<std::vector<int>> ratio_budget_neighbors(const std::vector<std::vector<double>>& costs,
                                                const std::vector<double>& prizes,
                                                double budget) {
    const int n = (int)costs.size();
    std::vector<std::vector<int>> picks(n);
    if (n == 0) return picks;
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    for (int i = 0; i < n; ++i)
        if ((int)costs[i].size() != n)
            throw std::runtime_error("costs must be square");

    const double EPS = 1e-12;

    struct Cand {
        int j;
        double ratio;
        double prize;
        double cost;
    };

    for (int i = 0; i < n; ++i) {
        // Build candidates with ratios relative to fixed node i
        std::vector<Cand> cand;
        cand.reserve(n-1);

        for (int j = 0; j < n; ++j) {
            if (j == i) continue;
            double w = costs[i][j];
            if (!std::isfinite(w)) continue;          // skip invalid edges

            double r;
            if (w <= EPS) {
                // Zero/near-zero cost: treat positive prize as +inf ratio
                if (prizes[j] > 0.0)      r = std::numeric_limits<double>::infinity();
                else if (prizes[j] < 0.0) r = -std::numeric_limits<double>::infinity();
                else                       r = 0.0;
            } else {
                r = prizes[j] / w;
            }
            cand.push_back(Cand{j, r, prizes[j], w});
        }

        // Sort by: higher ratio, then higher prize, then lower cost, then lower index
        std::sort(cand.begin(), cand.end(), [](const Cand& a, const Cand& b){
            if (a.ratio != b.ratio)     return a.ratio > b.ratio;
            if (a.prize != b.prize)     return a.prize > b.prize;
            if (a.cost  != b.cost)      return a.cost  < b.cost;
            return a.j < b.j;
        });

        double spent = 0.0;
        auto& out = picks[i];

        // Greedily include feasible items in that sorted order
        for (const auto& c : cand) {
            if (spent + c.cost <= budget + EPS) {
                out.push_back(c.j);
                spent += c.cost;
                // continue scanning—lower ratio items might still fit
            }
        }
        // Note: if nothing fits, out is empty for node i (analogous to {} in KNN).
    }

    return picks;
}


// Heuristic score: attach `node` to `cluster` paying min_dist,
// then greedily “fill” with node’s KNN neighbors until that cluster’s budget runs out.
// Only used for ranking; we actually add a single node when we pop the heap.
double getScore(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                const std::vector<bool>& selected,
                const std::vector<std::vector<int>>& knn,
                const std::vector<double>& remain_budget,
                int node, int cluster, double min_dist) {

    if (min_dist > remain_budget[cluster]) return -1.0;

    double sum_prize = prizes[node];
    double sum_cost  = min_dist;

    for (int nb : knn[node]) {
        if (selected[nb]) continue;
        double w = costs[node][nb];
        if (sum_cost + w <= remain_budget[cluster]) {
            sum_cost  += w;
            sum_prize += prizes[nb];
        } else break;
    }

    // prefer higher prize per unit cost; guard zero
    // return (sum_cost > 0.0) ? (sum_prize / sum_cost) : std::numeric_limits<double>::infinity();
    return (min_dist > 0.0) ? (sum_prize / min_dist) : std::numeric_limits<double>::infinity();
}

double getScore2(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                const std::vector<bool>& selected,
                const std::vector<std::vector<int>>& knn,
                const std::vector<double>& remain_budget,
                int node, int cluster, double min_dist) {

    if (min_dist > remain_budget[cluster]) return -1.0;

    double sum_prize = prizes[node];
    double sum_cost  = min_dist;
    std::vector<bool> local_selected =selected;

    int cur_node = node;
    double w = 0.0;
    while (sum_cost <= remain_budget[cluster]) {
        bool moved = false;
        for (int nb : knn[cur_node]) {
            if (local_selected[nb]) continue;
            w = costs[cur_node][nb];
            if (sum_cost + w <= remain_budget[cluster]) {
                sum_cost  += w;
                sum_prize += prizes[nb];
                local_selected[nb] = true;
                cur_node = nb;
                moved = true;
                break;  // only break after successfully moving
            }
        }
        if (!moved) break; // avoid infinite loop if stuck
    }

    // prefer higher prize per unit cost; guard zero
    // return (sum_cost > 0.0) ? (sum_prize / sum_cost) : std::numeric_limits<double>::infinity();
    // return (min_dist > 0.0) ? (sum_prize / min_dist) : std::numeric_limits<double>::infinity();
    return (min_dist > 0.0) ? (sum_prize) : std::numeric_limits<double>::infinity();
}

struct Entry {
    double score;
    double stamp_budget; // budget used when computing `score`
    double stamp_dist;   // minDist used when computing `score`
    int v;
    int cluster;
};

struct EntryLess {
    bool operator()(const Entry& a, const Entry& b) const {
        return a.score < b.score; // max-heap on score
    }
};


// Returns the clusters (each begins with s,t). Nodes are unique across clusters.
// Budget is per-cluster. Greedy best-first by the score above.
std::vector<std::vector<int>>
selectNClusters(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                const std::vector<std::vector<int>>& knn,
                int num_paths, int s, int t, double budget) {

    // validate_symmetric(costs);
    const int n = (int)costs.size();
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    if (s < 0 || s >= n || t < 0 || t >= n || s == t) throw std::runtime_error("bad s/t");
    if (num_paths <= 0) return {};

    const double EPS = 1e-12;

    std::vector<std::vector<int>> clusters(num_paths);
    std::vector<double> remain_budget(num_paths, budget-costs[s][t]);
    std::vector<bool> selected(n, false);

    // seed clusters with s,t; never selectable again
    selected[s] = true;
    selected[t] = true;

    // minDist[c][v] = current connection cost to add v to cluster c
    std::vector<std::vector<double>> minDist(num_paths, std::vector<double>(n, std::numeric_limits<double>::infinity()));
    for (int c = 0; c < num_paths; ++c) {
        clusters[c].push_back(s);
        clusters[c].push_back(t);
        for (int v = 0; v < n; ++v) {
            if (v == s || v == t) continue;
            minDist[c][v] = std::min(costs[s][v], costs[t][v]);
        }
    }

    std::priority_queue<Entry, std::vector<Entry>, EntryLess> pq;

    auto maybe_push = [&](int v, int c) {
        if (selected[v]) return;
        double d = minDist[c][v];
        if (d > remain_budget[c]) return;
        double sc = getScore2(costs, prizes, selected, knn, remain_budget, v, c, d);
        if (sc > 0.0) {
            pq.push(Entry{sc, remain_budget[c], d, v, c});
        }
    };

    // initialize heap
    for (int v = 0; v < n; ++v) {
        if (v == s || v == t) continue;
        for (int c = 0; c < num_paths; ++c) maybe_push(v, c);
    }

    // Greedy extraction
    while (!pq.empty()) {
        Entry e = pq.top(); pq.pop();

        if (selected[e.v]) continue;                         // already taken elsewhere
        if (e.stamp_budget < remain_budget[e.cluster] - EPS  // budget changed (got bigger)
            || e.stamp_budget > remain_budget[e.cluster] + EPS) {
            // budget changed – recompute a fresh entry
            maybe_push(e.v, e.cluster);
            continue;
        }
        // // If connection distance improved since we pushed, recompute with cheaper cost
        // if (minDist[e.cluster][e.v] + 1e-12 < e.stamp_dist) {
        //     maybe_push(e.v, e.cluster);
        //     continue;
        // }

        double need = minDist[e.cluster][e.v];
        if (need > remain_budget[e.cluster] + EPS) continue; // stale

        // Commit selection
        selected[e.v] = true;
        clusters[e.cluster].push_back(e.v);
        remain_budget[e.cluster] -= need;

        // Update that cluster's minDist and push new candidates
        for (int w = 0; w < n; ++w) {
            if (selected[w]) continue;
            double d = costs[e.v][w];
            if (d < minDist[e.cluster][w]) {
                minDist[e.cluster][w] = d;
            }
            maybe_push(w, e.cluster);
        }
    }

    return clusters;
}

// Variant: per-cluster starts (s_c for cluster c), shared t.
std::vector<std::vector<int>>
selectNClustersMultiS(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      const std::vector<std::vector<int>>& knn,
                      const std::vector<int>& starts,
                      int t, double budget) {

    // using std::vector;
    const int n = (int)costs.size();
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    if (n == 0) return {};

    const int num_paths = (int)starts.size();
    if (num_paths <= 0) return {};
    if (t < 0 || t >= n) throw std::runtime_error("bad t");

    // // Validate 'starts': in range, distinct, and != t
    // {
    //     std::vector<char> seen(n, 0);
    //     for (int c = 0; c < num_paths; ++c) {
    //         int s_c = starts[c];
    //         if (s_c < 0 || s_c >= n) throw std::runtime_error("bad s for some cluster");
    //         if (s_c == t) throw std::runtime_error("s_c must differ from t");
    //         if (seen[s_c]) throw std::runtime_error("duplicate s across clusters");
    //         seen[s_c] = 1;
    //     }
    // }

    const double EPS = 1e-12;

    std::vector<std::vector<int>> clusters(num_paths);
    std::vector<double> remain_budget(num_paths, 0.0);
    std::vector<bool> selected(n, false);

    // mark all starts and t as "taken" so they are not chosen later
    for (int c = 0; c < num_paths; ++c) selected[starts[c]] = true;
    selected[t] = true;

    // minDist[c][v] = current connection cost to add v to cluster c
    std::vector<std::vector<double>> minDist(num_paths, std::vector<double>(n, std::numeric_limits<double>::infinity()));

    // Seed each cluster with its own s_c and shared t
    for (int c = 0; c < num_paths; ++c) {
        int s_c = starts[c];
        clusters[c].push_back(s_c);
        clusters[c].push_back(t);

        // Baseline "must-pay" is s_c -> t (so routes can at least connect)
        double base = costs[s_c][t];
        if (budget + EPS < base)
            throw std::runtime_error("budget too small for some s_c->t");

        remain_budget[c] = budget - base;

        for (int v = 0; v < n; ++v) {
            if (v == s_c || v == t) continue;
            // Initial attach cost is to the nearer of s_c or t
            minDist[c][v] = std::min(costs[s_c][v], costs[t][v]);
        }
    }

    std::priority_queue<Entry, std::vector<Entry>, EntryLess> pq;

    auto maybe_push = [&](int v, int c) {
        if (selected[v]) return;
        double d = minDist[c][v];
        if (d > remain_budget[c]) return;
        double sc = getScore(costs, prizes, selected, knn, remain_budget, v, c, d);
        if (sc > 0.0) {
            pq.push(Entry{sc, remain_budget[c], d, v, c});
        }
    };

    // initialize heap
    for (int v = 0; v < n; ++v) {
        if (selected[v]) continue; // skips all s_c and t
        for (int c = 0; c < num_paths; ++c) maybe_push(v, c);
    }

    // Greedy extraction 
    while (!pq.empty()) {
        Entry e = pq.top(); pq.pop();

        if (selected[e.v]) continue; // already taken elsewhere

        // budget changed?
        if (e.stamp_budget < remain_budget[e.cluster] - EPS ||
            e.stamp_budget > remain_budget[e.cluster] + EPS) {
            maybe_push(e.v, e.cluster);
            continue;
        }

        double need = minDist[e.cluster][e.v];
        if (need > remain_budget[e.cluster] + EPS) continue; // stale

        // Commit selection into cluster e.cluster
        selected[e.v] = true;
        clusters[e.cluster].push_back(e.v);
        remain_budget[e.cluster] -= need;

        // Update that cluster's frontier and push new candidates for it
        for (int w = 0; w < n; ++w) {
            if (selected[w]) continue;
            double d = costs[e.v][w];
            if (d < minDist[e.cluster][w]) {
                minDist[e.cluster][w] = d;
            }
            maybe_push(w, e.cluster);
        }
    }

    return clusters;
}


std::vector<int> buildPath(const std::vector<std::vector<double>>& costs, 
                    const std::vector<int>& indices,
                    double budget) {
    Graph g;
    Graph::EdgeMap<double> weight(g);
    buildGraphUpdate(costs, indices, indices.size(), g, weight);

    std::vector<Edge> mst_edges;
    buildMST(g, weight, mst_edges);
    std::vector<int> mst_path = std::move(christofidesPathTwoFixed(g, costs, indices, indices.size(), mst_edges, 0, 1));
    std::vector<int> mst_path_convert;
    for (int node : mst_path) {
        mst_path_convert.push_back(indices[node]);
    }       
    
    //2-opt to remove cross
    // two_opt_st(mst_path_convert, costs);
    bool has_improve = two_opt_st(mst_path_convert, costs);
    while(has_improve) {
        // has_improve = false;
        has_improve = two_opt_st(mst_path_convert, costs);
    }
    return mst_path_convert;
}


std::pair<std::vector<std::vector<int>>, double> GCP2(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes, 
                    double budget, int s, int t, int num_paths) {
    // std::vector<std::vector<int>> knn = std::move(knn_budget_edges(costs, budget*num_paths));
    std::vector<std::vector<int>> knn = std::move(ratio_budget_neighbors(costs, prizes, budget));
    std::vector<std::vector<int>> clusters = selectNClusters(costs,prizes,knn,
                                                            num_paths, s, t, budget);
    // print_paths(clusters);
    std::vector<std::vector<int>> paths;                                             
    for(auto cluster : clusters) {
        std::vector<int> path = buildPath(costs, cluster, budget);
        // fix_cross_st_path(path, costs);
        bool has_improve = two_opt_st(path, costs);
        while(has_improve) {
            // has_improve = false;
            has_improve = two_opt_st(path, costs);
        }
        paths.push_back(path);
    }
    // print_paths(paths);
    repairMultiPaths(costs, prizes, paths, budget, s, t);
    // print_paths(paths);

    double total_prize = multiPathsSumPrize(prizes, paths);
    return {paths, total_prize};
}

std::pair<std::vector<std::vector<int>>, double> GCP2(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes, 
                    double budget, const std::vector<int>& s, int t, int num_paths) {
    std::vector<std::vector<int>> knn = std::move(knn_budget_edges(costs, budget*num_paths));
    // std::vector<std::vector<int>> knn = std::move(ratio_budget_neighbors(costs, prizes, budget*num_paths));
    std::vector<std::vector<int>> clusters = selectNClustersMultiS(costs,prizes,knn,
                                                                    s, t, budget);
    // print_paths(clusters);
    std::vector<std::vector<int>> paths;                                             
    for(auto cluster : clusters) {
        std::vector<int> path = buildPath(costs, cluster, budget);

        bool has_improve = two_opt_st(path, costs);
        while(has_improve) {
            // has_improve = false;
            has_improve = two_opt_st(path, costs);
        }
        paths.push_back(path);
    }
    // print_paths(paths);
    repairMultiPaths(costs, prizes, paths, budget, s, t);
    // print_paths(paths);

    double total_prize = multiPathsSumPrize(prizes, paths);
    return {paths, total_prize};
}