#include "GCP_overlap.h"
#include "GCP.h"
#include "LocalSearch.h"
#include "PathOperations.h"

#include <algorithm>
#include <cmath>
#include <functional>
#include <future>
#include <limits>
#include <numeric>
#include <queue>
#include <stdexcept>
#include <vector>

// Greedy selection

constexpr double kEpsilon = 1.0;   // same denominator guard as the prize code

struct ScoredNode {
    double score;
    int slot;
    long long version;
    bool operator<(const ScoredNode& o) const { return score < o.score; }
};

double scoreOf(CoverageScore kind, double marginal, double min_dist) {
    switch (kind) {
        case CoverageScore::Marginal:
            return marginal;
        case CoverageScore::Boosted: {
            const double boosted = std::min(1.0, std::pow(marginal, 4));
            return std::asin(boosted) / (min_dist + kEpsilon);
        }
        case CoverageScore::Ratio:
        default:
            return marginal / (min_dist + kEpsilon);
    }
}

void selectNodesByCoverageInRange(CoverageOracle& oracle,
                                  std::vector<int>& indices,
                                  std::pair<int, int> range) {
    const int n = static_cast<int>(indices.size());
    if (range.first < 0 || range.first >= range.second || range.first >= n) return;
    if (range.second > n) range.second = n;

    // Seed from the prefix so repeated calls with growing ranges extend one
    // nested selection, which is what the doubling loop assumes.
    oracle.reset();
    for (int i = 0; i < range.first; ++i) oracle.commit(indices[i]);

    // Work on the tail by value rather than swapping in place; keeping heap
    // entries tied to positions that move under you is how the stale-high bug
    // sneaks back in.
    const std::vector<int> tail(indices.begin() + range.first, indices.end());
    std::vector<char> used(tail.size(), 0);

    // Lazy greedy 
    std::priority_queue<ScoredNode> heap;
    long long version = 0;
    for (std::size_t i = 0; i < tail.size(); ++i)
        heap.push({oracle.marginal(tail[i]), static_cast<int>(i), version});

    std::vector<int> picked;
    const int wanted = range.second - range.first;
    while (static_cast<int>(picked.size()) < wanted && !heap.empty()) {
        const ScoredNode top = heap.top();
        heap.pop();
        if (used[top.slot]) continue;
        if (top.version != version) {
            heap.push({oracle.marginal(tail[top.slot]), top.slot, version});
            continue;
        }
        used[top.slot] = 1;
        picked.push_back(top.slot);
        oracle.commit(tail[top.slot]);
        ++version;   // every cached score is now stale
    }

    std::vector<int> rebuilt;
    rebuilt.reserve(tail.size());
    for (int slot : picked) rebuilt.push_back(tail[slot]);
    for (std::size_t i = 0; i < tail.size(); ++i)
        if (!used[i]) rebuilt.push_back(tail[i]);
    std::copy(rebuilt.begin(), rebuilt.end(), indices.begin() + range.first);
}

void selectNodesByCoverageScore(const std::vector<std::vector<double>>& costs,
                                CoverageOracle& oracle,
                                std::vector<int>& indices,
                                std::vector<bool>& selected,
                                std::vector<double>& minDist,
                                int num,
                                CoverageScore score) {
    const int n = static_cast<int>(selected.size());
    if (num <= 0 || n == 0) return;
    if (static_cast<int>(costs.size()) != n || static_cast<int>(minDist.size()) != n)
        throw std::invalid_argument(
            "selectNodesByCoverageScore: costs/minDist size mismatch.");

    // The caller owns selected/minDist/indices across calls; the oracle's
    // committed set is rebuilt from `indices` so the two cannot drift apart.
    oracle.reset();
    oracle.commitAll(indices);

    std::priority_queue<ScoredNode> heap;
    long long version = 0;
    for (int v = 0; v < n; ++v)
        if (!selected[v])
            heap.push({scoreOf(score, oracle.marginal(v), minDist[v]), v, version});

    while (num > 0 && !heap.empty()) {
        const ScoredNode top = heap.top();
        heap.pop();
        if (selected[top.slot]) continue;
        if (top.version != version) {
            // Rescore rather than trust it. For Boosted and Ratio the
            // numerator falls while the denominator also falls, so the score
            // is not monotone in either direction and a stale value is simply
            // wrong rather than merely conservative.
            heap.push({scoreOf(score, oracle.marginal(top.slot), minDist[top.slot]),
                       top.slot, version});
            continue;
        }

        const int u = top.slot;
        selected[u] = true;
        indices.push_back(u);
        oracle.commit(u);
        --num;
        ++version;

        for (int v = 0; v < n; ++v) {
            if (selected[v]) continue;
            minDist[v] = std::min(minDist[v], costs[u][v]);
            heap.push({scoreOf(score, oracle.marginal(v), minDist[v]), v, version});
        }
    }
}

void selectNodesByCoverageBoosting(const std::vector<std::vector<double>>& costs,
                                   CoverageOracle& oracle,
                                   std::vector<int>& indices,
                                   std::vector<bool>& selected,
                                   std::vector<double>& minDist,
                                   int num) {
    selectNodesByCoverageScore(costs, oracle, indices, selected, minDist, num,
                               CoverageScore::Boosted);
}

void selectNodesByCoverageRatio(const std::vector<std::vector<double>>& costs,
                                CoverageOracle& oracle,
                                std::vector<int>& indices,
                                std::vector<bool>& selected,
                                std::vector<double>& minDist,
                                int num) {
    selectNodesByCoverageScore(costs, oracle, indices, selected, minDist, num,
                               CoverageScore::Ratio);
}


// Path scoring and insertion
double pathCoverage(const CoverageOracle& oracle, const std::vector<int>& path) {
    return oracle.rewardOf(path);
}

void insertHighestMarginalCheapest(const std::vector<std::vector<double>>& costs,
                                   const CoverageOracle& oracle,
                                   std::vector<int>& path,
                                   double& path_cost,
                                   double budget,
                                   std::vector<bool>& visited) {
    if (path.size() < 2) return;

    while (true) {
        int best_node = -1;
        int best_position = -1;
        double best_marginal = 0.0;
        double best_extra = 0.0;

        for (int v = 0; v < static_cast<int>(visited.size()); ++v) {
            if (visited[v]) continue;

            // Marginal over the path as it stands, so a tile that only repeats
            // pixels already covered scores zero and is never inserted.
            const double gain = oracle.marginalOver(v, path);
            if (gain <= 0.0 || gain < best_marginal) continue;

            // Cheapest interior splice; positions 1..size-1 keep both endpoints.
            double cheapest = std::numeric_limits<double>::infinity();
            int position = -1;
            for (std::size_t i = 1; i < path.size(); ++i) {
                const double extra = costs[path[i - 1]][v] + costs[v][path[i]] -
                                     costs[path[i - 1]][path[i]];
                if (extra < cheapest) {
                    cheapest = extra;
                    position = static_cast<int>(i);
                }
            }
            if (position < 0 || path_cost + cheapest > budget) continue;

            if (gain > best_marginal ||
                (gain == best_marginal && cheapest < best_extra)) {
                best_marginal = gain;
                best_node = v;
                best_position = position;
                best_extra = cheapest;
            }
        }

        if (best_node < 0) break;
        path.insert(path.begin() + best_position, best_node);
        visited[best_node] = true;
        path_cost += best_extra;
    }
}


void binarySearchBestCoveragePathTwoFixed(
        const std::vector<std::vector<double>>& costs,
        const CoverageOracle& oracle,
        const std::vector<int>& indices,
        std::vector<int>& best_mst_path,
        int N_min, int N_max, double budget) {

    double best_coverage = 0.0;
    N_max += 1;
    while (N_max - N_min > 1) {
        const int N = (N_max + N_min) / 2;

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        if (computeMSTWeight(g, weight, mst_edges) > budget) {
            N_max = N;
            continue;
        }

        if (N == 2) {
            const double direct = costs[indices[0]][indices[1]];
            if (direct <= budget) best_mst_path = {indices[0], indices[1]};
            else                  best_mst_path = {indices[0]};
            return;
        }

        std::vector<int> mst_path =
            christofidesPathTwoFixed(g, costs, indices, N, mst_edges, 0, 1);
        std::vector<int> path;
        path.reserve(mst_path.size());
        for (int node : mst_path) path.push_back(indices[node]);

        while (two_opt_st(path, costs)) {}

        double path_cost = 0.0;
        for (std::size_t i = 1; i < path.size(); ++i)
            path_cost += costs[path[i - 1]][path[i]];

        if (path_cost > budget) {
            N_max = N;
            continue;
        }

        std::vector<bool> visited(costs.size(), false);
        for (int v : path) visited[v] = true;
        insertHighestMarginalCheapest(costs, oracle, path, path_cost, budget, visited);

        // Union coverage
        const double coverage = oracle.rewardOf(path);
        if (coverage > best_coverage) {
            best_coverage = coverage;
            best_mst_path = path;
        }
        N_min = N;
    }
}

std::vector<int> GCP_BestCoverageTwoFix(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        double budget, int s_idx, int t_idx) {

    const int n = static_cast<int>(costs.size());
    if (s_idx < 0 || s_idx >= n || t_idx < 0 || t_idx >= n ||
        s_idx == t_idx || budget <= 0)
        return {s_idx};
    if (costs[s_idx][t_idx] > budget) return {};

    CoverageOracle oracle(pixel_probs, member_pixels);

    std::vector<int> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    // Guard the swap order: if t_idx were 0 the second swap would clobber s.
    std::swap(indices[0], indices[s_idx]);
    std::swap(indices[1], indices[(t_idx == 0) ? s_idx : t_idx]);

    int N_pre = 2;
    int N = std::min(n, 2 * N_pre);
    while (true) {
        if (N < 3) break;

        selectNodesByCoverageInRange(oracle, indices, {N_pre, N});

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        if (computeMSTWeight(g, weight, mst_edges) > budget) break;
        if (N >= n) break;

        N_pre = N;
        N = std::min(n, N * 2);
    }

    int N_min = std::max(2, N / 3);
    int N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestCoveragePathTwoFixed(costs, oracle, indices, best_mst_path,
                                         N_min, N_max, budget);

    while (best_mst_path.empty()) {
        if (N_max <= 2) return {s_idx, t_idx};
        N_max = N_min;
        N_min = std::max(2, N_min / 2);
        binarySearchBestCoveragePathTwoFixed(costs, oracle, indices, best_mst_path,
                                             N_min, N_max, budget);
        if (N_min == N_max) return {s_idx, t_idx};
    }
    return best_mst_path;
}

// shared by GCP_CoverageBoostTwoFix and GCP_CoverageRatioTwoFix
std::vector<int> gcpCoverageScoredTwoFix(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        double budget, int s_idx, int t_idx, CoverageScore score) {

    const int n = static_cast<int>(costs.size());
    if (s_idx < 0 || s_idx >= n || t_idx < 0 || t_idx >= n ||
        s_idx == t_idx || budget <= 0)
        return {s_idx};
    if (costs[s_idx][t_idx] > budget) return {};

    CoverageOracle oracle(pixel_probs, member_pixels);

    std::vector<int> indices;
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());

    selected[s_idx] = true;
    indices.push_back(s_idx);
    selected[t_idx] = true;
    indices.push_back(t_idx);

    // Keep the minDist seeding from the prize version. The heap that used to
    // be seeded alongside it is gone: it is rebuilt inside the selector, since
    // a heap carried across calls is what lets a stale-high marginal survive.
    for (int v = 0; v < n; ++v)
        if (v != s_idx && v != t_idx)
            minDist[v] = std::min(costs[s_idx][v], costs[t_idx][v]);

    int N_pre = 2;
    int N = std::min(n, 2 * N_pre);
    while (true) {
        if (N < 3) break;

        selectNodesByCoverageScore(costs, oracle, indices, selected, minDist,
                                   N - N_pre, score);

        Graph g;
        Graph::EdgeMap<double> weight(g);
        buildGraphUpdate(costs, indices, N, g, weight);

        std::vector<Edge> mst_edges;
        buildMST(g, weight, mst_edges);
        if (computeMSTWeight(g, weight, mst_edges) > budget) break;
        if (N >= n) break;

        N_pre = N;
        N = std::min(n, N * 2);
    }

    int N_min = std::max(2, N / 3);
    int N_max = N;
    std::vector<int> best_mst_path;
    binarySearchBestCoveragePathTwoFixed(costs, oracle, indices, best_mst_path,
                                         N_min, N_max, budget);

    while (best_mst_path.empty()) {
        if (N_max <= 2) return {s_idx, t_idx};
        N_max = N_min;
        N_min = std::max(2, N_min / 2);
        binarySearchBestCoveragePathTwoFixed(costs, oracle, indices, best_mst_path,
                                             N_min, N_max, budget);
        if (N_min == N_max) return {s_idx, t_idx};
    }
    return best_mst_path;
}

std::vector<int> GCP_CoverageBoostTwoFix(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        double budget, int s_idx, int t_idx) {
    return gcpCoverageScoredTwoFix(costs, pixel_probs, member_pixels, budget,
                                   s_idx, t_idx, CoverageScore::Boosted);
}

std::vector<int> GCP_CoverageRatioTwoFix(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        double budget, int s_idx, int t_idx) {
    return gcpCoverageScoredTwoFix(costs, pixel_probs, member_pixels, budget,
                                   s_idx, t_idx, CoverageScore::Ratio);
}

std::vector<int> GCP_Coverage(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& pixel_probs,
        const std::vector<std::vector<int>>& member_pixels,
        double budget, int s_idx, int t_idx) {

    auto f1 = std::async(std::launch::async, GCP_BestCoverageTwoFix,
                         std::cref(costs), std::cref(pixel_probs),
                         std::cref(member_pixels), budget, s_idx, t_idx);
    auto f2 = std::async(std::launch::async, GCP_CoverageBoostTwoFix,
                         std::cref(costs), std::cref(pixel_probs),
                         std::cref(member_pixels), budget, s_idx, t_idx);
    auto f3 = std::async(std::launch::async, GCP_CoverageRatioTwoFix,
                         std::cref(costs), std::cref(pixel_probs),
                         std::cref(member_pixels), budget, s_idx, t_idx);

    std::vector<std::vector<int>> paths{f1.get(), f2.get(), f3.get()};

    CoverageOracle oracle(pixel_probs, member_pixels);
    std::vector<int> best;
    double best_coverage = -1.0;
    for (auto& p : paths) {
        if (p.size() < 2) continue;

        // Every GCP constructor already returns a budget-feasible s-t path.
        // Two-opt only changes the order and reduces its travel cost, so no
        // separate PixelOverlapPath repair helper is required here.
        while (two_opt_st(p, costs)) {}

        const double coverage = oracle.rewardOf(p);
        if (coverage > best_coverage) {
            best_coverage = coverage;
            best = p;
        }
    }
    return best;
}
