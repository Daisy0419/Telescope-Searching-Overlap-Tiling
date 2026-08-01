#include "ReadData.h"
#include "PathSplitter.h"
#include "PathOperations.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>
#include <iostream>
#include <tuple>
#include <functional>
#include <unordered_set>



double route_cost(const Instance& I, const std::vector<int>& path) {
    double C = 0.0;
    for (size_t i = 0; i + 1 < path.size(); ++i) C += I.costs[path[i]][path[i+1]];
    return C;
}
double route_prize(const Instance& I, const std::vector<int>& path) {
    double P = 0.0;
    for (int v : path) if (v != I.s && v != I.t) P += I.prizes[v];
    return P;
}


// 1) Single-route best subsequence DP (budget discretized).
//    Input : I, B, perm (giant route without s,t), bins
//    Output: (feasible, best_route [s..t], true_cost, total_prize)

int to_bin(double c, double scale) {
    if (c <= 0.0) return 0;
    double z = std::ceil(c * scale - 1e-12);
    if (z < 0.0) z = 0.0;
    if (z > static_cast<double>(std::numeric_limits<int>::max()/4))
        z = static_cast<double>(std::numeric_limits<int>::max()/4);
    return static_cast<int>(z);
}


std::tuple<bool, std::vector<int>, double, double>
best_subsequence_dp(const Instance& I, double B,
                    const std::vector<int>& perm,
                    int bins) {
    // Build W = [s, perm..., t]
    std::vector<int> W; W.reserve(perm.size() + 2);
    W.push_back(I.s);
    W.insert(W.end(), perm.begin(), perm.end());
    W.push_back(I.t);
    const int m = (int)W.size();
    if (m < 2) return {false, {}, 0.0, 0.0};

    const int K = std::max(1, bins);
    const double scale = (B > 0.0) ? (double)K / B : 1e12;

    // Bin costs for forward jumps i->j
    std::vector<std::vector<int>> Cbin(m, std::vector<int>(m, std::numeric_limits<int>::max()/4));
    for (int i = 0; i < m; ++i)
        for (int j = i + 1; j < m; ++j)
            Cbin[i][j] = to_bin(I.costs[W[i]][W[j]], scale);

    // Reward on arrival
    std::vector<double> R(m, 0.0);
    for (int j = 1; j + 1 < m; ++j) R[j] = I.prizes[W[j]];

    const double NEG = -1e300;
    std::vector<std::vector<double>> F(m, std::vector<double>(K+1, NEG));
    std::vector<std::vector<int>>    Prev(m, std::vector<int>(K+1, -1));
    std::vector<std::vector<int>>    PrevB(m, std::vector<int>(K+1, -1));
    for (int b = 0; b <= K; ++b) F[0][b] = 0.0;

    for (int j = 1; j < m; ++j) {
        for (int b = 0; b <= K; ++b) {
            double best = NEG; int best_i = -1, best_bp = -1;
            for (int i = 0; i < j; ++i) {
                int cb = Cbin[i][j];
                if (cb <= b) {
                    double prev = F[i][b - cb];
                    if (prev > NEG/2) {
                        double cand = prev + R[j];
                        if (cand > best) { best = cand; best_i = i; best_bp = b - cb; }
                    }
                }
            }
            F[j][b] = best; Prev[j][b] = best_i; PrevB[j][b] = best_bp;
        }
    }
    // pick best arrival at t
    int best_b = 0; double best_val = NEG;
    for (int b = 0; b <= K; ++b)
        if (F[m-1][b] > best_val) { best_val = F[m-1][b]; best_b = b; }
    if (best_val <= NEG/2) return {false, {}, 0.0, 0.0};

    // backtrack
    std::vector<int> nodes;
    for (int j = m-1, b = best_b; j >= 0; ) {
        nodes.push_back(j);
        int pi = Prev[j][b], pb = PrevB[j][b];
        if (pi < 0) break;
        j = pi; b = pb;
    }
    std::reverse(nodes.begin(), nodes.end());

    std::vector<int> route; route.reserve(nodes.size());
    for (int idx : nodes) route.push_back(W[idx]);
    double Ctrue = route_cost(I, route);
    double Ptrue = route_prize(I, route);
    bool feasible = (Ctrue <= B + 1e-12);
    return {feasible, route, Ctrue, Ptrue};
}


// team split: split into ≤ m disjoint INDEX RANGES, each range uses the
// best-subsequence DP inside it (so the resulting route can skip vertices).
// Complexity: Precompute O(n^2) subranges with a single-route DP (bins);
// then the same O(m n^2) outer DP .
std::pair<std::vector<std::vector<int>>, double>
split_segment_dp_with_subseq(const Instance& I, double B,
                             const std::vector<int>& perm,
                             int m_routes,
                             int bins_for_subseq) {
    const int n = (int)perm.size();
    std::vector<std::vector<int>> routes;
    double total_prize = 0.0;
    if (n == 0 || m_routes <= 0) return {routes, total_prize};

    // Precompute, for every subrange [i..j], the best subsequence route & prize.
    // i,j are 1-based indices into 'perm' for consistency with your code.
    std::vector<std::vector<char>>   ok(n+1, std::vector<char>(n+1, 0));
    std::vector<std::vector<double>> bestPrize(n+1, std::vector<double>(n+1, 0.0));
    // Store the actual precomputed route so we can output it directly on backtrack.
    std::vector<std::vector<std::vector<int>>> bestRoute(n+1,
        std::vector<std::vector<int>>(n+1));

    for (int i = 1; i <= n; ++i) {
        std::vector<int> sub;
        sub.reserve(n - i + 1);
        for (int j = i; j <= n; ++j) {
            sub.push_back(perm[j-1]);  // grow subrange [i..j]
            auto [feas, r, cst, pr] = best_subsequence_dp(I, B, sub, int(bins_for_subseq/m_routes));
            ok[i][j] = feas ? 1 : 0;
            if (feas) {
                bestPrize[i][j] = pr;
                bestRoute[i][j] = std::move(r); // r is [s, ... , t]
            } else {
                bestPrize[i][j] = 0.0;
                bestRoute[i][j].clear();
            }
        }
    }

    // Outer DP: same as segment DP but uses bestPrize/ok for each [i..j]
    const double NEG = -1e300;
    std::vector<std::vector<double>> DP(m_routes+1, std::vector<double>(n+1, NEG));
    for (int j = 0; j <= n; ++j) DP[0][j] = 0.0;

    for (int k = 1; k <= m_routes; ++k) {
        DP[k][0] = 0.0;
        for (int j = 1; j <= n; ++j) {
            double best = DP[k][j-1]; // skip j
            for (int i = 1; i <= j; ++i) {
                if (!ok[i][j]) continue;
                double cand = DP[k-1][i-1] + bestPrize[i][j];
                if (cand > best) best = cand;
            }
            DP[k][j] = best;
        }
    }

    // Backtrack to recover chosen disjoint ranges (and their internal subsequences)
    int k = m_routes, j = n;
    auto equal_eps = [](double a,double b){
        return std::abs(a-b) <= 1e-9 * (1.0 + std::max(std::abs(a),std::abs(b)));
    };
    std::vector<std::pair<int,int>> chosen;
    while (k > 0 && j > 0) {
        if (equal_eps(DP[k][j], DP[k][j-1])) { --j; continue; }
        int best_i = -1;
        for (int i = 1; i <= j; ++i) {
            if (!ok[i][j]) continue;
            if (equal_eps(DP[k][j], DP[k-1][i-1] + bestPrize[i][j])) { best_i = i; break; }
        }
        if (best_i < 0) { --j; continue; }
        chosen.emplace_back(best_i, j);
        j = best_i - 1; --k;
    }
    std::reverse(chosen.begin(), chosen.end());

    // Build final routes & total prize from the precomputed best routes of each [i..j]
    for (auto [i, j2] : chosen) {
        const std::vector<int>& r = bestRoute[i][j2]; // already [s, ... , t]
        routes.push_back(r);
        total_prize += route_prize(I, r);
    }
    return {routes, total_prize};
}


// Team split – contiguous segment DP (O(m·n²))
//    Output: (routes, total_prize)
std::pair<std::vector<std::vector<int>>, double>
split_segment_dp(const Instance& I, double B,
                 const std::vector<int>& perm,
                 int m_routes) {

    const int n = (int)perm.size();
    std::vector<std::vector<int>> routes;
    double total_prize = 0.0;
    if (n == 0 || m_routes <= 0) return {routes, total_prize};

    // precompute feasibility and segment prizes
    std::vector<std::vector<char>> feasible(n+1, std::vector<char>(n+1, 0));
    std::vector<std::vector<double>> segPrize(n+1, std::vector<double>(n+1, 0.0));
    std::vector<double> pref(n+1, 0.0);
    for (int i = 1; i <= n; ++i) pref[i] = pref[i-1] + I.prizes[perm[i-1]];
    auto prize_ij = [&](int i, int j){ return pref[j] - pref[i-1]; };

    for (int i = 1; i <= n; ++i) {
        double acc = I.costs[I.s][perm[i-1]];
        for (int j = i; j <= n; ++j) {
            if (j > i) acc += I.costs[perm[j-2]][perm[j-1]];
            double c = acc + I.costs[perm[j-1]][I.t];
            feasible[i][j] = (c <= B + 1e-12) ? 1 : 0;
            segPrize[i][j] = prize_ij(i,j);
        }
    }

    const double NEG = -1e300;
    std::vector<std::vector<double>> DP(m_routes+1, std::vector<double>(n+1, NEG));
    for (int j = 0; j <= n; ++j) DP[0][j] = 0.0;

    for (int k = 1; k <= m_routes; ++k) {
        DP[k][0] = 0.0;
        for (int j = 1; j <= n; ++j) {
            double best = DP[k][j-1]; // skip j
            for (int i = 1; i <= j; ++i) {
                if (!feasible[i][j]) continue;
                double cand = DP[k-1][i-1] + segPrize[i][j];
                if (cand > best) best = cand;
            }
            DP[k][j] = best;
        }
    }

    // backtrack segments
    int k = m_routes, j = n;
    auto equal_eps = [](double a,double b){
        return std::abs(a-b) <= 1e-9*(1.0 + std::max(std::abs(a),std::abs(b)));
    };
    std::vector<std::pair<int,int>> segs;
    while (k > 0 && j > 0) {
        if (equal_eps(DP[k][j], DP[k][j-1])) { --j; continue; }
        int best_i = -1;
        for (int i = 1; i <= j; ++i) {
            if (!feasible[i][j]) continue;
            if (equal_eps(DP[k][j], DP[k-1][i-1] + segPrize[i][j])) { best_i = i; break; }
        }
        if (best_i < 0) { --j; continue; }
        segs.emplace_back(best_i, j);
        j = best_i - 1; --k;
    }
    std::reverse(segs.begin(), segs.end());

    for (auto [i, j2] : segs) {
        std::vector<int> r; r.reserve(j2 - i + 3);
        r.push_back(I.s);
        for (int idx = i; idx <= j2; ++idx) r.push_back(perm[idx-1]);
        r.push_back(I.t);
        total_prize += route_prize(I, r);
        routes.push_back(std::move(r));
    }
    return {routes, total_prize};
}


// Team split – interval-graph DP (Dang–Guibadj–Moukrim, O(m·n))
//    Output: (routes, total_prize)
std::pair<std::vector<std::vector<int>>, double>
split_interval_graph(const Instance& I, double B,
                     const std::vector<int>& perm,
                     int m_routes) {

    const int n = (int)perm.size();
    std::vector<std::vector<int>> routes;
    if (n == 0 || m_routes <= 0) return {routes, 0.0};

    // (1) longest feasible prefix (saturated tour) for each start i
    std::vector<int>    lmax(n, -1);      // -1 means: no feasible tour from i
    std::vector<int>    succ(n, n);       // index in [0..n], n = past-the-end
    std::vector<double> P(n, 0.0);

    for (int i = 0; i < n; ++i) {
        double acc = I.costs[I.s][perm[i]];
        int j = i, last = i - 1;          // last feasible index in [i..n-1]
        while (j < n && acc + I.costs[perm[j]][I.t] <= B + 1e-12) {
            last = j;
            ++j;
            if (j < n) acc += I.costs[perm[j-1]][perm[j]];
        }
        if (last >= i) {
            lmax[i] = last - i;           // length >= 0
            succ[i] = last + 1;           // ∈ [i+1 .. n]
            double p = 0.0; for (int k = i; k <= last; ++k) p += I.prizes[perm[k]];
            P[i] = p;
        } else {
            // no feasible saturated tour starting at i
            lmax[i] = -1;
            succ[i] = i + 1;              // arbitrary valid index; won't be used
            P[i]    = -1e300;             // NEG → “take” never chosen
        }
    }

    // (2) DP: C[i][k] = best prize using k routes starting at index i
    const double NEG = -1e300;
    std::vector<std::vector<double>> C(n + 1, std::vector<double>(m_routes + 1, NEG));
    for (int k = 0; k <= m_routes; ++k) C[n][k] = 0.0;  // past-the-end = 0

    for (int i = n - 1; i >= 0; --i) {
        C[i][0] = 0.0;                                   // no routes left
        for (int k = 1; k <= m_routes; ++k) {
            double take = (lmax[i] >= 0) ? (P[i] + C[succ[i]][k-1]) : NEG;
            double skip = C[i+1][k];
            C[i][k] = std::max(take, skip);
        }
    }

    // (3) backtrack (take saturated vs skip)
    int i = 0, k = m_routes; double total_prize = C[0][m_routes];
    while (i < n && k > 0) {
        if (C[i][k] == C[i+1][k]) { ++i; continue; }     // skip
        int last = i + lmax[i];                          // take
        std::vector<int> r; r.reserve(last - i + 3);
        r.push_back(I.s);
        for (int p = i; p <= last; ++p) r.push_back(perm[p]);
        r.push_back(I.t);
        routes.push_back(std::move(r));
        i = succ[i]; --k;
    }
    return {routes, total_prize};
}


std::pair<std::vector<std::vector<int>>, double>
naive_split(const Instance& I, double B,
            const std::vector<int>& perm,
            int m_routes) {
    std::vector<std::vector<int>> routes;

    std::vector<int> path{ I.s };
    double route_cost = 0.0;
    int cur = I.s;

    for (size_t i = 0; i < perm.size(); ++i) {
        int v = perm[i];
        if (v == I.s || v == I.t) continue;

        // cost to go s…cur -> v, then be able to close v -> t
        double leg = I.costs[cur][v];
        double close_extra = leg + I.costs[v][I.t];

        if (route_cost + close_extra <= B + 1e-12) {
            // fits current route
            path.push_back(v);
            route_cost += leg;
            cur = v;
        } else {
            // close current route with t
            route_cost += I.costs[cur][I.t];
            path.push_back(I.t);
            routes.push_back(path);

            // start a new route from s
            path.assign(1, I.s);
            cur = I.s;
            route_cost = 0.0;

            // try to place v into the new route as first visit
            leg = I.costs[cur][v];
            close_extra = leg + I.costs[v][I.t];
            if (close_extra <= B + 1e-12) {
                path.push_back(v);
                route_cost += leg;
                cur = v;
            } 
        }
    }
    // close the last route
    route_cost += I.costs[cur][I.t];
    path.push_back(I.t);
    routes.push_back(path);

    // enforce m_routes: truncate or pad with empty {s,t} routes
    if ((int)routes.size() > m_routes) {
        routes = std::vector<std::vector<int>>(routes.begin(),
                                               routes.begin() + m_routes);
    } else if ((int)routes.size() < m_routes) {
        for (int k = (int)routes.size(); k < m_routes; ++k)
            routes.push_back({ I.s, I.t });
    }

    // total prize of produced routes (exclude s,t)
    double total_prize = 0.0;
    for (const auto& r : routes)
        for (size_t i = 1; i + 1 < r.size(); ++i)
            total_prize += I.prizes[r[i]];

    return { routes, total_prize };
}



std::pair<std::vector<std::vector<int>>, double>
naive_split(const std::vector<std::vector<double>>& costs,
            const std::vector<double>& prizes,
            int s, int t,
            double B,
            const std::vector<int>& perm,
            int m_routes) {
    std::vector<std::vector<int>> routes;

    std::vector<int> path{ s };
    double route_cost = 0.0;
    int cur = s;

    for (size_t i = 0; i < perm.size(); ++i) {
        int v = perm[i];
        if (v == s || v == t) continue;

        // cost to go s…cur -> v, then be able to close v -> t
        double leg = costs[cur][v];
        double close_extra = leg + costs[v][t];

        if (route_cost + close_extra <= B + 1e-12) {
            // fits current route
            path.push_back(v);
            route_cost += leg;
            cur = v;
        } else {
            // close current route with t
            route_cost += costs[cur][t];
            path.push_back(t);
            routes.push_back(path);

            // start a new route from s
            path.assign(1, s);
            cur = s;
            route_cost = 0.0;

            // try to place v into the new route as first visit
            leg = costs[cur][v];
            close_extra = leg + costs[v][t];
            if (close_extra <= B + 1e-12) {
                path.push_back(v);
                route_cost += leg;
                cur = v;
            } 
        }
    }
    // close the last route
    route_cost += costs[cur][t];
    path.push_back(t);
    routes.push_back(path);

    // enforce m_routes: truncate or pad with empty {s,t} routes
    if ((int)routes.size() > m_routes) {
        routes = std::vector<std::vector<int>>(routes.begin(),
                                               routes.begin() + m_routes);
    } else if ((int)routes.size() < m_routes) {
        for (int k = (int)routes.size(); k < m_routes; ++k)
            routes.push_back({ s, t });
    }

    // total prize of produced routes (exclude s,t)
    double total_prize = 0.0;
    for (const auto& r : routes)
        for (size_t i = 1; i + 1 < r.size(); ++i)
            total_prize += prizes[r[i]];

    return {routes, total_prize};
}


// std::pair<std::vector<std::vector<int>>, double>
// greedy_split (const std::vector<std::vector<double>>& costs,
//               const std::vector<double>& prizes,
//               const std::vector<int>& starts, int t,
//               double B,
//               const std::vector<int>& perm, // includes one start at front and the common sink t at back
//               int m_routes) {

//     const int S = (int)starts.size();
//     std::vector<std::vector<int>> routes;
//     routes.reserve(m_routes);

//     // Sanity: need at least [s, t] in perm
//     if (perm.size() < 2) {
//         // no path possible; pad empty {s_i, t} routes
//         std::vector<std::vector<int>> empty_routes;
//         for (int i = 0; i < m_routes; ++i) {
//             int si = (i < S ? i : 0);
//             empty_routes.push_back({ starts[si], t });
//         }
//         return { empty_routes, 0.0 };
//     }

//     // Extract s and t from perm; if t param disagrees, prefer perm.back()
//     int s0 = perm.front();
//     int t_node = perm.back();
//     if (t_node != t) t_node = t; // or keep perm.back(); choose one policy; here we honor function arg 't'

//     // Mark the start equal to s0 as used (if present in 'starts')
//     std::vector<bool> used(S, false);
//     for (int i = 0; i < S; ++i) if (starts[i] == s0) { used[i] = true; break; }

//     // Helper: pick nearest (unused if possible) start to customer v
//     auto pick_start_idx = [&](int v)->int {
//         int best = -1; double bestd = std::numeric_limits<double>::infinity();
//         for (int i = 0; i < S; ++i) if (!used[i]) {
//             double d = costs[starts[i]][v];
//             if (d < bestd) { bestd = d; best = i; }
//         }
//         if (best != -1) return best;
//         // all used: allow reuse; pick nearest overall
//         best = 0; bestd = std::numeric_limits<double>::infinity();
//         for (int i = 0; i < S; ++i) {
//             double d = costs[starts[i]][v];
//             if (d < bestd) { bestd = d; best = i; }
//         }
//         return best;
//     };

//     // Build first route starting at s0
//     std::vector<int> path;
//     path.reserve(2 + perm.size());
//     path.push_back(s0);

//     double route_cost = 0.0;
//     int cur = s0;

//     // Pack customers in order: perm[1]..perm[n-2] (skip front s0 and back t)
//     const std::size_t n = perm.size();
//     for (std::size_t i = 1; i + 1 < n; ++i) {
//         int v = perm[i];

//         // attempt to append v and still be able to close v -> t_node
//         double leg = costs[cur][v];
//         double close_extra = leg + costs[v][t_node];

//         if (route_cost + close_extra <= B + 1e-12) {
//             // fits current route
//             path.push_back(v);
//             route_cost += leg;
//             cur = v;
//         } else {
//             // close current route to t_node
//             route_cost += costs[cur][t_node];
//             path.push_back(t_node);
//             routes.push_back(path);

//             // start a new route for v: pick start nearest to v
//             int s_idx = pick_start_idx(v);
//             used[s_idx] = true;
//             int s_new = starts[s_idx];

//             path.clear();
//             path.push_back(s_new);
//             route_cost = 0.0;
//             cur = s_new;

//             // try to place v as first visit
//             double first_leg = costs[cur][v];
//             double first_close = first_leg + costs[v][t_node];
//             if (first_close <= B + 1e-12) {
//                 path.push_back(v);
//                 route_cost += first_leg;
//                 cur = v;
//             } else {
//                 // v cannot fit even alone with s_new…v…t_node under budget B:
//                 // treat as remainder (skip). If you want to record remainders, collect them in a vector here.
//                 continue;
//             }
//         }
//     }

//     // Close the last (possibly s-only) route to t_node
//     route_cost += costs[cur][t_node];
//     path.push_back(t_node);
//     routes.push_back(path);

//     // Enforce exactly m_routes: truncate or pad with empty {s_i, t_node} routes
//     if ((int)routes.size() > m_routes) {
//         routes.resize(m_routes);
//     } else if ((int)routes.size() < m_routes) {
//         int base = 0;
//         while ((int)routes.size() < m_routes) {
//             int si = (base < S ? base : 0);
//             routes.push_back({ starts[si], t_node });
//             ++base;
//         }
//     }

//     // Total prize (exclude depots and t_node)
//     double total_prize = 0.0;
//     for (const auto& r : routes) {
//         for (std::size_t i = 1; i + 1 < r.size(); ++i)
//             total_prize += prizes[r[i]];
//     }

//     return { routes, total_prize };
// }


// std::pair<std::vector<std::vector<int>>, double>
// greedy_split(const std::vector<std::vector<double>>& costs,
//              const std::vector<double>& prizes,
//              const std::vector<int>& starts, int t,
//              double B,
//              const std::vector<int>& perm,  // customers in visit order (no s, no t)
//              int m_routes) {

//     const int S = (int)starts.size();
//     std::vector<bool> used(S, false);
//     std::vector<std::vector<int>> routes;

//     // Helper: pick an unused start nearest to node v (falls back to nearest overall if all used)
//     auto pick_start_idx = [&](int v)->int {
//         int best = -1; double bestd = std::numeric_limits<double>::infinity();
//         // prefer unused
//         for (int i = 0; i < S; ++i) if (!used[i]) {
//             double d = costs[starts[i]][v];
//             if (d < bestd) { bestd = d; best = i; }
//         }
//         if (best != -1) return best;
//         // all used — allow reuse, pick nearest overall
//         bestd = std::numeric_limits<double>::infinity(); best = 0;
//         for (int i = 0; i < S; ++i) {
//             double d = costs[starts[i]][v];
//             if (d < bestd) { bestd = d; best = i; }
//         }
//         return best;
//     };

//     // Early out: no customers
//     if (perm.empty()) {
//         // create empty {s,t} routes
//         for (int k = 0; k < m_routes; ++k) {
//             int si = k < S ? k : 0;
//             routes.push_back({ starts[si], t });
//         }
//         return { routes, 0.0 };
//     }

//     // Initialize first route: choose start depot closest to first customer
//     int first = perm.front();
//     int s_idx = pick_start_idx(first);
//     used[s_idx] = true;
//     int s = starts[s_idx];

//     std::vector<int> path;
//     path.reserve(2 + perm.size());
//     path.push_back(s);

//     double route_cost = 0.0;
//     // go s -> first
//     route_cost += costs[s][first];
//     path.push_back(first);
//     int cur = first;

//     // Greedy pack along perm
//     for (std::size_t i = 1; i < perm.size(); ++i) {
//         int v = perm[i];

//         // try to extend current path with v and still be able to close to t
//         double leg = costs[cur][v];
//         double close_extra = leg + costs[v][t];

//         if (route_cost + close_extra <= B + 1e-12) {
//             // fits current route
//             path.push_back(v);
//             route_cost += leg;
//             cur = v;
//         } else {
//             // close current route to t
//             route_cost += costs[cur][t];
//             path.push_back(t);
//             routes.push_back(path);

//             // start a new route (if allowed)
//             // pick start depot nearest to next customer v
//             int next_idx = pick_start_idx(v);
//             used[next_idx] = true;
//             int s2 = starts[next_idx];

//             path.clear();
//             path.push_back(s2);
//             route_cost = 0.0;
//             // place v as first visit if it fits with closing to t
//             double first_leg = costs[s2][v];
//             double first_close_extra = first_leg + costs[v][t];
//             if (first_close_extra <= B + 1e-12) {
//                 path.push_back(v);
//                 route_cost += first_leg;
//                 cur = v;
//             } else {
//                 // cannot place v in a fresh route under budget — skip or handle as remainder
//                 // Here we simply skip v (remainder). You can store it somewhere if needed.
//                 continue;
//             }
//         }
//     }

//     // close the last route
//     route_cost += costs[cur][t];
//     path.push_back(t);
//     routes.push_back(path);

//     // Enforce exactly m_routes: truncate or pad with empty {s,t} routes
//     if ((int)routes.size() > m_routes) {
//         routes.resize(m_routes);
//     } else if ((int)routes.size() < m_routes) {
//         int base = 0;
//         while ((int)routes.size() < m_routes) {
//             int si = (base < S ? base : 0);
//             routes.push_back({ starts[si], t });
//             ++base;
//         }
//     }

//     // total prize (exclude depots and t)
//     double total_prize = 0.0;
//     for (const auto& r : routes) {
//         for (std::size_t i = 1; i + 1 < r.size(); ++i) total_prize += prizes[r[i]];
//     }

//     return { routes, total_prize };
// }

std::pair<std::vector<std::vector<int>>, double>
greedy_split(const std::vector<std::vector<double>>& costs,
             const std::vector<double>& prizes,
             const std::vector<int>& starts, int t,
             double B,
             const std::vector<int>& perm, // includes one start at front and common sink t at back
             int m_routes) {

    const int S = (int)starts.size();
    std::vector<std::vector<int>> routes;
    routes.reserve(m_routes);

    if (perm.size() < 2 || S == 0 || m_routes <= 0) {
        return { {}, 0.0 };
    }

    const int t_node = t; 

    // Build first route at starts[0]
    int r = 0; 
    std::vector<int> path;
    path.reserve(2 + perm.size());
    path.push_back(starts[r]);

    double route_cost = 0.0;
    int cur = starts[r];

    const std::size_t n = perm.size();

    // Pack customers in the order of perm[1]..perm[n-2] (skip front and back)
    for (std::size_t i = 1; i + 1 < n; ++i) {
        int v = perm[i];

        // Try appending v and still closing v -> t within budget B
        double leg = costs[cur][v];
        double close_extra = leg + costs[v][t_node];

        if (route_cost + close_extra <= B + 1e-12) {
            // fits current route
            path.push_back(v);
            route_cost += leg;
            cur = v;
        } else {
            // close current route to t
            route_cost += costs[cur][t_node];
            path.push_back(t_node);
            routes.push_back(path);

            // move to next route (keep start order unchanged)
            ++r;
            if (r >= m_routes) {
                // no more routes allowed; skip remaining customers (they become remainder)
                // If you prefer, break here and pad/truncate below.
                break;
            }

            // start new route at starts[r]
            path.clear();
            path.push_back(starts[r]);
            route_cost = 0.0;
            cur = starts[r];

            // try to place v as the first visit in the new route
            double first_leg = costs[cur][v];
            double first_close = first_leg + costs[v][t_node];
            if (first_close <= B + 1e-12) {
                path.push_back(v);
                route_cost += first_leg;
                cur = v;
            } else {
                // v cannot fit even alone with starts[r]…v…t under B -> skip (remainder)
                // (Optionally collect remainders in a vector if needed)
                continue;
            }
        }
    }

    // Close the current (possibly start-only) route if we are still within m_routes
    if ((int)routes.size() < m_routes) {
        route_cost += costs[cur][t_node];
        path.push_back(t_node);
        routes.push_back(path);
        ++r;
    }

    // Pad with empty {starts[r], t} routes to reach exactly m_routes (sequence preserved)
    while (r < m_routes) {
        int s_fixed = starts[r];
        routes.push_back({ s_fixed, t_node });
        ++r;
    }

    // Total prize (exclude depots and t)
    double total_prize = 0.0;
    for (const auto& rt : routes) {
        for (std::size_t i = 1; i + 1 < rt.size(); ++i)
            total_prize += prizes[rt[i]];
    }

    return { routes, total_prize };
}

static double optimistic_delta_cost(const std::vector<std::vector<double>>& costs,
                                    const std::vector<int>& starts,
                                    int t, int v) {
    double best = std::numeric_limits<double>::infinity();
    for (int s : starts) {
        double d = costs[s][v] + costs[v][t];
        if (d < best) best = d;
    }
    return best;
}

std::pair<std::vector<std::vector<int>>, double>
greedy_split_prune_until_feasible(const std::vector<std::vector<double>>& costs,
                                  const std::vector<double>& prizes,
                                  const std::vector<int>& starts, int t,
                                  double B,
                                  const std::vector<int>& perm_in, // includes one start at front and common t at back
                                  int m_routes) {
    // Guard: need at least [s0, t]
    if (perm_in.size() < 2) {
        std::vector<std::vector<int>> empty_routes;
        for (int i = 0; i < m_routes; ++i) {
            int si = (i < (int)starts.size() ? i : 0);
            empty_routes.push_back({ starts[si], t });
        }
        return { empty_routes, 0.0 };
    }

    const int n_all = (int)perm_in.size();
    const int s0 = perm_in.front();
    const int t_node = t;

    // Initial active set = all customers in perm_in (skip front s0 and back t)
    std::vector<int> active;
    active.reserve(n_all);
    for (std::size_t i = 1; i + 1 < perm_in.size(); ++i) active.push_back(perm_in[i]);

    // Helper to build working permutation [s0, active in original order, t]
    auto build_perm_work = [&](const std::vector<int>& active_list)->std::vector<int> {
        std::unordered_set<int> keep(active_list.begin(), active_list.end());
        std::vector<int> out; out.reserve(2 + active_list.size());
        out.push_back(s0);
        for (std::size_t i = 1; i + 1 < perm_in.size(); ++i) {
            int v = perm_in[i];
            if (keep.find(v) != keep.end()) out.push_back(v);
        }
        out.push_back(t_node);
        return out;
    };

    // Track best found (highest prize) in case we stop early
    std::pair<std::vector<std::vector<int>>, double> best_sol;
    best_sol.second = -1.0;

    // Limit iterations to number of customers (each step can drop at most one)
    const int max_iters = (int)active.size() + 1;
    for (int it = 0; it < max_iters; ++it) {
        // Build current working permutation and split
        std::vector<int> perm_work = build_perm_work(active);
        auto sol = greedy_split(costs, prizes, starts, t_node, B, perm_work, m_routes);
        std::vector<std::vector<int>>& routes = sol.first;
        repairMultiPaths(costs, prizes, routes, B, starts,t);
        double total_prize = sol.second;

        // Record best
        if (total_prize > best_sol.second) best_sol = sol;

        // Gather served customers from routes (exclude first and last of each route)
        std::unordered_set<int> served;
        for (const auto& r : routes) {
            for (std::size_t i = 1; i + 1 < r.size(); ++i) served.insert(r[i]);
        }

        // Find unserved among active (in perm order)
        std::vector<int> unserved;
        unserved.reserve(active.size());
        for (int v : active) if (served.find(v) == served.end()) unserved.push_back(v);

        if (unserved.empty()) {
            // Feasible: everyone in 'active' is routed
            return sol; // done
        }

        // Pick one to drop: smallest prize / optimistic_delta_cost
        double worst_score = std::numeric_limits<double>::infinity();
        int worst_idx = -1;
        for (std::size_t i = 0; i < unserved.size(); ++i) {
            int v = unserved[i];
            double delta = optimistic_delta_cost(costs, starts, t_node, v);
            double eps = 1e-12;
            double score = prizes[v] / (delta + eps);
            if (score < worst_score) {
                worst_score = score;
                worst_idx = (int)i;
            }
        }
        if (worst_idx < 0) break; // safety

        // Remove that node from 'active'
        int victim = unserved[worst_idx];
        auto it_rm = std::find(active.begin(), active.end(), victim);
        if (it_rm != active.end()) active.erase(it_rm);
        else break; // shouldn't happen
    }

    // If we exit, return the best we saw
    return best_sol.second >= 0.0 ? best_sol
                                  : std::pair<std::vector<std::vector<int>>, double>({/*routes*/}, 0.0);
}

