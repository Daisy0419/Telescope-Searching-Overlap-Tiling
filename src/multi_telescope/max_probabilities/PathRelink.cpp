#include "LocalSearch.h"
#include "PathOperations.h"

#include <vector>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <limits>
#include <unordered_set>
#include <cmath>
#include <cassert>
#include <random> 
#include <iostream>
#include <atomic>


std::vector<std::vector<int>>
construct_paths(const std::vector<std::vector<double>>& costs,
                const std::vector<double>& prizes,
                int start, int end, int nPath,
                double budgetT, double greediness) {

    const double EPS = 1e-12;
    const int n = (int)costs.size();

    if (costs[start][end] > budgetT + EPS) {
        std::vector<std::vector<int>> triv(nPath, std::vector<int>{start, end});
        return triv;
    }

    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double alpha = (greediness < 0.0 ? U01(rng) : std::clamp(greediness, 0.0, 1.0));

    // init routes as [s, t]
    std::vector<std::vector<int>> paths(nPath, std::vector<int>{start, end});

    // per-route costs: c(s,t)
    std::vector<double> pathCosts(nPath, costs[start][end]);

    std::vector<char> used(n, 0);
    used[start] = 1; used[end] = 1;

    struct Cand { int v, r, pos; double dt, h; };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue;           
            for (int r = 0; r < nPath; ++r) {
                auto [pos, dt] = cheapestInsertionOnePath(costs, paths[r], v);
                if (pos == -1) continue;
                if (pathCosts[r] + dt <= budgetT + EPS) {
                    double denom = std::max(dt, EPS);
                    double h = prizes[v] / denom; 
                    C.push_back({v, r, pos, dt, h});
                }
            }
        }
    };

    int steps = 0;
    while (steps++ < 10000) {
        std::vector<Cand> C; buildCandidates(C);
        if (C.empty()) break;

        // Value-based RCL: threshold = minH + alpha * (maxH - minH)
        double minH = C.front().h, maxH = C.front().h;
        for (const auto& c : C) { minH = std::min(minH, c.h); maxH = std::max(maxH, c.h); }
        const double thr = minH + alpha * (maxH - minH);

        std::vector<int> RCL; RCL.reserve(C.size());
        for (int i = 0; i < (int)C.size(); ++i)
            if (C[i].h + EPS >= thr) RCL.push_back(i);
        if (RCL.empty()) break;

        std::uniform_int_distribution<int> pick(0, (int)RCL.size() - 1);
        const Cand& mv = C[RCL[pick(rng)]];

        // Apply insertion 
        auto& route = paths[mv.r];
        route.insert(route.begin() + mv.pos, mv.v);
        pathCosts[mv.r] += mv.dt;
        used[mv.v] = 1;

    }

    return paths;
}


// Multi-start constructor: start_indices[r] is the start for route r, shared end.
std::vector<std::vector<int>>
construct_paths_multiS(const std::vector<std::vector<double>>& costs,
                       const std::vector<double>& prizes,
                       const std::vector<int>& start_indices,
                       int end,
                       int nPath,               
                       double budgetT,
                       double greediness) {
    const double EPS = 1e-12;
    const int n = (int)costs.size();
    if ((int)prizes.size() != n) throw std::runtime_error("prizes size mismatch");
    if ((int)start_indices.size() != nPath) throw std::runtime_error("start_indices.size() != nPath");
    if (end < 0 || end >= n) throw std::runtime_error("bad end index");

    {
        std::vector<char> seen(n, 0);
        for (int s : start_indices) {
            if (s < 0 || s >= n) throw std::runtime_error("bad start index");
            if (s == end) throw std::runtime_error("start cannot equal end");
            if (seen[s]) throw std::runtime_error("duplicate start in start_indices");
            seen[s] = 1;
        }
    }

    // RNG and alpha (sample once per construct run)
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double alpha = (greediness < 0.0 ? U01(rng) : std::clamp(greediness, 0.0, 1.0));

    // Initialize routes as [s_r, end] and base costs = c(s_r, end)
    std::vector<std::vector<int>> paths(nPath);
    std::vector<double> pathCosts(nPath, 0.0);
    for (int r = 0; r < nPath; ++r) {
        int s_r = start_indices[r];
        paths[r] = { s_r, end };
        pathCosts[r] = costs[s_r][end]; 
    }

    std::vector<char> used(n, 0);
    for (int s_r : start_indices) used[s_r] = 1;
    used[end] = 1;

    struct Cand { int v, r, pos; double dt, h; };

    auto cheapestInsertionOnePath = [&](const std::vector<int>& route, int v) -> std::pair<int,double> {
        int best_pos = -1; double best_dt = std::numeric_limits<double>::infinity();
        if (route.size() < 2) return { -1, best_dt };
        for (std::size_t j = 0; j + 1 < route.size(); ++j) {
            int a = route[j], b = route[j+1];
            double dt = costs[a][v] + costs[v][b] - costs[a][b];
            if (std::isnan(dt)) continue;
            if (dt < best_dt) { best_dt = dt; best_pos = (int)j + 1; }
        }
        return { best_pos, best_dt };
    };

    auto buildCandidates = [&](std::vector<Cand>& C){
        C.clear();
        for (int v = 0; v < n; ++v) {
            if (used[v]) continue; // only unvisited
            for (int r = 0; r < nPath; ++r) {
                auto [pos, dt] = cheapestInsertionOnePath(paths[r], v);
                if (pos == -1) continue;
                if (pathCosts[r] + dt <= budgetT + EPS) {
                    double denom = std::max(dt, EPS);
                    double h = prizes[v] / denom; // value-per-extra-cost
                    C.push_back({v, r, pos, dt, h});
                }
            }
        }
    };

    int steps = 0;
    while (steps++ < 10000) {
        std::vector<Cand> C; C.reserve(256);
        buildCandidates(C);
        if (C.empty()) break;

        // Value-based RCL: threshold = minH + alpha * (maxH - minH)
        double minH = C.front().h, maxH = C.front().h;
        for (const auto& c : C) { minH = std::min(minH, c.h); maxH = std::max(maxH, c.h); }
        const double thr = minH + alpha * (maxH - minH);

        std::vector<int> RCL; RCL.reserve(C.size());
        for (int i = 0; i < (int)C.size(); ++i)
            if (C[i].h + EPS >= thr) RCL.push_back(i);
        if (RCL.empty()) break;

        std::uniform_int_distribution<int> pick(0, (int)RCL.size() - 1);
        const Cand& mv = C[RCL[pick(rng)]];

        // Apply insertion on route mv.r
        auto& route = paths[mv.r];
        route.insert(route.begin() + mv.pos, mv.v);
        pathCosts[mv.r] += mv.dt;
        used[mv.v] = 1;
    }

    return paths;
}



// //randomized construct with greediness
// std::vector<std::vector<int>> construct_paths(const std::vector<std::vector<double>>& costs,  
//                                             const std::vector<double>& prizes,
//                                             int start, int end, int nPath,
//                                             double budgetT, double greediness) {

//     const double EPS = 1e-12;
//     int n = int(costs.size());

//     // init RNG
//     std::mt19937_64 rng(std::random_device{}());
//     std::uniform_real_distribution<double> dist(0, 1);

//     // init routes as [s_r, t_r]
//     std::vector<std::vector<int>> paths(nPath);
//     for (int r = 0; r < nPath; ++r) paths[r] = {start, end};

//     // per-route cost
//     std::vector<double> pathCosts(nPath, 0.0);
//     for (int r = 0; r < nPath; ++r) pathCosts[r] = costs[start][end];

//     // global "used" (internal nodes only); mark all endpoints as used to prevent re-use
//     std::vector<bool> visited(n, false);
//     visited[start] = true;
//     visited[end] = true;

//     struct Candidate {
//         int node; double score; int path_id; int pos; double dt; 
//     };

//     auto buildCandidates = [&](std::vector<Candidate>& C){
//         C.clear();
//         for (int node = 0; node < n; ++node) {
//             if (visited[node]) continue; 
//             for (int p = 0; p < nPath; ++p) {
//                 auto [pos, dt] = cheapestInsertionOnePath(costs, paths[p], node);
//                 if (pos == -1) continue;
//                 if (pathCosts[p] + dt <= budgetT + EPS) {
//                     double denom = std::max(dt, EPS);
//                     double score = prizes[node] / denom;     // maximize S/Δt
//                     C.push_back({node, score, p, pos, dt});
//                 }
//             }
//         }
//     };

//     int steps = 0;
//     while (steps++ < 1000) {
//         std::vector<Candidate> C;
//         buildCandidates(C);
//         if (C.empty()) break;

//         // value-based RCL: threshold = minH + alpha*(maxH - minH)
//         double minH = C.front().score, maxH = C.front().score;
//         for (const auto& c : C) { 
//             minH = std::min(minH, c.score); 
//             maxH = std::max(maxH, c.score); 
//         }

//         greediness = dist(rng); 
//         double thr = minH + greediness * (maxH - minH);

//         std::vector<int> restrictedCandList; 
//         restrictedCandList.reserve(C.size());
//         for (int i = 0; i < (int)C.size(); ++i) 
//             if (C[i].score + EPS >= thr) restrictedCandList.push_back(i);
//         if (restrictedCandList.empty()) break;

//         std::uniform_int_distribution<int> pick(0, (int)restrictedCandList.size()-1);
//         const Candidate& mv = C[restrictedCandList[pick(rng)]];

//         // apply
//         paths[mv.path_id].insert(paths[mv.path_id].begin() + mv.pos, mv.node);
//         pathCosts[mv.path_id] += mv.dt;
//         visited[mv.node] = 1;
//     }

//     return paths;
// }


// Internal-node mask across all routes; 
// size nNodes; 1 if internal node is visited
// std::vector<bool> getVisited(const std::vector<std::vector<int>>& paths, int nNodes) {
//     std::vector<bool> visited(nNodes, false);
//     for (const auto& r : paths)
//         for (int i = 1; i+1 < (int)r.size(); ++i)
//             visited[r[i]] = true;
//     return visited;
// }

// Similarity on internal-node sets:  2|A∩B| / (|A| + |B|)
double solutionSimilarity(const std::vector<std::vector<int>>& A,
                        const std::vector<std::vector<int>>& B,
                        int nNodes) {

    auto a = getVisited(A, nNodes);
    auto b = getVisited(B, nNodes);
    int inter = 0, na = 0, nb = 0;
    for (int v = 0; v < nNodes; ++v) {
        na += a[v] ? 1 : 0;
        nb += b[v] ? 1 : 0;
        inter += (a[v] && b[v]) ? 1 : 0;
    }
    if (na + nb == 0) return 1.0;
    return 2.0 * inter / double(na + nb);
}

// // Best insertion on a specific s–t route: insert at index 'pos' (between pos-1 and pos), pos ∈ [1, n-1]
// std::pair<int,double> bestInsertionOnRoute(const std::vector<std::vector<double>>& costs,
//                                                          const std::vector<int>& route,
//                                                          int v) {
//     const int n = (int)route.size();
//     if (n < 2) return {-1, std::numeric_limits<double>::infinity()};
//     int bestPos = -1; double bestDt = std::numeric_limits<double>::infinity();
//     for (int pos = 1; pos <= n-1; ++pos) {
//         int a = route[pos-1], b = route[pos];
//         double dt = costs[a][v] + costs[v][b] - costs[a][b];
//         if (dt < bestDt) { bestDt = dt; bestPos = pos; }
//     }
//     return {bestPos, bestDt};
// }

// // Cheapest feasible insertion across routes
// std::tuple<int,int,double> cheapestInsertionAll(const std::vector<std::vector<double>>& costs,
//                                             const std::vector<std::vector<int>>& paths,
//                                             const std::vector<double>& pathCosts,
//                                             double budgetT,
//                                             int v) {
//     int bestRoute = -1, bestPos = -1; double bestDt = std::numeric_limits<double>::infinity();
//     for (int r = 0; r < (int)paths.size(); ++r) {
//         auto [pos, dt] = cheapestInsertionOnePath(costs, paths[r], v);
//         if (pos == -1) continue;
//         if (pathCosts[r] + dt <= budgetT + 1e-12 && dt < bestDt) {
//             bestDt = dt; bestRoute = r; bestPos = pos;
//         }
//     }
//     if (bestRoute == -1) return {-1,-1,0.0};
//     return {bestRoute, bestPos, bestDt};
// }

// // 2-Opt per route
// bool twoOptAllRoutes(const std::vector<std::vector<double>>& costs,
//                     std::vector<std::vector<int>>& paths,
//                     std::vector<double>& pathCosts,
//                     double budgetT) {

//     bool any = false;
//     for (int r = 0; r < (int)paths.size(); ++r) {
//         auto& route = paths[r];
//         int n = (int)route.size();
//         if (n <= 4) continue;

//         bool improved = true;
//         while (improved) {
//             improved = false;
//             double bestDelta = -1e-12; int bi=-1, bk=-1;
//             const double oldT = pathCosts[r];

//             for (int i = 1; i <= n-3; ++i) {
//                 for (int k = i+1; k <= n-2; ++k) {
//                     std::reverse(route.begin()+i, route.begin()+k+1);
//                     double newT = pathCost(costs, route);
//                     std::reverse(route.begin()+i, route.begin()+k+1);
//                     double d = newT - oldT;
//                     if (d < bestDelta && newT <= budgetT + 1e-12) {
//                         bestDelta = d; bi=i; bk=k;
//                     }
//                 }
//             }
//             if (bi != -1) {
//                 std::reverse(route.begin()+bi, route.begin()+bk+1);
//                 pathCosts[r] = pathCost(costs, route);
//                 improved = true; any = true;
//             }
//         }
//     }
//     return any;
// }


void localSearchLoop(const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<int>>& paths,
                    std::vector<double>& pathCosts,
                    const std::vector<double>& prizes,
                    double budgetT) {

    while (true) {
        bool changed = false;
        changed |= twoOptAllRoutes(costs, paths, pathCosts, budgetT);
        // std::cout << "twoOptAllRoutes: \n";
        // verify_routes(costs, prizes, paths[0].front(), paths[0].back(), paths, budgetT);
        changed |= swapNodesBetweenPaths(costs, paths, pathCosts,budgetT);
        // std::cout << "swapNodesBetweenPaths: \n";
        // verify_routes(costs, prizes, paths[0].front(), paths[0].back(), paths, budgetT);
        changed |= insertNodesToPaths(costs, paths, pathCosts, prizes, budgetT);
        // std::cout << "insertNodesToPaths: \n";
        // verify_routes(costs, prizes, paths[0].front(), paths[0].back(), paths, budgetT);
        changed |= replaceNodesInPaths(costs, paths, pathCosts, prizes, budgetT);
        // std::cout << "replaceNodesInPaths: \n";
        // verify_routes(costs, prizes, paths[0].front(), paths[0].back(), paths, budgetT);
        if (!changed) break;
    }
}


bool repairPathWorstEfficiency(const std::vector<std::vector<double>>& costs,
                                std::vector<int>& path,
                                double& pathCost,
                                const std::vector<double>& prizes,
                                double budgetT) {
    bool removedAny = false;
    const double EPS = 1e-12;

    while (pathCost > budgetT + 1e-12) {
        const int n = (int)path.size();
        int bestIdx = -1;
        double bestScore = -1.0; // larger = worse (save/prize)

        for (int i = 1; i+1 < n; ++i) {                 // internal only
            int v = path[i];
            double prize = std::max(prizes[v], EPS);
            int a = path[i-1], b = path[i+1];
            double save = costs[a][v] + costs[v][b] - costs[a][b];
            double eff  = save / prize;                 // "worst" efficiency
            if (eff > bestScore) { bestScore = eff; bestIdx = i; }
        }
        if (bestIdx == -1) break;

        int a = path[bestIdx-1], v = path[bestIdx], b = path[bestIdx+1];
        double save = costs[a][v] + costs[v][b] - costs[a][b];
        path.erase(path.begin() + bestIdx);
        pathCost -= save;
        removedAny = true;

        if ((int)path.size() <= 1) break;
    }
    return removedAny;
}

bool repairAllRoutes(const std::vector<std::vector<double>>& costs,
                    std::vector<std::vector<int>>& paths,
                    std::vector<double>& pathCosts,
                    const std::vector<double>& prizes,
                    double budgetT) {
    bool any = false;
    for (int r = 0; r < (int)paths.size(); ++r)
        if (pathCosts[r] > budgetT + 1e-12)
            any |= repairPathWorstEfficiency(costs, paths[r], pathCosts[r], prizes, budgetT);
    return any;
}


struct Solution {
    std::vector<std::vector<int>> paths;
    std::vector<double> pathCosts;
    double score = 0.0;
};

// linearized order of internal nodes in guide (any order is fine; this uses route order)
static inline std::vector<int> guidingOrder(const std::vector<std::vector<int>>& guide) {
    std::vector<int> seq;
    for (const auto& r : guide)
        for (int i = 1; i+1 < (int)r.size(); ++i)
            seq.push_back(r[i]);
    return seq;
}


// Check if every route is infeasible
inline bool allRoutesInfeasible(const std::vector<double>& pc, double T) {
    for (double t : pc) if (t <= T + 1e-12) return false;
    return true;
}

Solution pathRelink(const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes,
                    const std::vector<std::vector<int>>& startPaths,
                    const std::vector<std::vector<int>>& guidePaths,
                    double budgetT) {
    const int n = (int)prizes.size();

    Solution cur{startPaths,
                 multiPathsCost(costs, startPaths),
                 multiPathsSumPrize(prizes, startPaths)};
    Solution best = cur;

    // toAdd = internal nodes in guide but not currently present (internal) in start
    auto visited = getVisited(cur.paths, n);
    std::vector<int> seq = guidingOrder(guidePaths);
    std::vector<int> toAdd;
    toAdd.reserve(seq.size());
    for (int v : seq) if (!visited[v]) toAdd.push_back(v);

    while (!toAdd.empty()) {
        bool insertedThisRound = false;

        // Insert as long as at least one route is still feasible; allow the chosen route to become infeasible
        while (!toAdd.empty() && !allRoutesInfeasible(cur.pathCosts, budgetT)) {
            int v = toAdd.front();
            toAdd.erase(toAdd.begin());

            int bestR = -1, bestPos = -1; double bestDt = std::numeric_limits<double>::infinity();

            // only consider routes that are FEASIBLE BEFORE insertion (paper’s rule)
            for (int r = 0; r < (int)cur.paths.size(); ++r) {
                if (cur.pathCosts[r] > budgetT + 1e-12) continue;
                auto [pos, dt] = cheapestInsertionOnePath(costs, cur.paths[r], v);
                if (pos == -1) continue;
                if (dt < bestDt) { bestDt = dt; bestPos = pos; bestR = r; }
            }

            if (bestR == -1) break; // all routes already infeasible → go repair

            // perform insertion (even if it makes this route infeasible)
            cur.paths[bestR].insert(cur.paths[bestR].begin() + bestPos, v);
            cur.pathCosts[bestR] += bestDt;        // <-- keep times in sync
            visited[v] = 1;
            insertedThisRound = true;
        }

        // Repair infeasibility, then full LS; keep best along the path
        repairMultiPaths(costs,prizes, cur.paths, budgetT, cur.paths[0].front(), cur.paths[0].back());
        // repairAllRoutes(costs, cur.paths, cur.pathCosts, prizes, budgetT);
        localSearchLoop(costs, cur.paths, cur.pathCosts, prizes, budgetT);

        cur.score = multiPathsSumPrize(prizes, cur.paths);
        if (cur.score > best.score) best = cur;

        if (!insertedThisRound) break; // nothing more could be added; stop
    }

    return best;
}

//Elite Pool 

// struct Elite {
//     std::vector<std::vector<int>> paths;
//     std::vector<double> pathCosts;
//     double score = 0.0;
//     std::atomic<int> age{0};
// };
// #include <vector>
// #include <atomic>
// #include <utility>

struct Elite {
    std::vector<std::vector<int>> paths;
    std::vector<double>            pathCosts;
    double                         score = 0.0;
    std::atomic<int>               age{0};   // atomic age

    Elite() = default;

    Elite(std::vector<std::vector<int>> p,
          std::vector<double> c,
          double s,
          int a = 0)
        : paths(std::move(p)),
          pathCosts(std::move(c)),
          score(s),
          age(a) {}

    // Copy
    Elite(const Elite& other)
        : paths(other.paths),
          pathCosts(other.pathCosts),
          score(other.score),
          age(other.age.load(std::memory_order_relaxed)) {}

    Elite& operator=(const Elite& other) {
        if (this != &other) {
            paths = other.paths;
            pathCosts = other.pathCosts;
            score = other.score;
            age.store(other.age.load(std::memory_order_relaxed),
                      std::memory_order_relaxed);
        }
        return *this;
    }

    // Move
    Elite(Elite&& other) noexcept
        : paths(std::move(other.paths)),
          pathCosts(std::move(other.pathCosts)),
          score(other.score),
          age(other.age.load(std::memory_order_relaxed)) {}

    Elite& operator=(Elite&& other) noexcept {
        if (this != &other) {
            paths = std::move(other.paths);
            pathCosts = std::move(other.pathCosts);
            score = other.score;
            age.store(other.age.load(std::memory_order_relaxed),
                      std::memory_order_relaxed);
        }
        return *this;
    }
};


class ElitePool {
public:
    ElitePool(int pool_size)
        : K(pool_size), no_improve_iters(0),
          best_score(-std::numeric_limits<double>::infinity()) {
        pool.reserve(K);
    }

    // Single-shot update per iteration:
    //  - 'best_from_PR' = best solution found across ALL links this iteration
    void updateElitePool(const Solution& best_from_PR) {
        const int threshold = std::max(10, no_improve_iters / 10);
        pool.erase(std::remove_if(pool.begin(), pool.end(),
                     [&](const Elite& e){ return e.age >= threshold; }),
                   pool.end());

        //Consider inserting/replacing with the PR-best candidate
        Elite cand{ best_from_PR.paths, best_from_PR.pathCosts, best_from_PR.score, 0 };

        if ((int)pool.size() < K) {
            pool.push_back(std::move(cand));
        } else {
            auto worst_it = std::min_element(pool.begin(), pool.end(),
                             [](const Elite& a, const Elite& b){ return a.score < b.score; });
            if (worst_it != pool.end() && cand.score > worst_it->score) {
                *worst_it = std::move(cand);
            }
        }

        // 4)track global improvement here
        if (best_from_PR.score > best_score + 1e-12) {
            best_score = best_from_PR.score;
            no_improve_iters = 0;
        } else {
            ++no_improve_iters;
        }
    }

    int size() const { return (int)pool.size(); }
    std::vector<Elite> pool;
    int no_improve_iters = 0;
    int K;
    double best_score;
};


// linkToElites mutates 'pool' by aging used elites (side effect)
Solution linkToElites(const std::vector<std::vector<double>>& costs,
                      const std::vector<double>& prizes,
                      const Solution& current,
                      std::vector<Elite>& pool,   // non-const: ages inside
                      double budgetT,
                      double skipSimAbove = 0.99) {
    Solution best = current;
    const int n = (int)prizes.size();

    for (auto& E : pool) {
        double sim = solutionSimilarity(current.paths, E.paths, n);
        // std::cout << "solutionSimilarity: " << sim << std::endl;
        if (sim >= skipSimAbove) continue;

        // Mark “used once” for this iteration
        // E.age ++;
        E.age.fetch_add(1); 

        // PR both directions; keep the best found on either path
        Solution s1 = pathRelink(costs, prizes, current.paths, E.paths, budgetT);
        if (s1.score > best.score) best = s1;

        Solution s2 = pathRelink(costs, prizes, E.paths, current.paths, budgetT);
        if (s2.score > best.score) best = s2;
    }
    return best;
}

// Calls linkToElites (ages used elites), then updates the pool (purge/add/replace)
// Returns the best solution after PR (or the input 'current' if PR didn’t improve it).
inline Solution link_and_update(const std::vector<std::vector<double>>& costs,
                                const std::vector<double>& prizes,
                                const Solution& current,
                                ElitePool& elites,                 // owns 'pool' and counters
                                double budgetT,
                                double simSkip = 0.95) {
    // 1) Run PR against all eligible elites (linkToElites will age used elites)
    Solution best_from_PR = linkToElites(costs, prizes, current,
                                         elites.pool /* side-effect: ages used ones */,
                                         budgetT, simSkip);

    // 2) Update the elite pool per paper (purge by age, add/replace worst-by-score,
    //    and update no_improve_iters / best_score inside the pool)
    elites.updateElitePool(best_from_PR);

    // 3) Keep the better of current and PR-best as the iteration’s result
    return (best_from_PR.score > current.score ? best_from_PR : current);
}


std::pair<std::vector<std::vector<int>>, double> PROptimization(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               int start, int end) {
    // --- Tunables ---
    const int   POOL_SIZE        = 5;     // paper commonly uses 5
    const int   MAX_ITERS        = 300;  // hard cap
    const int   STOP_NO_IMPROVE  = 30;   // paper's "slow PR" ~300; use 10 for "fast"
    const double SIM_SKIP        = 0.95;  // skip PR if similarity >= this

    ElitePool elites(POOL_SIZE);

    // RNG for GRASP greediness alpha
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // Track global best
    Solution bestGlobal{ {}, {}, -std::numeric_limits<double>::infinity() };
    int no_improve_streak = 0;

    for (int it = 0; it < MAX_ITERS; ++it) {
        // ---- Construct (randomized GRASP) ----
        const double greediness = U01(rng); // alpha ~ U[0,1)
        auto construct = construct_paths(costs, prizes, start, end, nPath, budgetT, greediness);

        // ---- Local Search to a local optimum ----
        auto pathCosts = multiPathsCost(costs, construct);
        auto pathPrize = multiPathsSumPrize(prizes, construct);
        std::cout << "iteration: " << it << std::endl;
        std::cout << "construct prize: " << pathPrize << std::endl;
        std::cout << "global best prize: " << bestGlobal.score << std::endl;


        // return {construct, pathPrize};

        localSearchLoop(costs, construct, pathCosts, prizes, budgetT);
        pathCosts = multiPathsCost(costs, construct); // re-read after LS

        Solution current{ construct, pathCosts, multiPathsSumPrize(prizes, construct) };

        // ---- Path Relinking (both directions vs all eligible elites; ages used) ----
        Solution afterPR = linkToElites(costs, prizes, current, elites.pool, budgetT, SIM_SKIP);
        repairMultiPaths(costs,prizes,afterPR.paths, budgetT, start, end);
        // repairAllRoutes(costs, afterPR.paths, afterPR.pathCosts, prizes, budgetT);
        std::cout << "construct prize after local search: " << afterPR.score << std::endl;

        // ---- Update elite pool exactly per paper (purge/add/replace; update pool's counters) ----
        elites.updateElitePool(afterPR);

        // ---- Keep iteration incumbent and global best ----
        const Solution& incumbent = (afterPR.score >= current.score ? afterPR : current);
        if (incumbent.score > bestGlobal.score + 1e-12) {
            bestGlobal = incumbent;
            no_improve_streak = 0;
        } else {
            ++no_improve_streak;
        }

        // Optional early stop (paper's “no-improvement” rule)
        if (no_improve_streak >= STOP_NO_IMPROVE) break;
    }

    // twoOptAllRoutes(costs, bestGlobal.paths, bestGlobal.pathCosts, budgetT);
    // insertHighValueCheapestMulti(costs, prizes, bestGlobal.paths, budgetT, start, end);

    bestGlobal.score = multiPathsSumPrize(prizes, bestGlobal.paths);
    return {bestGlobal.paths, bestGlobal.score};
}


std::pair<std::vector<std::vector<int>>, double> PROptimization_parallel(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               int start, int end) {
    // --- Tunables ---
    const int   POOL_SIZE        = 5;     // paper commonly uses 5
    const int   MAX_ITERS        = 300;  // hard cap
    const int   STOP_NO_IMPROVE  = 10;   // paper's "slow PR" ~300; use 10 for "fast"
    const double SIM_SKIP        = 0.98;  // skip PR if similarity >= this
    int num_threads = 40;
    std::vector<Solution> pr_solutions(num_threads);

    ElitePool elites(POOL_SIZE);

    // RNG for GRASP greediness alpha
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // Track global best
    Solution bestGlobal{ {}, {}, -std::numeric_limits<double>::infinity() };
    int no_improve_streak = 0;

    for (int it = 0; it < MAX_ITERS; ++it) {
        #pragma omp parallel for num_threads(num_threads)
        for (int thread = 0; thread < num_threads; thread++) {
            // ---- Construct (randomized GRASP) ----
            const double greediness = U01(rng); 
            auto construct = construct_paths(costs, prizes, start, end, nPath, budgetT, greediness);
            auto pathCosts = multiPathsCost(costs, construct);
            auto pathPrize = multiPathsSumPrize(prizes, construct);
            // return {construct, pathPrize};
            localSearchLoop(costs, construct, pathCosts, prizes, budgetT);
            pathCosts = multiPathsCost(costs, construct); // re-read after LS

            Solution current{ construct, pathCosts, multiPathsSumPrize(prizes, construct) };
            Solution afterPR = linkToElites(costs, prizes, current, elites.pool, budgetT, SIM_SKIP);
            repairMultiPaths(costs,prizes,afterPR.paths, budgetT, start, end);
            pr_solutions[thread] = afterPR;
            // pr_soluions[thread] = (afterPR.score >= current.score ? afterPR : current);
            // repairAllRoutes(costs, afterPR.paths, afterPR.pathCosts, prizes, budgetT);
        }

        for(auto& s : pr_solutions)
            elites.updateElitePool(s);

        // ---- Keep iteration incumbent and global best ----
        int best_idx = std::distance(
                    pr_solutions.begin(),
                    std::max_element(pr_solutions.begin(), pr_solutions.end(),
                                    [](const Solution& a, const Solution& b) {
                                        return a.score < b.score;
                                    }));

        Solution& incumbent = pr_solutions[best_idx];
        if (incumbent.score > bestGlobal.score + 1e-12) {
            bestGlobal = incumbent;
            no_improve_streak = 0;
        } else {
            ++no_improve_streak;
        }

        // Optional early stop (paper's “no-improvement” rule)
        if (no_improve_streak >= STOP_NO_IMPROVE) break;
    }

    // twoOptAllRoutes(costs, bestGlobal.paths, bestGlobal.pathCosts, budgetT);
    // insertHighValueCheapestMulti(costs, prizes, bestGlobal.paths, budgetT, start, end);

    bestGlobal.score = multiPathsSumPrize(prizes, bestGlobal.paths);
    return {bestGlobal.paths, bestGlobal.score};
}


std::pair<std::vector<std::vector<int>>, double> PROptimization_parallel(const std::vector<std::vector<double>>& costs,
                               const std::vector<double>& prizes,
                               double budgetT, int nPath,
                               const std::vector<int>& starts, int end) {
    // --- Tunables ---
    const int   POOL_SIZE        = 5;     // paper commonly uses 5
    const int   MAX_ITERS        = 300;  // hard cap
    const int   STOP_NO_IMPROVE  = 10;   // paper's "slow PR" ~300; use 10 for "fast"
    const double SIM_SKIP        = 0.98;  // skip PR if similarity >= this
    int num_threads = 20;
    std::vector<Solution> pr_solutions(num_threads);

    ElitePool elites(POOL_SIZE);

    // RNG for GRASP greediness alpha
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<double> U01(0.0, 1.0);

    // Track global best
    Solution bestGlobal{ {}, {}, -std::numeric_limits<double>::infinity() };
    int no_improve_streak = 0;

    for (int it = 0; it < MAX_ITERS; ++it) {
        #pragma omp parallel for num_threads(num_threads)
        for (int thread = 0; thread < num_threads; thread++) {
            // ---- Construct (randomized GRASP) ----
            const double greediness = U01(rng); 
            auto construct = construct_paths_multiS(costs, prizes, starts, end, nPath, budgetT, greediness);
            auto pathCosts = multiPathsCost(costs, construct);
            auto pathPrize = multiPathsSumPrize(prizes, construct);
            // return {construct, pathPrize};
            localSearchLoop(costs, construct, pathCosts, prizes, budgetT);
            pathCosts = multiPathsCost(costs, construct); // re-read after LS

            Solution current{ construct, pathCosts, multiPathsSumPrize(prizes, construct) };
            Solution afterPR = linkToElites(costs, prizes, current, elites.pool, budgetT, SIM_SKIP);
            repairMultiPaths(costs,prizes,afterPR.paths, budgetT, starts, end);
            pr_solutions[thread] = afterPR;
            // pr_soluions[thread] = (afterPR.score >= current.score ? afterPR : current);
            // repairAllRoutes(costs, afterPR.paths, afterPR.pathCosts, prizes, budgetT);
        }

        for(auto& s : pr_solutions)
            elites.updateElitePool(s);

        // ---- Keep iteration incumbent and global best ----
        int best_idx = std::distance(
                    pr_solutions.begin(),
                    std::max_element(pr_solutions.begin(), pr_solutions.end(),
                                    [](const Solution& a, const Solution& b) {
                                        return a.score < b.score;
                                    }));

        Solution& incumbent = pr_solutions[best_idx];
        if (incumbent.score > bestGlobal.score + 1e-12) {
            bestGlobal = incumbent;
            no_improve_streak = 0;
        } else {
            ++no_improve_streak;
        }

        // Optional early stop (paper's “no-improvement” rule)
        if (no_improve_streak >= STOP_NO_IMPROVE) break;
    }

    // twoOptAllRoutes(costs, bestGlobal.paths, bestGlobal.pathCosts, budgetT);
    // insertHighValueCheapestMulti(costs, prizes, bestGlobal.paths, budgetT, start, end);

    bestGlobal.score = multiPathsSumPrize(prizes, bestGlobal.paths);
    return {bestGlobal.paths, bestGlobal.score};
}
