#include "ReadData.h"
#include "ParticleSwarmInspired.h"

#include <vector>
#include <functional>
#include <random>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <cmath>
#include <iostream>

// // O(m·n) interval-graph split over saturated tours (Dang–Guibadj–Moukrim, 2013)
// std::pair<std::vector<std::vector<int>>, double>
// split_interval_graph_eval(const Instance& I, double B,
//                           const std::vector<int>& perm, int m_routes)
// {
//     const int n = (int)perm.size();
//     std::vector<std::vector<int>> routes;
//     double total_prize = 0.0;
//     if (n == 0 || m_routes <= 0) return {routes, total_prize};

//     // 1) Longest feasible prefix (saturated tour) starting at each i
//     std::vector<int>    lmax(n, -1);   // length of [i..last]
//     std::vector<int>    succ(n, n);    // next start index or n (past-the-end)
//     std::vector<double> P(n, 0.0);     // prize of saturated block i

//     for (int i = 0; i < n; ++i) {
//         double acc = I.costs[I.s][perm[i]];
//         int j = i, last = i - 1;
//         while (j < n && acc + I.costs[perm[j]][I.t] <= B + 1e-12) {
//             last = j;
//             ++j;
//             if (j < n) acc += I.costs[perm[j-1]][perm[j]];
//         }
//         if (last >= i) {
//             lmax[i] = last - i;
//             succ[i] = last + 1;
//             double pr = 0.0; for (int k = i; k <= last; ++k) pr += I.prizes[perm[k]];
//             P[i] = pr;
//         } else {
//             lmax[i] = -1;
//             succ[i] = i + 1;
//             P[i]    = -1e300; // NEG → never “take”
//         }
//     }

//     // 2) DP: C[i][k] = max prize using k routes starting at index i
//     const double NEG = -1e300;
//     std::vector<std::vector<double>> C(n + 1, std::vector<double>(m_routes + 1, NEG));
//     for (int k = 0; k <= m_routes; ++k) C[n][k] = 0.0;

//     for (int i = n - 1; i >= 0; --i) {
//         C[i][0] = 0.0;
//         for (int k = 1; k <= m_routes; ++k) {
//             double take = (lmax[i] >= 0) ? (P[i] + C[succ[i]][k - 1]) : NEG;
//             double skip = C[i + 1][k];
//             C[i][k] = (take > skip) ? take : skip;
//         }
//     }

//     // 3) Backtrack
//     total_prize = C[0][m_routes];
//     int i = 0, k = m_routes;
//     while (i < n && k > 0) {
//         if (C[i][k] == C[i + 1][k]) { ++i; continue; }  // skip i
//         const int last = i + lmax[i];                   // take block [i..last]
//         std::vector<int> r; r.reserve(last - i + 3);
//         r.push_back(I.s);
//         for (int p = i; p <= last; ++p) r.push_back(perm[p]);
//         r.push_back(I.t);
//         routes.push_back(std::move(r));
//         i = succ[i]; --k;
//     }
//     return {routes, total_prize};
// }

std::vector<int> idch_heuristic(const Instance& I, std::mt19937& rng) {
    // Giant tour = all customers except s,t; sort by prize, add a little randomness.
    std::vector<int> perm; perm.reserve(I.size - 2);
    for (int v = 1; v + 1 < I.size; ++v) perm.push_back(v);
    std::stable_sort(perm.begin(), perm.end(),
        [&](int a, int b){ return I.prizes[a] > I.prizes[b]; });
    int k = std::min<int>(6, (int)perm.size());
    std::shuffle(perm.end() - k, perm.end(), rng);
    return perm;
}


// PSOiA crossover-like update (Sec. 3.5)
// std::vector<int> psoa_update(const std::vector<int>& cur,
//                                     const std::vector<int>& pbest,
//                                     const std::vector<int>& gbest,
//                                     double w, double c1, double c2,
//                                     std::mt19937& rng)
// {
//     const int n = (int)cur.size();
//     std::uniform_real_distribution<double> U01(0.0,1.0);
//     const double r1 = U01(rng), r2 = U01(rng);

//     const int l_cur = (int)std::floor(std::clamp(w, 0.0, 1.0) * n);
//     int l_p = 0, l_g = 0;
//     if (c1*r1 + c2*r2 > 1e-12) {
//         int rem = n - l_cur;
//         l_p = (int)std::floor(rem * (c1*r1) / (c1*r1 + c2*r2));
//         l_g = rem - l_p;
//     } else {
//         l_p = 0; l_g = n - l_cur;
//     }

//     auto extract = [&](const std::vector<int>& src, int L,
//                        std::unordered_set<int>& used)->std::vector<int>
//     {
//         std::vector<int> out;
//         if (L <= 0) return out;
//         std::uniform_int_distribution<int> Upos(0, n-1);
//         int r = Upos(rng);
//         for (int i = r; i < n && (int)out.size() < L; ++i)
//             if (!used.count(src[i])) { used.insert(src[i]); out.push_back(src[i]); }
//         for (int i = r - 1; i >= 0 && (int)out.size() < L; --i)
//             if (!used.count(src[i])) { used.insert(src[i]); out.insert(out.begin(), src[i]); }
//         return out;
//     };

//     std::unordered_set<int> used;
//     std::vector<int> A = extract(cur,   l_cur, used, rng);
//     std::vector<int> B = extract(pbest, l_p,   used, rng);
//     std::vector<int> C = extract(gbest, l_g,   used, rng);

//     std::array<std::vector<int>*,3> chunks = {&A,&B,&C};
//     std::shuffle(chunks.begin(), chunks.end(), rng);

//     std::vector<int> child; child.reserve(n);
//     for (auto* v : chunks) child.insert(child.end(), v->begin(), v->end());
//     for (int v : cur) if (!used.count(v)) child.push_back(v); // safety
//     return child;
// }

std::vector<int> psoa_update(const std::vector<int>& cur,
                             const std::vector<int>& pbest,
                             const std::vector<int>& gbest,
                             double w, double c1, double c2,
                             std::mt19937& rng) {
    const int n = (int)cur.size();
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    const double r1 = U01(rng), r2 = U01(rng);

    const int l_cur = (int)std::floor(std::max(0.0, std::min(1.0, w)) * n);
    int l_p = 0, l_g = 0;
    if (c1*r1 + c2*r2 > 1e-12) {
        const int rem = n - l_cur;
        l_p = (int)std::floor(rem * (c1*r1) / (c1*r1 + c2*r2));
        l_g = rem - l_p;
    } else {
        l_p = 0; l_g = n - l_cur;
    }

    // Capture rng from outer scope; no rng parameter here
    auto extract = [&](const std::vector<int>& src, int L,
                       std::unordered_set<int>& used) -> std::vector<int>
    {
        std::vector<int> out;
        if (L <= 0) return out;
        std::uniform_int_distribution<int> Upos(0, n - 1);
        int r = Upos(rng);
        for (int i = r; i < n && (int)out.size() < L; ++i)
            if (!used.count(src[i])) { used.insert(src[i]); out.push_back(src[i]); }
        for (int i = r - 1; i >= 0 && (int)out.size() < L; --i)
            if (!used.count(src[i])) { used.insert(src[i]); out.insert(out.begin(), src[i]); }
        return out;
    };

    std::unordered_set<int> used;
    // Call with three arguments (rng is captured)
    std::vector<int> A = extract(cur,   l_cur, used);
    std::vector<int> B = extract(pbest, l_p,   used);
    std::vector<int> C = extract(gbest, l_g,   used);

    std::array<std::vector<int>*, 3> chunks = { &A, &B, &C };
    std::shuffle(chunks.begin(), chunks.end(), rng);

    std::vector<int> child; child.reserve(n);
    for (auto* v : chunks) child.insert(child.end(), v->begin(), v->end());
    for (int v : cur) if (!used.count(v)) child.push_back(v); // fill any missing
    return child;
}

// Local search stub (can replace with your own)
void local_search_stub(const Instance& I, std::vector<int>& perm) {
    (void)I; (void)perm; 
}


std::pair<std::vector<std::vector<int>>, double>
run_psoia_team_op(const Instance& I, double B, int m_routes,
                  int swarm_size, int k_stop, unsigned seed,
                  const SplitterFn& splitter,
                  double ph, double d,
                  double w0, double c1, double c2) {
    // Optional vertices are 1..size-2
    const int N = I.size;
    const int nOpt = std::max(0, N - 2);
    std::vector<int> base(nOpt); std::iota(base.begin(), base.end(), 1);

    // quick feasibility
    if (I.costs[I.s][I.t] > B + 1e-12) {
        std::cerr << "[warn] s->t cost exceeds B; no feasible route.\n";
        return {{}, 0.0};
    }

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U01(0.0,1.0);

    struct Particle {
        std::vector<int> pos;
        std::vector<int> lbest;
        double lbest_fit = -1e300;
        std::vector<std::vector<int>> lbest_routes;
        double lbest_travel = 0.0; // sum of travel costs (for similarity)
    };
    auto routes_travel_sum = [&](const std::vector<std::vector<int>>& R)->double {
        double c=0.0;
        for (const auto& r : R)
            for (size_t i=0;i+1<r.size();++i) c += I.costs[r[i]][r[i+1]];
        return c;
    };

    // itermax = k_stop * nOpt / m_routes  (paper)
    const int itermax = std::max(1, k_stop * std::max(1, nOpt) / std::max(1, m_routes));
    std::cout << "itermax : " << itermax << "\n";

    // Init
    std::vector<Particle> S(swarm_size);
    for (int i = 0; i < swarm_size; ++i) {
        S[i].pos = base;
        std::shuffle(S[i].pos.begin(), S[i].pos.end(), rng);
        if (i < std::min(5, swarm_size)) S[i].pos = idch_heuristic(I, rng);
        auto [rts, fit] = splitter(I, B, S[i].pos, m_routes);
        S[i].lbest = S[i].pos; S[i].lbest_fit = fit; S[i].lbest_routes = rts;
        S[i].lbest_travel = routes_travel_sum(rts);
    }

    // Global best = best lbest
    auto best_it = std::max_element(S.begin(), S.end(),
        [](const Particle& a, const Particle& b){ return a.lbest_fit < b.lbest_fit; });
    std::vector<int> gbest = best_it->lbest;
    auto [gbest_routes, gbest_fit] = splitter(I, B, gbest, m_routes);
    double gbest_travel = routes_travel_sum(gbest_routes);

    int no_improve = 0;
    double w = w0;

    while (no_improve < itermax) {
        bool any_improve = false;
        const double pm = 1.0 - double(no_improve) / double(itermax);

        // Find worst lbest (for Rule 1)
        int worst_idx = int(std::min_element(S.begin(), S.end(),
            [](const Particle& a, const Particle& b){ return a.lbest_fit < b.lbest_fit; }) - S.begin());

        for (int x = 0; x < swarm_size; ++x) {
            // Move-out
            if (U01(rng) < ph) S[x].pos = idch_heuristic(I, rng);

            // PSOiA crossover-like update
            S[x].pos = psoa_update(S[x].pos, S[x].lbest, gbest, w, c1, c2, rng);

            // Local search (optional)
            if (U01(rng) < pm) local_search_stub(I, S[x].pos);

            // Evaluate
            auto [rts, fit] = splitter(I, B, S[x].pos, m_routes);
            double trav = routes_travel_sum(rts);

            // Similarity-based lbest update
            if (fit > S[worst_idx].lbest_fit + 1e-12) {
                bool replaced = false;
                for (int y = 0; y < swarm_size; ++y) {
                    if (std::abs(S[y].lbest_fit - fit) <= 1e-12 &&
                        std::abs(S[y].lbest_travel - trav) <= d) {
                        S[y].lbest = S[x].pos;
                        S[y].lbest_fit = fit;
                        S[y].lbest_routes = rts;
                        S[y].lbest_travel = trav;
                        replaced = true; any_improve = true;
                        break;
                    }
                }
                if (!replaced) { // replace worst
                    S[worst_idx].lbest = S[x].pos;
                    S[worst_idx].lbest_fit = fit;
                    S[worst_idx].lbest_routes = rts;
                    S[worst_idx].lbest_travel = trav;
                    any_improve = true;
                }
            }
        }

        // Global best
        auto nb = std::max_element(S.begin(), S.end(),
            [](const Particle& a, const Particle& b){ return a.lbest_fit < b.lbest_fit; });
        if (nb->lbest_fit > gbest_fit + 1e-12) {
            gbest = nb->lbest;
            gbest_routes = nb->lbest_routes;
            gbest_fit = nb->lbest_fit;
            gbest_travel = nb->lbest_travel;
            any_improve = true;
        }

        no_improve = any_improve ? 0 : (no_improve + 1);
        w = std::max(0.1, 0.9 * w); // cooling
    }

    return {gbest_routes, gbest_fit};
}
