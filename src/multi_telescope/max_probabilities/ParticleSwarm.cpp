#include "ReadData.h"
#include "PathSplitter.h" 
#include "ParticleSwarm.h"
#include "LocalSearch.h"

#include <algorithm>
#include <numeric>
#include <random>
#include <utility>
#include <vector>
#include <functional>
#include <iostream>


std::vector<int> order_from_keys_excluding_st(const std::vector<double>& keys, int N) {
    std::vector<int> base(N-2); std::iota(base.begin(), base.end(), 1);
    std::vector<int> idx(N-2);  std::iota(idx.begin(),  idx.end(),  0);
    std::stable_sort(idx.begin(), idx.end(),
        [&](int a, int b){ return keys[a] < keys[b]; });
    std::vector<int> order; order.reserve(N-2);
    for (int id : idx) order.push_back(base[id]);
    return order;
}

//build node list and order from keys
std::vector<int> customers_excluding_st(const int N,
                       const std::vector<int>& starts,
                       int t) {
    std::vector<char> banned(N, 0);
    for (int s : starts) if (s >= 0 && s < N) banned[s] = 1;
    if (t >= 0 && t < N) banned[t] = 1;

    std::vector<int> cust; cust.reserve(N);
    for (int v = 0; v < N; ++v) if (!banned[v]) cust.push_back(v);
    return cust;
}

static inline std::vector<int>
order_from_keys_on_customers(const std::vector<double>& keys,
                             const std::vector<int>& customers)
{
    const int D = (int)customers.size();
    std::vector<int> idx(D); std::iota(idx.begin(), idx.end(), 0);
    std::stable_sort(idx.begin(), idx.end(),
        [&](int a, int b){ return keys[a] < keys[b]; });

    std::vector<int> perm; perm.reserve(D);
    for (int id : idx) perm.push_back(customers[id]);
    return perm;
}

//greedy decoder for multiple starts + common t
std::pair<std::vector<std::vector<int>>, double>
decode_greedy_multi(const std::vector<std::vector<double>>& costs,
                    const std::vector<double>& prizes,
                    const std::vector<int>& starts,
                    int t, double B,
                    const std::vector<int>& perm) {
    std::vector<int> avail = starts;             
    std::vector<std::vector<int>> routes; routes.reserve(starts.size());
    std::vector<char> used(costs.size(), 0); 

    auto take_closest_start_pos = [&](int v_first)->int {
        int best = -1; double bc = std::numeric_limits<double>::infinity();
        for (int i = 0; i < (int)avail.size(); ++i) {
            double c = costs[avail[i]][v_first];
            if (c < bc) { bc = c; best = i; }
        }
        return best;
    };

    int i = 0; 
    const int M = (int)perm.size();

    // Greedily build routes while we have node and free starts
    while (i < M && !avail.empty()) {
        while (i < M && used[perm[i]]) ++i;
        if (i >= M) break;

        const int v_first = perm[i];
        int pos = take_closest_start_pos(v_first);
        if (pos < 0) break;

        const int s = avail[pos];
        const double s_to_first = costs[s][v_first];
        const double first_to_t = costs[v_first][t];

        if (s_to_first + first_to_t > B) {
            ++i;
            continue;
        }

        // Grow segment as far as budget allows
        std::vector<int> segment;
        segment.push_back(v_first);
        used[v_first] = 1;

        double acc_edges = 0.0;
        int j = i;

        while (j + 1 < M) {
            int u = perm[j];
            int v = perm[j + 1];
            if (used[v]) { ++j; continue; }

            double next_edge = costs[u][v];
            double v_to_t    = costs[v][t];
            if (s_to_first + (acc_edges + next_edge) + v_to_t <= B) {
                acc_edges += next_edge;
                segment.push_back(v);
                used[v] = 1;
                ++j;
            } else {
                break;
            }
        }

        std::vector<int> route;
        route.reserve(2 + (int)segment.size());
        route.push_back(s);
        for (int v : segment) route.push_back(v);
        route.push_back(t);
        routes.push_back(std::move(route));

        avail.erase(avail.begin() + pos);

        i = j + 1;
    }

    for (int s : avail) {
        if (costs[s][t] <= B) routes.push_back({s, t});
        else routes.push_back({s}); 
    }

    // Score = sum prizes
    double total_prize = 0.0;
    for (int v = 0; v < (int)used.size(); ++v) if (used[v]) total_prize += prizes[v];

    return {routes, total_prize};
}

//PSO multiple starts + common t
std::pair<std::vector<std::vector<int>>, double>
run_pso_team_op_multi(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes,
                    double budget,
                    int m_routes,
                    const std::vector<int>& starts,
                    int t, 
                    int swarm_size, int iters,
                    double vmax, unsigned seed) {

    const int N = (int)prizes.size();

    // Feasibility check for padding routes
    for (int s : starts) {
        if (s < 0 || s >= N) {
            std::cerr << "[warn] start id out of range: " << s << "\n";
            return {};
        }
    }
    if (t < 0 || t >= N) {
        std::cerr << "[warn] t id out of range\n";
        return {};
    }

    // Customer set (exclude starts and t)
    const std::vector<int> customers = customers_excluding_st(N, starts, t);
    const int D = (int)customers.size();

    if (D == 0) {
        std::vector<std::vector<int>> routes;
        routes.reserve(starts.size());
        for (int s : starts) {
            if (costs[s][t] <= budget) routes.push_back({s, t});
            else routes.push_back({s}); 
        }
        return {routes, 0.0};
    }

    // PSO params (Clerc–Kennedy)
    const double c1 = 2.05, c2 = 2.05, chi = 0.72984;

    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> U01(0.0, 1.0);
    std::uniform_real_distribution<double> Uv(-vmax, vmax);

    auto decode = [&](const std::vector<double>& keys)
        -> std::pair<std::vector<std::vector<int>>, double>
    {
        std::vector<int> perm = order_from_keys_on_customers(keys, customers);
        return decode_greedy_multi(costs, prizes, starts, t, budget, perm);
    };

    // Swarm arrays
    std::vector<std::vector<double>> X(swarm_size, std::vector<double>(D));
    std::vector<std::vector<double>> V(swarm_size, std::vector<double>(D));
    std::vector<std::vector<double>> Pbest(swarm_size, std::vector<double>(D));
    std::vector<double> PbestFit(swarm_size, -1e300);

    std::vector<std::vector<int>> gbest_routes;
    double gbest_fit = -1e300;
    std::vector<double> gbest(D, 0.0);

    // Init
    for (int i = 0; i < swarm_size; ++i) {
        for (int d = 0; d < D; ++d) {
            X[i][d] = U01(rng);
            V[i][d] = 0.1 * Uv(rng);
        }
        auto [routes, fit] = decode(X[i]);
        Pbest[i] = X[i];
        PbestFit[i] = fit;
        if (fit > gbest_fit) { gbest_fit = fit; gbest = X[i]; gbest_routes = routes; }
    }

    // Iterate
    for (int it = 0; it < iters; ++it) {
        for (int i = 0; i < swarm_size; ++i) {
            // velocity & position
            for (int d = 0; d < D; ++d) {
                double r1 = U01(rng), r2 = U01(rng);
                double vv = V[i][d] + c1*r1*(Pbest[i][d] - X[i][d])
                                     + c2*r2*(gbest[d]     - X[i][d]);
                vv = chi * vv;
                if (vv >  vmax) vv =  vmax;
                if (vv < -vmax) vv = -vmax;
                V[i][d] = vv;

                double xx = X[i][d] + vv;
                if (xx < 0.0) xx = 0.0;
                if (xx > 1.0) xx = 1.0;
                X[i][d] = xx;
            }

            // evaluate via greedy decoder
            auto [routes, fit] = decode(X[i]);

            if (fit > PbestFit[i]) { PbestFit[i] = fit; Pbest[i] = X[i]; }
            if (fit > gbest_fit)   { gbest_fit   = fit; gbest   = X[i]; gbest_routes = routes; }
        }
    }

    for (auto& r : gbest_routes) {
        while (r.size() > 2 && two_opt_st(r, costs)) {}
    }

    return {gbest_routes, gbest_fit};
}

