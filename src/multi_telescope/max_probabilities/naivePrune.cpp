#include <vector>
#include <limits>
#include <cmath>
#include <algorithm>

static inline double cycle_cost(const std::vector<std::vector<double>>& C,
                                const std::vector<int>& cyc) {
    if (cyc.size() < 2) return 0.0;
    double s = 0.0;
    for (size_t i = 0; i + 1 < cyc.size(); ++i) s += C[cyc[i]][cyc[i+1]];
    // close the cycle
    s += C[cyc.back()][cyc.front()];
    return s;
}

// Prune a *cycle* by removing contiguous segments (possibly wrapping head+tail)
// to bring cost under 'budget'. 'cyc' is a list of distinct vertices in cycle order
// (no duplicated last==first needed). 'total_cost' is updated to the final cycle cost.
// If the cycle cannot be reduced further (no positive-savings segment), the function stops.
void pruneCycleBySegmentRatio(const std::vector<std::vector<double>>& costs,
                              const std::vector<double>& prizes,
                              std::vector<int>& cyc,
                              double budget,
                              double& total_cost)
{
    const double eps = 1e-12;
    if (cyc.size() < 2) { total_cost = 0.0; return; }

    // Ensure total_cost is initialized
    if (!(total_cost > 0.0) || !std::isfinite(total_cost)) {
        total_cost = cycle_cost(costs, cyc);
    }

    // Repeat until within budget or no positive-savings removal exists
    while (total_cost > budget + eps && cyc.size() > 2) {
        const size_t n = cyc.size();

        // Duplicate the vertex list to linearize wrap-around segments
        std::vector<int> W(2 * n);
        for (size_t k = 0; k < 2 * n; ++k) W[k] = cyc[k % n];

        // Prefix sums of prizes over W
        std::vector<double> prefP(2 * n + 1, 0.0);
        for (size_t k = 0; k < 2 * n; ++k) {
            prefP[k + 1] = prefP[k] + prizes[W[k]];
        }

        // Prefix sums of *edge* costs along W (k -> k+1), up to 2n-2 edges
        // We only need edges within windows of length <= n-1.
        std::vector<double> prefE(2 * n, 0.0);
        for (size_t k = 0; k + 1 < 2 * n; ++k) {
            prefE[k + 1] = prefE[k] + costs[W[k]][W[k + 1]];
        }

        auto segPrize = [&](size_t i2, size_t j2) -> double {
            // prizes of nodes W[i2..j2], inclusive
            return prefP[j2 + 1] - prefP[i2];
        };
        auto internalSegEdgeCost = [&](size_t i2, size_t j2) -> double {
            // sum of edges (W[i2]→W[i2+1]) + ... + (W[j2-1]→W[j2]); zero if i2==j2
            return prefE[j2] - prefE[i2];
        };
        auto segDeltaCost = [&](size_t i2, size_t j2) -> double {
            // Remove W[i2..j2] (length L = j2-i2+1, with L <= n-1), reconnect prev→next in the original cycle.
            const int prev = cyc[( (i2 % n) + n - 1 ) % n];       // cyc index wraps
            const int next = cyc[( (j2 % n) + 1 ) % n];
            const int first = W[i2];
            const int last  = W[j2];

            const double old_boundary = costs[prev][first] + costs[last][next];
            const double old_internal = internalSegEdgeCost(i2, j2);
            const double new_edge     = costs[prev][next];

            return (old_boundary + old_internal) - new_edge;  // savings
        };

        const double over = total_cost - budget;

        // Search the best (smallest prize / delta) removable segment with positive savings
        double best_ratio = std::numeric_limits<double>::infinity();
        size_t best_i2 = 0, best_j2 = 0;
        bool found = false;

        // i2 is the start index in W for the segment; j2 goes up to i2 + (n - 2)
        for (size_t i2 = 0; i2 < n; ++i2) {
            const size_t j2_max = i2 + (n - 2); // cannot remove all vertices
            double best_local = std::numeric_limits<double>::infinity();
            size_t best_local_j2 = i2;
            bool local_found = false;

            for (size_t j2 = i2; j2 <= j2_max; ++j2) {
                const double delta = segDeltaCost(i2, j2);
                if (delta <= eps) continue;  // no savings or negative savings

                const double pr = segPrize(i2, j2);
                const double ratio = pr / delta;

                if (ratio < best_local) {
                    best_local = ratio;
                    best_local_j2 = j2;
                    local_found = true;
                }
                // Early stop if this growing segment can already cover the overage
                if (delta >= over - eps) break;
            }

            if (local_found && best_local < best_ratio) {
                best_ratio = best_local;
                best_i2 = i2;
                best_j2 = best_local_j2;
                found = true;
            }
        }

        if (!found || !std::isfinite(best_ratio)) {
            // No positive-savings removal exists; stop (cannot meet budget)
            break;
        }

        // Apply removal of W[best_i2..best_j2] in the *cyc* (mod n) indexing
        const double delta_cost = segDeltaCost(best_i2, best_j2);
        if (delta_cost <= eps) break; // defensive

        // Build a mask of removed vertices in cyc-index space
        std::vector<char> removed(n, 0);
        for (size_t k = best_i2; k <= best_j2; ++k) {
            removed[k % n] = 1;
        }
        const size_t new_n = n - (best_j2 - best_i2 + 1);

        // Reconstruct the cycle starting from 'next' = cyc[(best_j2+1)%n]
        std::vector<int> new_cyc;
        new_cyc.reserve(new_n);
        size_t cur = (best_j2 + 1) % n;
        for (size_t cnt = 0; cnt < new_n; ++cnt) {
            if (!removed[cur]) new_cyc.push_back(cyc[cur]);
            cur = (cur + 1) % n;
        }

        cyc.swap(new_cyc);
        total_cost -= delta_cost;

        if (cyc.size() < 2) { total_cost = 0.0; break; }
    }

    // Optional: resync cost exactly (defensive)
    double exact = cycle_cost(costs, cyc);
    if (std::abs(exact - total_cost) > 1e-9) {
        total_cost = exact;
    }
}