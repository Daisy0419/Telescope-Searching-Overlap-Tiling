#include "recursive_greedy_op.h"
#include "SubmodularSelection.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <vector>


constexpr long long kInfeasible = -1;
constexpr std::size_t kMemoCap = 500000;

struct MemoKey {
    long long budget;
    int a, b, depth;
    std::vector<char> committed;
    std::vector<char> blocked;

    bool operator==(const MemoKey& o) const {
        return budget == o.budget && a == o.a && b == o.b &&
               depth == o.depth && committed == o.committed &&
               blocked == o.blocked;
    }
};

// Hashes an exact memoization key. Hash collisions are harmless because
// unordered_map also calls MemoKey::operator== before reusing a result.
struct MemoHash {
    std::size_t operator()(const MemoKey& key) const {
        std::size_t seed = 0;
        combine(seed, std::hash<long long>{}(key.budget));
        combine(seed, std::hash<int>{}(key.a));
        combine(seed, std::hash<int>{}(key.b));
        combine(seed, std::hash<int>{}(key.depth));
        for (char value : key.committed)
            combine(seed, std::hash<unsigned>{}(
                static_cast<unsigned char>(value)));
        for (char value : key.blocked)
            combine(seed, std::hash<unsigned>{}(
                static_cast<unsigned char>(value)));
        return seed;
    }

private:
    static void combine(std::size_t& seed, std::size_t value) {
        // A conventional polynomial combination is sufficient here because
        // unordered_map still verifies equality after a hash match.
        seed = seed * 31U + value;
    }
};

class RecursiveGreedySolver {
public:
    RecursiveGreedySolver(const std::vector<std::vector<double>>& costs,
                          const std::vector<double>& pixel_probs,
                          const std::vector<std::vector<int>>& member_pixels,
                          const CoverageOracle& coverage_oracle,
                          int start, int end, double budget,
                          const RecursiveGreedyParams& params)
        : costs_(costs),
          member_pixels_(member_pixels),
          coverage_oracle_(coverage_oracle),
          params_(params),
          node_count_(static_cast<int>(costs.size())),
          pixel_count_(static_cast<int>(pixel_probs.size())),
          start_(start),
          end_(end) {

        // Arc lengths round UP and the budget rounds DOWN, so anything
        // feasible under the quantised lengths is feasible under the real
        // costs. ceil() also preserves the triangle inequality.
        const long double q =
            static_cast<long double>(params_.time_quantum);
        length_.assign(node_count_, std::vector<long long>(node_count_, 0));
        for (int i = 0; i < node_count_; ++i) {
            for (int j = 0; j < node_count_; ++j) {
                if (i == j) continue;
                const long double scaled =
                    static_cast<long double>(costs_[i][j]) / q;
                if (scaled > static_cast<long double>(
                                 std::numeric_limits<long long>::max() / 4))
                    throw std::overflow_error(
                        "recursive greedy: quantised edge cost overflow.");
                length_[i][j] =
                    static_cast<long long>(std::ceil(scaled));
            }
        }

        const long double scaled_budget =
            static_cast<long double>(budget) / q;
        if (scaled_budget > static_cast<long double>(
                                std::numeric_limits<long long>::max() / 4))
            throw std::overflow_error(
                "recursive greedy: quantised budget overflow.");
        budget_ = static_cast<long long>(std::floor(scaled_budget));

        // Round the PIXEL weights rather than f. Weighted coverage over
        // integer weights is still exactly monotone submodular, so the
        // paper's implicit-rounding argument for a general oracle is not
        // needed here.
        double total = 0.0;
        double smallest_positive = std::numeric_limits<double>::infinity();
        for (double p : pixel_probs) {
            total += p;
            if (p > 0.0) smallest_positive = std::min(smallest_positive, p);
        }
        double scale = (total > 0.0)
            ? static_cast<double>(params_.reward_levels) / total
            : 0.0;

        // A HEALPix map spreads probability over thousands of pixels, so any
        // single pixel is a tiny share of the total. Scaling by the total
        // alone rounds every weight to zero, which silently flattens the
        // objective and makes the trivial s-t arc optimal. Raise the scale
        // until the smallest positive pixel is worth at least one unit,
        // then clamp so the total stays well inside int64.
        if (std::isfinite(smallest_positive) && smallest_positive * scale < 1.0)
            scale = 1.0 / smallest_positive;
        if (total > 0.0 && total * scale > 1e15) scale = 1e15 / total;

        weight_.assign(pixel_count_, 0);
        reward_scaling_exact_ = true;
        for (int p = 0; p < pixel_count_; ++p) {
            const double scaled_probability = pixel_probs[p] * scale;
            weight_[p] = std::max<long long>(
                0, static_cast<long long>(std::llround(scaled_probability)));
            const double rounding_tolerance =
                1e-12 * std::max(1.0, std::abs(scaled_probability));
            if (std::abs(scaled_probability -
                         static_cast<double>(weight_[p])) >
                rounding_tolerance)
                reward_scaling_exact_ = false;
        }
        reward_units_ = 0;
        for (long long w : weight_) reward_units_ += w;
        if (reward_units_ <= 0)
            throw std::invalid_argument(
                "recursive greedy: all pixel weights rounded to zero; the "
                "objective would be flat. Check that pixel_probs is nonzero.");

        cover_count_.assign(pixel_count_, 0);
        committed_count_.assign(node_count_, 0);
        committed_.assign(node_count_, 0);
        blocked_.assign(node_count_, 0);

        deadline_set_ = params_.time_limit > 0.0;
        started_ = std::chrono::steady_clock::now();
    }

    long long calls() const { return calls_; }
    long long rewardUnits() const { return reward_units_; }
    bool rewardScalingExact() const { return reward_scaling_exact_; }
    bool timedOut() const { return timed_out_; }

    double realCost(const std::vector<int>& path) const {
        double c = 0.0;
        for (std::size_t k = 0; k + 1 < path.size(); ++k)
            c += costs_[path[k]][path[k + 1]];
        return c;
    }

    std::vector<int> solve(int depth) {
        long long value = 0;
        return recurse(start_, end_, budget_, depth, value);
    }

    double trueReward(const std::vector<int>& path) const {
        return coverage_oracle_.rewardOf(path);
    }

    // RG's line-1 feasibility test assumes no s-t path beats the direct arc,
    // and the forbidden-node pruning below needs shortcutting to be free.
    // Both are the paper's WLOG triangle inequality; the sink padding is what
    // restores it on this graph.
    bool verifyMetric(double tol = 1e-6) const {
        for (int i = 0; i < node_count_; ++i)
            for (int k = 0; k < node_count_; ++k) {
                if (i == k) continue;
                for (int j = 0; j < node_count_; ++j) {
                    if (j == i || j == k) continue;
                    if (costs_[i][k] > costs_[i][j] + costs_[j][k] + tol)
                        return false;
                }
            }
        return true;
    }

private:
    bool outOfTime() {
        if (!deadline_set_ || timed_out_) return timed_out_;
        if ((++clock_probe_ & 0x3FF) != 0) return false;
        timed_out_ = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started_).count() >
            params_.time_limit;
        return timed_out_;
    }

    void applyCoverage(const std::vector<int>& nodes, long long& gain) {
        gain = 0;
        for (int v : nodes) {
            if (committed_count_[v]++ == 0) committed_[v] = 1;
            for (int p : member_pixels_[v])
                if (cover_count_[p]++ == 0) gain += weight_[p];
        }
    }

    void undoCoverage(const std::vector<int>& nodes) {
        for (int v : nodes) {
            if (committed_count_[v] <= 0)
                throw std::logic_error(
                    "recursive greedy: invalid committed-node rollback.");
            if (--committed_count_[v] == 0) committed_[v] = 0;

            for (int p : member_pixels_[v]) {
                if (cover_count_[p] <= 0)
                    throw std::logic_error(
                        "recursive greedy: invalid pixel-coverage rollback.");
                --cover_count_[p];
            }
        }
    }

    long long marginal(const std::vector<int>& nodes) {
        long long g = 0;
        applyCoverage(nodes, g);
        undoCoverage(nodes);
        return g;
    }
    void setBlocked(int v, bool on) {
        if (static_cast<bool>(blocked_[v]) == on) return;
        blocked_[v] = on ? 1 : 0;
    }

    struct Frame {
        std::vector<int> path;
        long long value = kInfeasible;
    };

    // Minimal budgets at which the value of RG(a, v, ., depth) steps up.
    // This is exactly the set {min b : value(b) >= level} that the paper's
    // a-loop enumerates, found in O(#levels * log B) probes instead of
    // O(A * log B). It relies on the same monotonicity in the budget that
    // makes the paper's binary search valid.
    void collectSplits(int a, int v, long long lo, long long hi, int depth,
                       std::vector<long long>& splits) {
        long long value_lo = 0, value_hi = 0;
        recurse(a, v, lo, depth, value_lo);
        recurse(a, v, hi, depth, value_hi);
        if (value_lo == value_hi || outOfTime()) {
            splits.push_back(lo);
            return;
        }
        if (hi - lo <= 1) {
            splits.push_back(lo);
            splits.push_back(hi);
            return;
        }
        const long long mid = lo + (hi - lo) / 2;
        collectSplits(a, v, lo, mid, depth, splits);
        collectSplits(a, v, mid, hi, depth, splits);
    }

    std::vector<int> middleCandidates(int a, int b, long long budget) {
        std::vector<int> cand;
        for (int v = 0; v < node_count_; ++v) {
            if (v == a || v == b || blocked_[v]) continue;
            if (length_[a][v] + length_[v][b] > budget) continue;
            cand.push_back(v);
        }
        if (params_.candidate_limit <= 0 ||
            static_cast<int>(cand.size()) <= params_.candidate_limit)
            return cand;

        // The one knob that voids the guarantee.
        std::vector<std::pair<double, int>> scored;
        scored.reserve(cand.size());
        for (int v : cand) {
            const double detour = static_cast<double>(
                length_[a][v] + length_[v][b] - length_[a][b]);
            scored.emplace_back(
                static_cast<double>(marginal({v})) / (detour + 1.0), v);
        }
        std::partial_sort(scored.begin(),
                          scored.begin() + params_.candidate_limit, scored.end(),
                          [](const auto& x, const auto& y) { return x.first > y.first; });
        scored.resize(params_.candidate_limit);
        cand.clear();
        for (const auto& s : scored) cand.push_back(s.second);
        return cand;
    }

    // RG-QP(a, b, budget, X, depth);
    std::vector<int> recurse(int a, int b, long long budget, int depth,
                             long long& value) {
        if (a == b) {
            value = marginal({a});
            return {a};
        }
        if (length_[a][b] > budget) {                 // step 1: Infeasible
            value = kInfeasible;
            return {};
        }

        // The exact selected-node and blocked-node states are part of the key.
        // This is deterministic and cannot confuse two residual-reward states.
        const MemoKey key{budget, a, b, depth, committed_, blocked_};
        auto hit = memo_.find(key);
        if (hit != memo_.end()) {
            value = hit->second.value;
            return hit->second.path;
        }
        ++calls_;

        std::vector<int> best{a, b};                  // step 2
        value = marginal(best);                       // step 4

        if (depth > 0 && !outOfTime()) {
            setBlocked(b, true);   // keep this frame's endpoint out of P1
            for (int v : middleCandidates(a, b, budget)) {
                if (outOfTime()) break;

                std::vector<long long> splits;
                // Leave at least the direct v-b cost for the second half.
                const long long largest_left_budget =
                    budget - length_[v][b];
                collectSplits(a, v, length_[a][v], largest_left_budget,
                              depth - 1, splits);
                std::sort(splits.begin(), splits.end());
                splits.erase(std::unique(splits.begin(), splits.end()),
                             splits.end());

                for (long long split : splits) {
                    long long first_value = 0;
                    std::vector<int> first =
                        recurse(a, v, split, depth - 1, first_value);
                    if (first.empty()) continue;

                    // Commit V(P1): it enters X, and everything except the
                    // junction v is barred from the P2 subtree so the halves
                    // can only share v. Lemma 3.2 survives because the optimal
                    // tail can be shortcut around the barred nodes at no extra
                    // length and no loss of marginal reward.
                    long long committed = 0;
                    applyCoverage(first, committed);
                    for (int u : first) if (u != v) setBlocked(u, true);
                    setBlocked(b, false);

                    long long second_value = 0;
                    std::vector<int> second =
                        recurse(v, b, budget - split, depth - 1, second_value);

                    setBlocked(b, true);
                    for (int u : first) if (u != v) setBlocked(u, false);
                    undoCoverage(first);

                    if (second.empty()) continue;
                    std::vector<int> joined = first;
                    joined.insert(joined.end(), second.begin() + 1, second.end());
                    const long long joined_value = marginal(joined);
                    if (joined_value > value) {
                        value = joined_value;
                        best = std::move(joined);
                    }
                }
            }
            setBlocked(b, false);
        }

        if (memo_.size() < kMemoCap && !timed_out_)
            memo_.emplace(key, Frame{best, value});
        return best;
    }

    const std::vector<std::vector<double>>& costs_;
    const std::vector<std::vector<int>>& member_pixels_;
    const CoverageOracle& coverage_oracle_;
    RecursiveGreedyParams params_;

    int node_count_, pixel_count_, start_, end_;
    std::vector<std::vector<long long>> length_;
    std::vector<long long> weight_;
    std::vector<int> cover_count_;
    std::vector<int> committed_count_;
    std::vector<char> committed_;
    std::vector<char> blocked_;
    std::unordered_map<MemoKey, Frame, MemoHash> memo_;

    long long budget_ = 0;
    long long calls_ = 0;
    long long reward_units_ = 0;
    bool reward_scaling_exact_ = false;

    bool deadline_set_ = false, timed_out_ = false;
    unsigned clock_probe_ = 0;
    std::chrono::steady_clock::time_point started_;
};


RecursiveGreedyResult recursiveGreedyOrienteering(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start, int end, double budget, double padding,
    const RecursiveGreedyParams& params) {

    RecursiveGreedyResult result;
    const int node_count = static_cast<int>(costs.size());

    if (node_count < 2)
        throw std::invalid_argument("recursive greedy: need at least two nodes.");
    if (start < 0 || start >= node_count || end < 0 || end >= node_count ||
        start == end)
        throw std::invalid_argument("recursive greedy: invalid start/end nodes.");
    if (static_cast<int>(member_pixels.size()) != node_count)
        throw std::invalid_argument(
            "recursive greedy: member_pixels must have one list per node.");
    if (!std::isfinite(budget) || budget < 0.0)
        throw std::invalid_argument("recursive greedy: budget must be finite and >= 0.");
    if (!std::isfinite(padding) || padding < 0.0 || padding > budget)
        throw std::invalid_argument(
            "recursive greedy: padding must be finite and in [0, budget].");
    if (params.max_depth < 1 || params.max_depth > 30 ||
        params.reward_levels < 1 ||
        !std::isfinite(params.time_quantum) || params.time_quantum <= 0.0 ||
        !std::isfinite(params.time_limit) || params.time_limit < 0.0)
        throw std::invalid_argument("recursive greedy: invalid parameters.");

    for (double probability : pixel_probs)
        if (!std::isfinite(probability) || probability < 0.0)
            throw std::invalid_argument(
                "recursive greedy: pixel probabilities must be finite and nonnegative.");

    for (int i = 0; i < node_count; ++i) {
        if (static_cast<int>(costs[i].size()) != node_count)
            throw std::invalid_argument("recursive greedy: cost matrix must be square.");
        for (int j = 0; j < node_count; ++j)
            if (!std::isfinite(costs[i][j]) || costs[i][j] < 0.0)
                throw std::invalid_argument(
                    "recursive greedy: costs must be finite and nonnegative.");
        for (int p : member_pixels[i])
            if (p < 0 || p >= static_cast<int>(pixel_probs.size()))
                throw std::out_of_range("recursive greedy: invalid compact pixel id.");
    }
    if (!member_pixels[start].empty() || !member_pixels[end].empty())
        throw std::invalid_argument(
            "recursive greedy: dummy start/end pixel lists must be empty.");

    // Use the same union-coverage data structure as the greedy,
    // metaheuristic, GCP, and reporting code. The recursive solver keeps only
    // an additional integer-valued residual state required by RG-QP's reward
    // thresholds and rollback operations.
    CoverageOracle coverage_oracle(pixel_probs, member_pixels);
    if (coverage_oracle.nodeCount() != node_count)
        throw std::invalid_argument(
            "recursive greedy: coverage oracle and graph sizes differ.");

    RecursiveGreedySolver solver(
        costs, pixel_probs, member_pixels, coverage_oracle,
        start, end, budget, params);

    result.metric_ok = solver.verifyMetric();
    result.guaranteed_arcs = 1 << std::max(0, params.max_depth - 1);
    result.ratio_bound = 1.0 / (1.0 + std::max(0, params.max_depth - 1));

    const int first_depth = params.iterative_deepening ? 1 : params.max_depth;
    double best_reward = -1.0;
    for (int depth = first_depth; depth <= params.max_depth; ++depth) {
        std::vector<int> candidate = solver.solve(depth);
        const bool complete = !solver.timedOut();
        const bool usable = candidate.size() >= 2 &&
                            candidate.front() == start && candidate.back() == end;
        if (usable) {
            const double reward = solver.trueReward(candidate);
            if (reward > best_reward) {
                best_reward = reward;
                result.path = std::move(candidate);
            }
        }
        if (!complete) break;
        // Only a pass that both finished and produced a usable route may
        // claim a depth; the retained reward is monotone in depth, so the
        // bound for this depth applies to whatever route we kept.
        if (usable) result.completed_depth = depth;
    }
    result.calls = solver.calls();
    result.timed_out = solver.timedOut();
    result.reward_units = solver.rewardUnits();

    if (result.path.size() < 2 || result.path.front() != start ||
        result.path.back() != end) {
        result.path.clear();
        if (params.verbose)
            std::cerr << "recursive greedy: no feasible s-t path within budget.\n";
        return result;
    }

    std::vector<char> seen(node_count, 0);
    for (int v : result.path) {
        if (seen[v])
            throw std::runtime_error("recursive greedy: produced a non-simple path.");
        seen[v] = 1;
    }

    result.reward = coverage_oracle.rewardOf(result.path);

    result.cost = solver.realCost(result.path);
    const double cost_tolerance =
        1e-9 * std::max({1.0, result.cost, budget});
    if (result.cost > budget + cost_tolerance)
        throw std::runtime_error(
            "recursive greedy: produced an over-budget path.");
    if (result.cost + cost_tolerance < padding)
        throw std::runtime_error(
            "recursive greedy: path cost is smaller than endpoint padding.");

    result.true_duration = result.cost - padding;
    result.guaranteed_arcs = 1 << std::max(0, result.completed_depth - 1);
    result.ratio_bound = 1.0 / (1.0 + std::max(0, result.completed_depth - 1));
    result.guarantee_holds =
        result.metric_ok && params.candidate_limit <= 0 &&
        result.completed_depth > 0 && solver.rewardScalingExact();

    if (params.verbose)
        std::cout << "recursive greedy: reward=" << result.reward
                  << ", duration=" << result.true_duration
                  << " (budget " << budget - padding << ")"
                  << ", tiles=" << result.path.size() - 2
                  << ", frames=" << result.calls
                  << ", depth=" << result.completed_depth
                  << ", guarantee=" << (result.guarantee_holds ? "yes" : "no")
                  << " (1/" << 1.0 / result.ratio_bound << " up to "
                  << result.guaranteed_arcs << " arcs)\n";
    return result;
}