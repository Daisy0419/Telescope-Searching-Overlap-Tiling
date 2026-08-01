
#include "greedy_mink.h"
#include "ILP_gurobi_mink.h"
#include "ReadData.h"
#include "helplers.h"
#include "BuildGraph.h"
#include "GCP.h"
#include "SimulatedAnnealing.h"

#include <chrono>
#include <fstream>
#include <filesystem> 
#include <numeric> 
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_set>
#include <limits>
#include <algorithm>
#include <optional>
#include <random>

std::vector<std::vector<int>>
convertPathsToRanks(const std::vector<std::vector<int>>& paths,
                    const std::vector<int>& ranks) {
    std::vector<std::vector<int>> paths_in_rank;
    paths_in_rank.reserve(paths.size());

    for (const auto& path : paths) {
        std::vector<int> p;
        p.reserve(path.size());

        for (int node : path) {
            if (node < 0 || node >= static_cast<int>(ranks.size())) {
                std::cerr << "Invalid node index in path: " << node << "\n";
                continue;
            }
            p.push_back(ranks[node]);
        }

        paths_in_rank.push_back(std::move(p));
    }

    return paths_in_rank;
}

void run_mink_unrooted(const std::string& tilefile,
                       const std::string& mapfile,
                       const std::string& out_file,
                       double budget, double dwell_zenith,
                       double w_max, double w_acc, double settle_time,
                       bool is_deepslow, int m,
                       double mipGap, double timeLimit) {

    std::vector<std::vector<double>> slew_costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    // auto tiling_start = std::chrono::high_resolution_clock::now();
    // int num_tiles;
    // m = 0;
    // auto start_indices = buildGraphMultiStarts(
    //     tilefile, mapfile, slew_costs, probability, ranks, dwell_times,
    //     dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    // auto tiling_end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;
    auto tiling_start = std::chrono::high_resolution_clock::now();

    int num_tiles;
    m = 0;

    auto start_indices = buildGraphMultiStarts(
        tilefile, mapfile, slew_costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    // Debug check: can each tile be covered by a single unrooted route?
    double max_dwell = 0.0;
    int worst = -1;

    for (int i = 0; i < num_tiles; ++i) {
        if (dwell_times[i] > max_dwell) {
            max_dwell = dwell_times[i];
            worst = i;
        }
    }

    std::cout << "[unrooted debug] num_tiles = " << num_tiles << "\n";
    std::cout << "[unrooted debug] max tile dwell = " << max_dwell
            << " at local tile node " << worst
            << ", rank = " << ranks[worst]
            << ", budget = " << budget << "\n";

    if (max_dwell > budget) {
        std::cout << "[unrooted debug] WARNING: at least one tile has dwell time > budget, "
                << "so full cover is infeasible for unrooted min-k.\n";
    }

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::vector<std::vector<int>> paths;
    // --- Greedy ---
    {
        auto wall_start = std::chrono::high_resolution_clock::now();

        std::cout << "*************run greedy min routes****************\n";
        int k = UnrootedMinKCover(slew_costs, dwell_times, budget, paths);

        auto wall_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = wall_end - wall_start;
        std::cout << "running time (wallclock): " << elapsed.count()
                  << " seconds\n";

        int covered_tiles = 0;
        for (const auto& path : paths) {
            covered_tiles += path.size();
        }
        bool isValid = (covered_tiles == num_tiles);
        auto paths_in_rank = convertPathsToRanks(paths, ranks);
        save_result_without_paths(out_file, "mink_unrooted_greedy", mapfile,
                                  budget, w_max, w_acc, num_tiles, covered_tiles, k, k,
                                  tiling_elapsed.count(), elapsed.count(), paths_in_rank, isValid);
    }

    // --- ILP ---
    {
        auto wall_start = std::chrono::high_resolution_clock::now();

        std::cout << "*************run ILP min routes****************\n";
        auto [ilp_paths, best_k] = MinKCoverUnrooted_ILP(
            budget, slew_costs, dwell_times, mipGap, timeLimit, paths);

        auto wall_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = wall_end - wall_start;
        std::cout << "running time (wallclock): " << elapsed.count()
                  << " seconds\n";

        int covered_tiles = 0;
        for (const auto& path : ilp_paths) {
            covered_tiles += path.size();
        }
        bool isValid = (covered_tiles == num_tiles);
        auto ilp_paths_in_rank = convertPathsToRanks(ilp_paths, ranks);
        save_result_without_paths(out_file, "mink_unrooted_ilp", mapfile,
                                  budget, w_max, w_acc, num_tiles, covered_tiles,
                                  ilp_paths.size(), best_k,
                                  tiling_elapsed.count(), elapsed.count(), ilp_paths_in_rank, isValid);
    }
}


std::vector<std::pair<double, double>> run_mink_rooted(const std::string& tilefile,
                     const std::string& mapfile,
                     const std::string& out_file,
                     double budget, double dwell_zenith,
                     double w_max, double w_acc, double settle_time,
                     bool is_deepslow, int m,
                     double mipGap, double timeLimit) {

    std::vector<std::vector<double>> slew_costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();

    int num_tiles;
    auto start_indices = buildGraphMultiStarts(
        tilefile, mapfile, slew_costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    std::cout << "[rooted debug] num_tiles = " << num_tiles << "\n";
    std::cout << "[rooted debug] num_starts = " << start_indices.size() << "\n";

    int unreachable_count = 0;

    for (int v = 0; v < num_tiles; ++v) {
        double best = std::numeric_limits<double>::infinity();
        int best_start = -1;

        for (int s : start_indices) {
            double c = slew_costs[s][v] + dwell_times[v];
            if (c < best) {
                best = c;
                best_start = s;
            }
        }

        if (best > budget) {
            ++unreachable_count;
            std::cout << "[rooted debug] tile local node " << v
                    << ", rank = " << ranks[v]
                    << " cannot be covered by any start. best cost = "
                    << best
                    << ", best start node = " << best_start
                    << ", best start rank = " << ranks[best_start]
                    << ", budget = " << budget << "\n";
        }
    }

    if (unreachable_count > 0) {
        std::cout << "[rooted debug] WARNING: " << unreachable_count
                << " tiles are individually unreachable under this budget.\n";
    }

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    // --- Greedy ---
    std::vector<std::vector<int>> paths;
    {
        auto wall_start = std::chrono::high_resolution_clock::now();

        std::cout << "*************run greedy min routes****************\n";
        
        int k = RootedOpenMinKCover(slew_costs, dwell_times, start_indices,
                                    budget, paths);

        auto wall_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = wall_end - wall_start;
        std::cout << "running time (wallclock): " << elapsed.count()
                  << " seconds\n";
        
        int covered_tiles = 0;
        for (const auto& path : paths) {
            covered_tiles += (path.size() - 1);
        }
        bool isValid = (covered_tiles == num_tiles);
        auto paths_in_rank = convertPathsToRanks(paths, ranks);
        save_result_without_paths(out_file, "mink_rooted_greedy", mapfile,
                                  budget, w_max, w_acc, num_tiles, covered_tiles, k, k,
                                  tiling_elapsed.count(), elapsed.count(), paths_in_rank, isValid);
    }

    // --- ILP ---
    {
        auto wall_start = std::chrono::high_resolution_clock::now();

        std::cout << "*************run ILP min routes****************\n";
        auto [ilp_paths, best_k] = MinKCoverRootedOpen_ILP(
            start_indices, budget, slew_costs, dwell_times,
            mipGap, timeLimit, paths);

        auto wall_end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = wall_end - wall_start;
        std::cout << "running time (wallclock): " << elapsed.count()
                  << " seconds\n";

        int covered_tiles = 0;
        for (const auto& path : ilp_paths) {
            covered_tiles += (path.size() - 1);
        }
        bool isValid = (covered_tiles == num_tiles);
        auto ilp_paths_in_rank = convertPathsToRanks(ilp_paths, ranks);
        save_result_without_paths(out_file, "mink_rooted_ILP", mapfile,
                                  budget, w_max, w_acc, num_tiles, covered_tiles,
                                  ilp_paths.size(), best_k,
                                  tiling_elapsed.count(), elapsed.count(), ilp_paths_in_rank, isValid);
    }

    std::vector<std::pair<double, double>> init_poses(start_indices.size());
    for (std::size_t i = 0; i < start_indices.size(); ++i) {
        init_poses[i] = {ra[start_indices[i]], dec[start_indices[i]]};
    }
    return init_poses;
}


void mink_annealing_rooted(const std::string& tilefile,
                           const std::string& mapfile,
                           const std::string& out_file,
                           double budget, double dwell_zenith,
                           double w_max, double w_acc, double settle_time,
                           bool is_deepslow, int m,
                           std::vector<std::pair<double, double>> init_poses) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();

    int num_tiles = 0;

    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time,
        is_deepslow, m, num_tiles, std::nullopt, init_poses);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::cout << "padding: " << padding << "\n";

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> all_paths;

    // Real tile nodes are exactly 0, ..., num_tiles - 1.
    std::unordered_set<int> covered_tiles;

    // Starts that have not been used yet.
    std::vector<int> remaining_starts = start_indices;

    const int N = static_cast<int>(costs.size());
    const double INF = 1e12;

    auto is_real_tile = [&](int node) -> bool {
        return 0 <= node && node < num_tiles;
    };

    auto count_new_tiles = [&](const std::vector<int>& path) -> int {
        int count = 0;
        for (int node : path) {
            if (is_real_tile(node) &&
                covered_tiles.find(node) == covered_tiles.end()) {
                ++count;
            }
        }
        return count;
    };

    while (static_cast<int>(covered_tiles.size()) < num_tiles &&
           !remaining_starts.empty()) {

        std::cout << "remaining real tiles: "
                  << (num_tiles - static_cast<int>(covered_tiles.size()))
                  << "    remaining starts: " << remaining_starts.size()
                  << "    budget: " << budget << "\n";

        // Build the forbidden set once per outer iteration:
        //   - all start nodes (we selectively unblock one at a time)
        //   - all already-covered tiles
        // Store original costs so we can restore per-candidate.
        std::unordered_set<int> all_starts_set(remaining_starts.begin(),
                                               remaining_starts.end());

        std::unordered_set<int> forbidden;
        for (int s : remaining_starts)
            forbidden.insert(s);
        for (int v : covered_tiles)
            forbidden.insert(v);

        // Apply forbidden costs once on a shared local matrix
        std::vector<std::vector<double>> base_costs = costs;
        for (int v : forbidden) {
            for (int i = 0; i < N; ++i) {
                if (i == v) continue;
                base_costs[v][i] = INF;
                base_costs[i][v] = INF;
            }
        }

        // Prize vector: uncovered real tiles get 1.0, everything else 0.0
        std::vector<double> cover_prize(probability.size(), 0.0);
        for (int v = 0; v < num_tiles; ++v) {
            if (covered_tiles.find(v) == covered_tiles.end()) {
                cover_prize[v] = 1.0;
            }
        }

        int best_start_pos = -1;
        std::vector<int> best_path;
        int best_new_count = 0;

        for (int pos = 0; pos < static_cast<int>(remaining_starts.size()); ++pos) {
            int current_start = remaining_starts[pos];

            // Unblock this start: restore its original row/column
            std::vector<std::vector<double>> local_costs = base_costs;
            for (int i = 0; i < N; ++i) {
                local_costs[current_start][i] = costs[current_start][i];
                local_costs[i][current_start] = costs[i][current_start];
            }
            // Re-block edges between this start and other forbidden nodes
            // (other starts and covered tiles should still be unreachable)
            for (int v : forbidden) {
                if (v == current_start) continue;
                local_costs[current_start][v] = INF;
                local_costs[v][current_start] = INF;
            }

            std::vector<int> path = simulated_annealing_optimization(
                local_costs, cover_prize, budget + padding,
                current_start, end_idx);

            int new_count = count_new_tiles(path);

            if (new_count > best_new_count) {
                best_new_count = new_count;
                best_start_pos = pos;
                best_path = std::move(path);
            }
        }

        // No remaining start can cover a new real tile.
        if (best_start_pos < 0 || best_new_count == 0) {
            std::cout << "[rooted SA] No start can cover any new tile. Stop.\n";
            break;
        }

        int chosen_start = remaining_starts[best_start_pos];

        // Output path: start rank + newly-covered real tile ranks (no dummies)
        std::vector<int> path_in_rank;
        path_in_rank.push_back(ranks[chosen_start]);

        for (int node : best_path) {
            if (!is_real_tile(node)) continue;
            if (covered_tiles.count(node)) continue;

            covered_tiles.insert(node);
            path_in_rank.push_back(ranks[node]);
        }

        all_paths.push_back(std::move(path_in_rank));

        // Consume this start
        remaining_starts.erase(remaining_starts.begin() + best_start_pos);
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;

    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    int covered_count = static_cast<int>(covered_tiles.size());
    int k = static_cast<int>(all_paths.size());
    bool isValid = (covered_count == num_tiles);

    save_result_without_paths(out_file, "mink_rooted_annealing", mapfile,
                              budget, w_max, w_acc,
                              num_tiles, covered_count, k, k,
                              tiling_elapsed.count(), elapsed.count(),
                              all_paths, isValid);
}

void mink_annealing_rooted_random(const std::string& tilefile,
                           const std::string& mapfile,
                           const std::string& out_file,
                           double budget, double dwell_zenith,
                           double w_max, double w_acc, double settle_time,
                           bool is_deepslow, int m,
                           std::vector<std::pair<double, double>> init_poses) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();

    int num_tiles = 0;

    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time,
        is_deepslow, m, num_tiles, std::nullopt, init_poses);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::cout << "padding: " << padding << "\n";

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> all_paths;

    std::unordered_set<int> covered_tiles;
    std::vector<int> remaining_starts = start_indices;

    const int N = static_cast<int>(costs.size());
    const double INF = 1e12;

    std::mt19937 rng(419);

    auto is_real_tile = [&](int node) -> bool {
        return 0 <= node && node < num_tiles;
    };

    auto count_new_tiles = [&](const std::vector<int>& path) -> int {
        int count = 0;
        for (int node : path) {
            if (is_real_tile(node) &&
                covered_tiles.find(node) == covered_tiles.end()) {
                ++count;
            }
        }
        return count;
    };

    while (static_cast<int>(covered_tiles.size()) < num_tiles &&
           !remaining_starts.empty()) {

        std::cout << "remaining real tiles: "
                  << (num_tiles - static_cast<int>(covered_tiles.size()))
                  << "    remaining starts: " << remaining_starts.size()
                  << "    budget: " << budget << "\n";

        // Forbidden set: all unused starts + all covered tiles.
        std::unordered_set<int> forbidden;
        for (int s : remaining_starts)
            forbidden.insert(s);
        for (int v : covered_tiles)
            forbidden.insert(v);

        std::vector<std::vector<double>> base_costs = costs;
        for (int v : forbidden) {
            for (int i = 0; i < N; ++i) {
                if (i == v) continue;
                base_costs[v][i] = INF;
                base_costs[i][v] = INF;
            }
        }

        // Prize vector: uncovered real tiles get 1.0, everything else 0.0
        std::vector<double> cover_prize(probability.size(), 0.0);
        for (int v = 0; v < num_tiles; ++v) {
            if (covered_tiles.find(v) == covered_tiles.end()) {
                cover_prize[v] = 1.0;
            }
        }

        // ----- pick ONE random root for this step -----
        std::uniform_int_distribution<int> pick(
            0, static_cast<int>(remaining_starts.size()) - 1);
        const int pos = pick(rng);
        const int current_start = remaining_starts[pos];

        // Unblock the chosen start, keep all other forbidden nodes blocked.
        std::vector<std::vector<double>> local_costs = base_costs;
        for (int i = 0; i < N; ++i) {
            local_costs[current_start][i] = costs[current_start][i];
            local_costs[i][current_start] = costs[i][current_start];
        }
        for (int v : forbidden) {
            if (v == current_start) continue;
            local_costs[current_start][v] = INF;
            local_costs[v][current_start] = INF;
        }

        std::vector<int> path = simulated_annealing_optimization(
            local_costs, cover_prize, budget + padding,
            current_start, end_idx);

        const int new_count = count_new_tiles(path);

        // Consume this start regardless of outcome (guarantees termination).
        remaining_starts.erase(remaining_starts.begin() + pos);

        if (new_count == 0) {
            std::cout << "[rooted SA] random start " << current_start
                      << " covered no new tile; discarding it.\n";
            continue;
        }

        // Commit: start rank + newly-covered real tile ranks.
        std::vector<int> path_in_rank;
        path_in_rank.push_back(ranks[current_start]);

        for (int node : path) {
            if (!is_real_tile(node)) continue;
            if (covered_tiles.count(node)) continue;

            covered_tiles.insert(node);
            path_in_rank.push_back(ranks[node]);
        }

        all_paths.push_back(std::move(path_in_rank));
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;

    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    int covered_count = static_cast<int>(covered_tiles.size());
    int k = static_cast<int>(all_paths.size());
    bool isValid = (covered_count == num_tiles);

    save_result_without_paths(out_file, "mink_rooted_annealing", mapfile,
                              budget, w_max, w_acc,
                              num_tiles, covered_count, k, k,
                              tiling_elapsed.count(), elapsed.count(),
                              all_paths, isValid);
}




void mink_annealing_rooted_best(const std::string& tilefile,
                                const std::string& mapfile,
                                const std::string& out_file,
                                double budget, double dwell_zenith,
                                double w_max, double w_acc, double settle_time,
                                bool is_deepslow, int m,
                                std::vector<std::pair<double, double>> init_poses) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();

    int num_tiles = 0;

    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time,
        is_deepslow, m, num_tiles, std::nullopt, init_poses);


    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::cout << "padding: " << padding << "\n";
    std::cout << "[rooted SA best debug] num_tiles = " << num_tiles << "\n";
    std::cout << "[rooted SA best debug] num_starts = " << start_indices.size() << "\n";
    std::cout << "[rooted SA best debug] end_idx = " << end_idx << "\n";
    std::cout << "[rooted SA best debug] graph size = " << probability.size() << "\n";

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> all_paths;

    // Real tile nodes are exactly 0, ..., num_tiles - 1.
    std::unordered_set<int> covered_tiles;

    // Each start can be used at most once.
    std::vector<int> remaining_starts = start_indices;

    const double INF = 1e12;

    auto is_real_tile = [&](int node) -> bool {
        return 0 <= node && node < num_tiles;
    };

    auto count_new_tiles = [&](const std::vector<int>& path) -> int {
        std::unordered_set<int> new_tiles;

        for (int node : path) {
            if (is_real_tile(node) &&
                covered_tiles.find(node) == covered_tiles.end()) {
                new_tiles.insert(node);
            }
        }

        return static_cast<int>(new_tiles.size());
    };

    auto new_probability_sum = [&](const std::vector<int>& path) -> double {
        std::unordered_set<int> new_tiles;
        double total = 0.0;

        for (int node : path) {
            if (is_real_tile(node) &&
                covered_tiles.find(node) == covered_tiles.end() &&
                new_tiles.find(node) == new_tiles.end()) {
                new_tiles.insert(node);
                total += probability[node];
            }
        }

        return total;
    };

    while (static_cast<int>(covered_tiles.size()) < num_tiles &&
           !remaining_starts.empty()) {

        std::cout << "remaining real tiles: "
                  << (num_tiles - static_cast<int>(covered_tiles.size()))
                  << "    remaining starts: " << remaining_starts.size()
                  << "    budget: " << budget << "\n";

        int best_start_pos = -1;
        std::vector<int> best_path;

        int best_new_count = 0;
        double best_new_prob = -1.0;

        // Try every unused start and choose the path that covers
        // the largest number of new real tiles.
        for (int pos = 0; pos < static_cast<int>(remaining_starts.size()); ++pos) {
            int current_start = remaining_starts[pos];

            // Cover objective:
            // uncovered real tiles have prize 1;
            // covered tiles, starts, and end node have prize 0.
            std::vector<double> cover_prize(probability.size(), 0.0);

            for (int v = 0; v < num_tiles; ++v) {
                if (covered_tiles.find(v) == covered_tiles.end()) {
                    cover_prize[v] = 1.0;
                }
            }

            // Make a local copy of the graph for this candidate start.
            // We forbid all other starts and all already-covered tiles.
            std::vector<std::vector<double>> local_costs = costs;

            auto forbid_node = [&](int v) {
                for (int i = 0; i < static_cast<int>(local_costs.size()); ++i) {
                    if (i == v) continue;
                    local_costs[v][i] = INF;
                    local_costs[i][v] = INF;
                }
            };

            // Forbid all other starts.
            for (int s : remaining_starts) {
                if (s != current_start) {
                    forbid_node(s);
                }
            }

            // Forbid already-covered real tiles.
            for (int v : covered_tiles) {
                forbid_node(v);
            }

            std::vector<int> path = simulated_annealing_optimization(
                local_costs, cover_prize,
                budget + padding,
                current_start, end_idx);

            int new_count = count_new_tiles(path);
            double new_prob = new_probability_sum(path);

            // Primary objective: cover more new tiles.
            // Tie-breaker: cover larger original probability mass.
            if (new_count > best_new_count ||
                (new_count == best_new_count && new_prob > best_new_prob)) {
                best_new_count = new_count;
                best_new_prob = new_prob;
                best_start_pos = pos;
                best_path = std::move(path);
            }
        }

        if (best_start_pos < 0 || best_new_count == 0) {
            std::cout << "[rooted SA best debug] No remaining start can cover a new tile. Stop.\n";
            break;
        }

        int chosen_start = remaining_starts[best_start_pos];

        // Convert to output ranks.
        // Important: only output the chosen start plus newly-covered real tiles.
        // Do not output other starts, dummy end, or already-covered tiles.
        std::vector<int> path_in_rank;
        path_in_rank.push_back(ranks[chosen_start]);

        std::unordered_set<int> new_tiles_this_path;

        for (int node : best_path) {
            if (!is_real_tile(node)) continue;
            if (covered_tiles.find(node) != covered_tiles.end()) continue;
            if (new_tiles_this_path.find(node) != new_tiles_this_path.end()) continue;

            new_tiles_this_path.insert(node);
            path_in_rank.push_back(ranks[node]);
        }

        for (int node : new_tiles_this_path) {
            covered_tiles.insert(node);
        }

        all_paths.push_back(std::move(path_in_rank));

        std::cout << "[rooted SA best debug] selected start rank = "
                  << ranks[chosen_start]
                  << ", new tiles covered = " << new_tiles_this_path.size()
                  << ", total covered = " << covered_tiles.size()
                  << " / " << num_tiles << "\n";

        // This start is now used.
        remaining_starts.erase(remaining_starts.begin() + best_start_pos);
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;

    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    int covered_count = static_cast<int>(covered_tiles.size());
    int k = static_cast<int>(all_paths.size());
    bool isValid = (covered_count == num_tiles);

    std::cout << "[rooted SA best debug] covered real tiles = "
              << covered_count << " / " << num_tiles << "\n";

    save_result_without_paths(out_file, "mink_rooted_annealing_best", mapfile,
                              budget, w_max, w_acc,
                              num_tiles, covered_count, k, k,
                              tiling_elapsed.count(), elapsed.count(),
                              all_paths, isValid);
}



void mink_annealing_unrooted(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow) {

    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;

    auto tiling_start = std::chrono::high_resolution_clock::now();
    int num_tiles;

    auto [start_idx, end_idx, padding] = buildGraphOrienteeringDummyRoot(
        tilefile, mapfile, costs, probability, ranks, dwell_times,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, num_tiles);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::cout << "padding: " << padding << "\n";

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> all_paths;
    int k = 0;
    while (probability.size() > 2) {  // start + end dummies remain
        std::cout << "remaining tiles: " << costs.size() << "    ";
        std::cout << "budget: " << budget << "\n";

        std::vector<int> path = simulated_annealing_optimization(costs, probability,
                                    budget + padding,
                                    start_idx, end_idx);

        // Build the list of nodes to remove (visited interior nodes)
        std::vector<int> nodes_to_remove;
        for (int node : path) {
            if (node != start_idx && node != end_idx)
                nodes_to_remove.push_back(node);
        }

        // No tiles visited — SA couldn't make progress, stop
        if (nodes_to_remove.empty()) break;

        // Convert path to ranks for output
        std::vector<int> path_in_rank;
        path_in_rank.reserve(path.size());
        for (int node : path)
            path_in_rank.push_back(ranks[node]);
        all_paths.push_back(std::move(path_in_rank));

        k += 1;

        std::vector<int> start_indices = {start_idx};
        auto [new_starts, new_end] = removeVisitedNodes(
                probability, costs, ranks, nodes_to_remove,
                start_indices, end_idx);

        start_idx = new_starts[0];
        end_idx = new_end;
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;
    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    int covered_tiles = 0;
    for (const auto& path : all_paths) {
        covered_tiles += static_cast<int>(path.size()) - 2;  // subtract start and end dummies
    }

    save_result_without_paths(out_file, "mink_unrooted_annealing", mapfile,
                            budget, w_max, w_acc, num_tiles, covered_tiles, k, k,
                            tiling_elapsed.count(), elapsed.count(), all_paths,
                            covered_tiles == num_tiles);
}


