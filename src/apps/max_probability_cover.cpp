#include "max_probability_cover.h"
#include "AntColony.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "ILP_gurobi_top.h"
#include "GCP.h"
#include "ILP_cplex_top.h"
#include "PathSplitter.h"
#include "ParticleSwarm.h"
#include "ParticleSwarmInspired.h"
#include "SimulatedAnnealing.h"
#include "PathRelink.h"
#include "MMkTC.h"
#include "helplers.h"
#include "PathOperations.h"
#include "PathCover.h"
#include "AnnealingMulti.h"
#include "BuildGraph.h"
#include "ILP_gurobi_op.h"

#include <chrono>
#include <fstream>
#include <filesystem> 
#include <numeric> 


void run_top(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m,
                            double accu_thr,
                            double time_limit) {
                                
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;


    // double time_limit = 2400;
    // double accu_thr = 0.01;
 
    int num_tiles;
 
    auto tiling_start = std::chrono::high_resolution_clock::now();
    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);
 
    std::cout << "padding: " << padding << "\n";

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;
 
    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;
    bool is_valid = false;
 
    std::vector<std::vector<int>> all_paths;
    double total_prize = 0.0;
 
    // team annealing
    std::cout << "*********running top annealing*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    // genetic_optimization(costs, budget, 0);
    auto [paths_top_annealing, prize_top_annealing] = TopAnnealing(costs, probability, budget + padding, m, start_indices, end_idx);
    end = std::chrono::high_resolution_clock::now();
    is_valid = verify_paths(costs, paths_top_annealing, start_indices, end_idx, costs.size(), budget + padding);
    // print_path(costs, probability, ranks, pso_path, padding);
    
    elapsed_seconds = end - start;
    std::cout << "total prize = " << prize_top_annealing << "\n";
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result_with_paths(out_file, "top_annealing", mapfile, paths_top_annealing, ranks, budget, w_max, w_acc, m, num_tiles, prize_top_annealing, prize_top_annealing, tiling_elapsed.count(), elapsed_seconds.count(), is_valid);


    // std::cout << "*********running Path Cover*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // auto [paths_pc, prize_pc] = PathCoverWraper(costs, probability, budget + padding, start_indices, end_idx, m);
    // end = std::chrono::high_resolution_clock::now();
    // is_valid = verify_paths_any_order(costs, paths_pc, start_indices, end_idx, costs.size(), budget + padding);
    // elapsed_seconds = end - start;
    // double pc = multiPathsSumPrize(probability, paths_pc);
    // if(std::abs(pc-prize_pc)>10e9) {
    //     std::cout << "wrong prize get\n";
    //     prize_pc = pc;
    // }
    // std::cout << "total prize = " << prize_pc << "\n";
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;

    // save_result_with_paths(out_file, "top_PathCover", mapfile, paths_pc, ranks, budget, w_max, w_acc, m, num_tiles, prize_pc, prize_pc, tiling_elapsed.count(), elapsed_seconds.count(), is_valid);


    std::cout << "*********running gurobi*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();
    auto [routes_gurobi, total_prize_gurobi, prize_bound] = gurobiSolveTeamSTBound(costs, probability, start_indices, end_idx, 
                                                                budget + padding, m, accu_thr, time_limit, {});
    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;
    is_valid = verify_paths(costs, routes_gurobi, start_indices, end_idx, costs.size(), budget + padding);
    
    std::cout << "gurobi best total prize = " << total_prize_gurobi << "\n";
    std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    save_result_with_paths(out_file, "top_ILP", mapfile, routes_gurobi, ranks, budget, w_max, w_acc, m, num_tiles, total_prize_gurobi, prize_bound, tiling_elapsed.count(), elapsed_seconds.count(), is_valid);


    // // path relinking
    // std::cout << "*********running path relinking*********" << std::endl;
    // start = std::chrono::high_resolution_clock::now();
    // auto [paths_pr, prize_pr] = PROptimization_parallel(costs, probability, budget + padding, m, start_indices, end_idx);
    // end = std::chrono::high_resolution_clock::now();
    // elapsed_seconds = end - start;
    // is_valid = verify_paths(costs, paths_pr, start_indices, end_idx, costs.size(), budget + padding);
    
    // std::cout << "total prize = " << prize_pr << "\n";
    // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
    // save_result_with_paths(out_file, "PathRelink", mapfile, paths_pr, ranks, budget, w_max, w_acc, m, prize_pr, elapsed_seconds.count(), is_valid);
                 
}


void run_greedy_multistarts_gcp(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();
    int num_tiles;

    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;

    std::cout << "padding: " << padding << "\n";

    auto wall_start = std::chrono::high_resolution_clock::now();

    std::vector<std::vector<int>> all_paths;
    double total_prize = 0.0;
    int nTiles = num_tiles;
    while (!start_indices.empty()) {
        std::cout << "remaining tiles: " << nTiles << "    ";
        std::cout << "budget: " << budget << "\n";

        // Run GCP from every remaining start
        std::vector<std::vector<int>> candidate_paths;
        candidate_paths.reserve(start_indices.size());

        for (int start_idx : start_indices) {
            std::vector<std::vector<double>> costs_single_s;
            auto [new_start, new_end] = extractSingleStartSubgraph(
                costs, start_indices, start_idx, end_idx, costs_single_s);

            std::vector<double> prob_single_s = extractSingleStartProbability(
                probability, start_indices);

            std::vector<int> path = GCP(costs_single_s, prob_single_s,
                                        budget + padding, new_start, new_end);

            // Remap start and end back to full-graph indices
            path.front() = start_idx;
            path.back() = end_idx;
            candidate_paths.push_back(std::move(path));
        }

        // Pick the best path by total prize
        int best_idx = -1;
        double best_prize = -1.0;

        for (int i = 0; i < static_cast<int>(candidate_paths.size()); ++i) {
            double p = pathPrize(probability, candidate_paths[i]);
            if (p > best_prize) {
                best_prize = p;
                best_idx   = i;
            }
        }

        if (best_idx < 0) break;
        total_prize += best_prize;

        const std::vector<int>& best_path = candidate_paths[best_idx];
        int best_start = start_indices[best_idx];

        print_path(costs, probability, ranks, best_path, padding);

        // Remove the chosen start from start_indices
        start_indices.erase(start_indices.begin() + best_idx);

        // Translate best_path to rank-space before the graph is rebuilt
        std::vector<int> path_in_rank;
        path_in_rank.reserve(best_path.size());
        for (int node : best_path)
            path_in_rank.push_back(ranks[node]);
        all_paths.push_back(std::move(path_in_rank));

        // Build the list of nodes to remove
        std::vector<int> nodes_to_remove;
        for (int node : best_path) {
            if (node != best_start && node != end_idx)
                nodes_to_remove.push_back(node);
        }
        nodes_to_remove.push_back(best_start);

        if (start_indices.empty()) break;

        auto [new_starts, new_end] = removeVisitedNodes(
                probability, costs, ranks, nodes_to_remove,
                start_indices, end_idx);

        start_indices = std::move(new_starts);
        end_idx = new_end;
        nTiles = (int)costs.size() - (int)start_indices.size() - 1;
    }

    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;
    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    save_result_with_paths2(out_file, "greedy+gcp", mapfile, all_paths,
                            budget, w_max, w_acc, m, num_tiles, total_prize, total_prize, tiling_elapsed.count(), elapsed.count());
}

void run_greedy_multistarts_annealing(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m) {
    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;
 
    auto tiling_start = std::chrono::high_resolution_clock::now();
    int num_tiles;

    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;
 
    std::cout << "padding: " << padding << "\n";
 
    auto wall_start = std::chrono::high_resolution_clock::now();
 
    std::vector<std::vector<int>> all_paths;
    double total_prize = 0.0;
 
    while (!start_indices.empty()) {
        std::cout << "remaining tiles: " << costs.size() << "    ";
        std::cout << "budget: " << budget << "\n";
        // Run GCP from every remaining start
        std::vector<std::vector<int>> candidate_paths;
        candidate_paths.reserve(start_indices.size());
 
        for (int start_idx : start_indices) {
            std::vector<int> path = simulated_annealing_optimization(costs, probability,
                                        budget + padding,
                                        start_idx, end_idx);
            candidate_paths.push_back(std::move(path));
        }
        std::cout << "budget: " << budget << "\n";

        // Pick the best path by total prize
        int best_idx = -1;
        double best_prize = -1.0;
 
        for (int i = 0; i < static_cast<int>(candidate_paths.size()); ++i) {
            double p = pathPrize(probability, candidate_paths[i]);
            if (p > best_prize) {
                best_prize = p;
                best_idx   = i;
            }
        }
 
        if (best_idx < 0) break;  // no valid path found
        total_prize += best_prize;
 
        const std::vector<int>& best_path = candidate_paths[best_idx];
        int best_start = start_indices[best_idx];

        print_path(costs, probability, ranks, best_path, padding);
 
        // Remove the chosen start from start_indices
        start_indices.erase(start_indices.begin() + best_idx);
 
        // Translate best_path to rank-space before the graph is rebuilt
        std::vector<int> path_in_rank;
        path_in_rank.reserve(best_path.size());
        for (int node : best_path)
            path_in_rank.push_back(ranks[node]);
        all_paths.push_back(std::move(path_in_rank));
 
        // Build the list of nodes to remove (visited interior nodes)
        // Exclude start and end from removal so other routes can still use them
        std::vector<int> nodes_to_remove;
        for (int node : best_path) {
            if (node != best_start && node != end_idx)
                nodes_to_remove.push_back(node);
        }
 
        // Also remove the used start 
        nodes_to_remove.push_back(best_start);
 
        // Shrink the graph and remap indices
        if (start_indices.empty()) break;
 
        auto [new_starts, new_end] = removeVisitedNodes(
                probability, costs, ranks, nodes_to_remove,
                start_indices, end_idx);

 
        start_indices = std::move(new_starts);
        end_idx= new_end;
    }
 
    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;
    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    save_result_with_paths2(out_file, "greedy+annealing", mapfile, all_paths, budget, w_max, w_acc, m, num_tiles, total_prize, total_prize, tiling_elapsed.count(), elapsed.count());
}


void run_greedy_multistarts_ilp(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m,
                            double mipGap, 
                            double timeLimit) {


    std::vector<std::vector<double>> costs;
    std::vector<double> probability;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<double> ra;
    std::vector<double> dec;

    auto tiling_start = std::chrono::high_resolution_clock::now();
    int num_tiles;
    // double mipGap = 0.001, timeLimit = 1800;
 
    auto [start_indices, end_idx, padding] = buildGraphTeamOrienteering(
        tilefile, mapfile, costs, probability, ranks, dwell_times, ra, dec,
        dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles);

    auto tiling_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> tiling_elapsed = tiling_end - tiling_start;
 
    std::cout << "padding: " << padding << "\n";
 
    auto wall_start = std::chrono::high_resolution_clock::now();
 
    std::vector<std::vector<int>> all_paths;
    double total_prize = 0.0;
 
    while (!start_indices.empty()) {
        std::cout << "remaining tiles: " << costs.size() << "    ";
        std::cout << "budget: " << budget << "\n";
        // Run GCP from every remaining start
        std::vector<std::vector<int>> candidate_paths;
        candidate_paths.reserve(start_indices.size());
 
        for (int start_idx : start_indices) {
            std::vector<int> path = gurobiSolveST_time(costs, probability,
                                        start_idx, end_idx, 
                                        budget + padding,
                                        mipGap, timeLimit);
            candidate_paths.push_back(std::move(path));
        }
 
        // Pick the best path by total prize
        int best_idx = -1;
        double best_prize = -1.0;
 
        for (int i = 0; i < static_cast<int>(candidate_paths.size()); ++i) {
            double p = pathPrize(probability, candidate_paths[i]);
            if (p > best_prize) {
                best_prize = p;
                best_idx   = i;
            }
        }
 
        if (best_idx < 0) break;  // no valid path found
        total_prize += best_prize;
 
        const std::vector<int>& best_path = candidate_paths[best_idx];
        int best_start = start_indices[best_idx];

        print_path(costs, probability, ranks, best_path, padding);
 
        // Remove the chosen start from start_indices
        start_indices.erase(start_indices.begin() + best_idx);
 
        // Translate best_path to rank-space before the graph is rebuilt
        std::vector<int> path_in_rank;
        path_in_rank.reserve(best_path.size());
        for (int node : best_path)
            path_in_rank.push_back(ranks[node]);
        all_paths.push_back(std::move(path_in_rank));
 
        // Build the list of nodes to remove (visited interior nodes)
        // Exclude start and end from removal so other routes can still use them
        std::vector<int> nodes_to_remove;
        for (int node : best_path) {
            if (node != best_start && node != end_idx)
                nodes_to_remove.push_back(node);
        }
 
        // Also remove the used start 
        nodes_to_remove.push_back(best_start);
 
        // Shrink the graph and remap indices
        if (start_indices.empty()) break;
 
        auto [new_starts, new_end] = removeVisitedNodes(
                probability, costs, ranks, nodes_to_remove,
                start_indices, end_idx);

 
        start_indices = std::move(new_starts);
        end_idx= new_end;
    }
 
    auto wall_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = wall_end - wall_start;
    std::cout << "running time (wallclock): " << elapsed.count()
              << " seconds\n";

    save_result_with_paths2(out_file, "greedy+ilp", mapfile, all_paths, budget, w_max, w_acc, m, num_tiles, total_prize, total_prize, tiling_elapsed.count(), elapsed.count());
}
