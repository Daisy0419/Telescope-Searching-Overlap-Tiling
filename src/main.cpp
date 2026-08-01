#include "BuildGraph.h"
#include "ILP_gurobi_path_cover_op.h"


#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>


double coveredProbabilityUnionSTPixelOverlap(
        const std::vector<int>& path_tiles,
        const std::vector<std::vector<int>>& member_pixels,
        const std::vector<double>& pixel_probs) {

    const int nTiles  = static_cast<int>(member_pixels.size());
    const int nPixels = static_cast<int>(pixel_probs.size());
    if (nPixels == 0) return 0.0;

    std::vector<char> seen(nPixels, 0);
    double total = 0.0;

    for (int tile : path_tiles) {
        if (tile < 0 || tile >= nTiles) {
            throw std::runtime_error(
                "coveredProbabilityUnionSTPixelOverlap: "
                "tile index out of range.");
        }
        for (int p : member_pixels[tile]) {
            if (p < 0 || p >= nPixels) {
                throw std::runtime_error(
                    "coveredProbabilityUnionSTPixelOverlap: "
                    "pixel index out of range.");
            }
            if (!seen[p]) {
                seen[p] = 1;
                total += pixel_probs[p];
            }
        }
    }

    return total;
}


// costs is the all-in edge matrix returned by buildGraph. Do not add dwell
// time again: each selected tile's dwell and the end padding are already on
// the corresponding edges.
double printPathSTPixelOverlap(
        const std::vector<std::vector<double>>& costs,
        const std::vector<std::vector<int>>& member_pixels,
        const std::vector<double>& pixel_probs,
        const std::vector<int>& rank,
        const std::vector<int>& best_path,
        int start_idx,
        int end_idx) {

    if (best_path.empty()) {
        std::cout << "No incumbent path was returned.\n";
        return 0.0;
    }
    if (best_path.front() != start_idx || best_path.back() != end_idx) {
        throw std::runtime_error(
            "printPathSTPixelOverlap: invalid path endpoints.");
    }

    double sum_probability = coveredProbabilityUnionSTPixelOverlap(
        best_path, member_pixels, pixel_probs);

    double total_cost = 0.0;
    for (std::size_t i = 1; i < best_path.size(); ++i) {
        total_cost += costs[best_path[i - 1]][best_path[i]];
    }

    std::cout << "Path(tile rank): ";
    for (std::size_t i = 1; i + 1 < best_path.size(); ++i) {
        const int tile = best_path[i];
        std::cout << (rank.empty() ? tile : rank.at(tile)) << " ";
    }
    std::cout << std::endl;

    const std::size_t num_selected_tiles =
        best_path.size() >= 2 ? best_path.size() - 2 : 0;
    std::cout << "Num Tiles in Path: " << num_selected_tiles << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}

void saveResultSTPixelOverlap(
        const std::string& filename,
        const std::string& method,
        const std::string& tiling_data,
        const std::string& likelihood_data,
        double budget,
        double slew_rate,
        const std::vector<std::vector<double>>& costs,
        const std::vector<std::vector<int>>& member_pixels,
        const std::vector<double>& pixel_probs,
        const std::vector<int>& rank,
        const std::pair<std::vector<std::vector<int>>, double>& result,
        int start_idx,
        int end_idx,
        double elapsed_time) {

    bool file_exists = std::filesystem::exists(filename);
    std::ofstream outfile(filename, std::ios::app);

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: "
                  << filename << std::endl;
        return;
    }

    if (!file_exists) {
        outfile
            << "Method,TilingData,LikelihoodData,Budget,SlewRate,NumTiles,"
               "SumProb,BestBound,TotalCost,TimeSec,Path\n";
    }

    if (result.first.size() > 1) {
        throw std::runtime_error(
            "saveResultSTPixelOverlap: the single-telescope solver "
            "returned more than one path.");
    }

    const std::vector<int> best_path =
        result.first.empty() ? std::vector<int>{} : result.first.front();
    double sum_probability = 0.0;
    double total_cost = 0.0;
    std::size_t num_selected_tiles = 0;

    if (!best_path.empty()) {
        if (best_path.front() != start_idx || best_path.back() != end_idx) {
            throw std::runtime_error(
                "saveResultSTPixelOverlap: invalid path endpoints.");
        }
        sum_probability = coveredProbabilityUnionSTPixelOverlap(
            best_path, member_pixels, pixel_probs);
        for (std::size_t i = 1; i < best_path.size(); ++i) {
            total_cost += costs[best_path[i - 1]][best_path[i]];
        }
        num_selected_tiles =
            best_path.size() >= 2 ? best_path.size() - 2 : 0;
    }

    outfile << method << ","
            << std::filesystem::path(tiling_data).filename().string() << ","
            << std::filesystem::path(likelihood_data).filename().string()
            << ","
            << budget << ","
            << slew_rate << ","
            << num_selected_tiles << ","
            << sum_probability << ","
            << result.second << ","
            << total_cost << ","
            << elapsed_time << ",";

    // Save only real tile ranks; omit the dummy start and end nodes.
    for (std::size_t i = 1; i + 1 < best_path.size(); ++i) {
        const int tile = best_path[i];
        outfile << (rank.empty() ? tile : rank.at(tile));
        if (i + 2 < best_path.size()) {
            outfile << " ";
        }
    }
    outfile << "\n";
}


void testPathsCoverSTPixelOverlapConfigured(
                                            std::string tilefile,
                                            std::string mapfile,
                                            std::string out_file,
                                            std::string overlap,
                                            double budget,
                                            double dwell_time,
                                            double w_max,
                                            double w_acc,
                                            double settle_time,
                                            bool is_deepslow,
                                            int init_pos_idx,
                                            double accu_thr,
                                            double time_limit) {

    std::vector<std::vector<double>> costs;
    std::vector<double> tile_probs;
    std::vector<int> ranks;
    std::vector<double> dwell_times;
    std::vector<std::vector<int>> member_pixels;
    std::vector<double> pixel_probs;

    int num_tiles = 0;

    const auto [start_idx, end_idx, padding] = buildGraphwithPixelsOP(
        tilefile, mapfile,
        costs, dwell_times, tile_probs, ranks,
        member_pixels, pixel_probs,
        dwell_time, w_max, w_acc, settle_time,
        is_deepslow, num_tiles, init_pos_idx);

    if (costs[start_idx][end_idx] > budget + 1e-9) {
        throw std::runtime_error(
            "testPathsCoverSTPixelOverlapConfigured: budget is smaller "
            "than the direct start-to-end padding cost.");
    }

    std::cout << "start: " << start_idx
              << ", end: " << end_idx
              << ", retained tiles: " << num_tiles
              << ", padding: " << padding << std::endl;

    std::cout
        << "*********ILP Path Cover (Continuous Time, Pixel Overlap)*********" << std::endl;

    const auto solve_start = std::chrono::high_resolution_clock::now();

    const std::pair<std::vector<std::vector<int>>, double> result =
        gurobiPathCoverPixelOverlapContinuousTime(costs, pixel_probs, member_pixels,
                                                start_idx, end_idx, budget,
                                                accu_thr, time_limit);


    // result = gurobiPathCoverPixelOverlapMTZ(
    //               costs, pixel_probs, member_pixels,
    //               start_idx, end_idx, budget,
    //               accu_thr, time_limit);

    const auto solve_end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double> elapsed_seconds = solve_end - solve_start;

    if (result.first.size() > 1) {
        throw std::runtime_error(
            "The single-telescope ILP returned more than one path.");
    }

    const std::vector<int> best_path = result.first.empty() ? std::vector<int>{} : result.first.front();

    std::cout << "Best Bound: " << result.second << std::endl;

    printPathSTPixelOverlap( costs, member_pixels, pixel_probs, ranks, best_path, start_idx, end_idx);

    std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

    saveResultSTPixelOverlap(out_file,
        overlap + "path+cover+continuous-time",
        tilefile, mapfile,
        budget, w_max,
        costs, member_pixels, pixel_probs,
        ranks, result,
        start_idx, end_idx,
        elapsed_seconds.count());
}


int main(int argc, char** argv) {
    const auto total_start =
        std::chrono::high_resolution_clock::now();

    std::string out_file = "../results/out.csv";
    std::string tilefile =
        "../data/tilings_small/6.9x6.9_tiling.csv";
    std::string mapfile =
        "../data/test_maps/GW200316_215756_50.txt";

    double budget = 60.0;
    double w_max = 10.0;
    double w_acc = 10.0;
    double dwell_zenith = 2.0;
    bool is_deepslow = false;

    // Fixed solver/graph parameters.
    const double settle_time = 0.0;
    const int source_tile = -2;
    const double mip_gap = 0.001;
    const double time_limit = 36000.0;

    if (argc > 9) {
        std::cerr
            << "Usage: " << argv[0]
            << " [MAP_FILE TILE_FILE BUDGET W_MAX W_ACC "
               "DWELL_ZENITH IS_DEEPSLOW OUTPUT_CSV]\n";
        return 2;
    }

    try {
        if (argc > 1) mapfile = argv[1];
        if (argc > 2) tilefile = argv[2];
        if (argc > 3) budget = std::stod(argv[3]);
        if (argc > 4) w_max = std::stod(argv[4]);
        if (argc > 5) w_acc = std::stod(argv[5]);
        if (argc > 6) dwell_zenith = std::stod(argv[6]);
        if (argc > 7) {
            is_deepslow = std::stoi(argv[7]) != 0;
        }
        if (argc > 8) out_file = argv[8];

        if (!std::isfinite(budget) || budget < 0.0 ||
            !std::isfinite(w_max) || w_max <= 0.0 ||
            !std::isfinite(w_acc) || w_acc <= 0.0 ||
            !std::isfinite(dwell_zenith) || dwell_zenith < 0.0) {
            throw std::invalid_argument(
                "Invalid budget, motion, or dwell parameter.");
        }

        const std::filesystem::path output_path(out_file);
        if (!output_path.parent_path().empty()) {
            std::filesystem::create_directories(
                output_path.parent_path());
        }

        testPathsCoverSTPixelOverlapConfigured(
            tilefile,
            mapfile,
            out_file,
            "overlap+",
            budget,
            dwell_zenith,
            w_max,
            w_acc,
            settle_time,
            is_deepslow,
            source_tile,
            mip_gap,
            time_limit);

        const auto total_end =
            std::chrono::high_resolution_clock::now();

        const std::chrono::duration<double> total_elapsed =
            total_end - total_start;

        std::cout << "Total running time: "
                  << total_elapsed.count()
                  << " seconds" << std::endl;

    } catch (const std::exception& error) {
        std::cerr << "Pixel-overlap path-cover failed: "
                  << error.what() << std::endl;
        return 1;
    }

    return 0;
}