
#include "max_probability_cover.h"
#include "min_telescope_cover.h"
#include "BuildGraph.h"
#include "ILP_gurobi_path_cover_op.h"


#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>



// Distinct name: the existing extractSingleStartProbability remains untouched.
// This helper is only needed when converting a multi-start scalar-prize graph;
// the pixel-overlap ILP itself uses pixel_probs directly.
std::vector<double> extractSingleStartProbabilitySTPixelOverlapMinimal(
        const std::vector<double>& probability_full,
        const std::vector<int>& start_indices) {

    if (probability_full.size() < start_indices.size() + 1) {
        throw std::runtime_error(
            "extractSingleStartProbabilitySTPixelOverlapMinimal: "
            "inconsistent graph sizes.");
    }

    int num_tiles = static_cast<int>(probability_full.size())
                  - static_cast<int>(start_indices.size()) - 1;
    std::vector<double> prob_out(num_tiles + 2, 0.0);

    for (int i = 0; i < num_tiles; ++i)
        prob_out[i] = probability_full[i];

    // prob_out[num_tiles]   = 0.0  (start)
    // prob_out[num_tiles+1] = 0.0  (end)
    return prob_out;
}


// Distinct name: the existing coveredProbabilityUnion remains untouched.
double coveredProbabilityUnionSTPixelOverlapMinimal(
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
                "coveredProbabilityUnionSTPixelOverlapMinimal: "
                "tile index out of range.");
        }
        for (int p : member_pixels[tile]) {
            if (p < 0 || p >= nPixels) {
                throw std::runtime_error(
                    "coveredProbabilityUnionSTPixelOverlapMinimal: "
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
double printPathSTPixelOverlapMinimal(
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
            "printPathSTPixelOverlapMinimal: invalid path endpoints.");
    }

    double sum_probability = coveredProbabilityUnionSTPixelOverlapMinimal(
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


// Distinct name: the existing save_result_overlap remains untouched.
void saveResultSTPixelOverlapMinimal(
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
        const PixelOverlapILPResult& result,
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
               "SumProb,TotalCost,TimeSec,IncumbentObj,BestBound,MIPGap,"
               "GurobiStatus,SolutionCount,Path\n";
    }

    const std::vector<int>& best_path = result.path;
    double sum_probability = 0.0;
    double total_cost = 0.0;
    std::size_t num_selected_tiles = 0;

    if (!best_path.empty()) {
        if (best_path.front() != start_idx || best_path.back() != end_idx) {
            throw std::runtime_error(
                "saveResultSTPixelOverlapMinimal: invalid path endpoints.");
        }
        sum_probability = coveredProbabilityUnionSTPixelOverlapMinimal(
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
            << total_cost << ","
            << elapsed_time << ","
            << result.incumbent_objective << ","
            << result.best_bound << ","
            << result.relative_gap << ","
            << result.gurobi_status << ","
            << result.solution_count << ",";

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


enum class PixelOverlapILPFormulation {
    MTZ,
    ContinuousTime,
};


// Configurable entry point used by main.  Its distinct name preserves the
// original testPathsCoverSTPixelOverlapMinimal interface below.
// budget must already include the padding returned by buildGraph.
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
        double time_limit,
        PixelOverlapILPFormulation formulation) {

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

    // The direct s->t edge consists of the terminal padding. Since budget is
    // already padded, it must at least permit this edge.
    if (costs[start_idx][end_idx] > budget + 1e-9) {
        throw std::runtime_error(
            "testPathsCoverSTPixelOverlapMinimal: budget is smaller than "
            "the direct start-to-end padding cost.");
    }

    std::cout << "start: " << start_idx
              << ", end: " << end_idx
              << ", retained tiles: " << num_tiles
              << ", padding: " << padding << std::endl;

    std::chrono::high_resolution_clock::time_point start;
    std::chrono::high_resolution_clock::time_point end;
    std::chrono::duration<double> elapsed_seconds;

    const bool use_continuous_time =
        formulation == PixelOverlapILPFormulation::ContinuousTime;
    const std::string formulation_name =
        use_continuous_time ? "Continuous Time" : "MTZ";
    const std::string method_suffix =
        use_continuous_time ? "continuous-time" : "mtz";

    std::cout << "*********ILP Path Cover (" << formulation_name
              << ", Pixel Overlap)*********" << std::endl;
    start = std::chrono::high_resolution_clock::now();

    PixelOverlapILPResult result = use_continuous_time
        ? gurobiPathCoverPixelOverlapContinuousTime(
              costs, pixel_probs, member_pixels,
              start_idx, end_idx, budget,
              accu_thr, time_limit)
        : gurobiPathCoverPixelOverlapMTZ(
              costs, pixel_probs, member_pixels,
              start_idx, end_idx, budget,
              accu_thr, time_limit);

    end = std::chrono::high_resolution_clock::now();
    elapsed_seconds = end - start;

    std::cout << "Incumbent Objective: "
              << result.incumbent_objective << std::endl;
    std::cout << "Best Bound: " << result.best_bound << std::endl;
    std::cout << "MIP Gap: " << result.relative_gap << std::endl;
    std::cout << "Gurobi Status: " << result.gurobi_status << std::endl;

    printPathSTPixelOverlapMinimal(
        costs, member_pixels, pixel_probs,
        ranks, result.path, start_idx, end_idx);

    std::cout << "running time (wallclock): "
              << elapsed_seconds.count() << " seconds" << std::endl;

    saveResultSTPixelOverlapMinimal(
        out_file, overlap + "path+cover+" + method_suffix,
        tilefile, mapfile, budget, w_max,
        costs, member_pixels, pixel_probs,
        ranks, result, start_idx, end_idx,
        elapsed_seconds.count());
}


// Preserve the earlier helper exactly for any existing callers.
void testPathsCoverSTPixelOverlapMinimal(
        std::string tilefile,
        std::string mapfile,
        std::string out_file,
        std::string overlap,
        double budget,
        double dwell_time,
        double w_max,
        double w_acc,
        double settle_time,
        bool is_deepslow) {

    testPathsCoverSTPixelOverlapConfigured(
        std::move(tilefile), std::move(mapfile), std::move(out_file),
        std::move(overlap), budget, dwell_time,
        w_max, w_acc, settle_time, is_deepslow,
        0, 1e-6, 3600.0, PixelOverlapILPFormulation::MTZ);
}


int main(int argc, char** argv) {
    const auto total_start = std::chrono::high_resolution_clock::now();

    // Defaults match the existing command-line program.
    std::string out_file = "../results5/out.csv";
    std::string tilefile = "../data/tilings_small/6.9x6.9_tiling.csv";
    std::string mapfile =
        "../data/LIGO_maps/txt_files_6.9x6.9/small_maps/"
        "GW200316_215756_50.txt";

    int source_tile = -2;  // Any negative value selects the deterministic seed.
    double w_max = 10.0;
    double w_acc = 10.0;
    double dwell_zenith = 2.0;
    double settle_time = 0.0;
    bool is_deepslow = false;
    double budget = 60.0;  // Already includes terminal padding.
    int m = 4;             // Existing method id: 4 means Gurobi.

    double mipGap = 0.001;
    double timeLimit = 36000.0;
    std::string formulation = "mtz";

    if (argc > 14) {
        std::cerr
            << "Usage: " << argv[0]
            << " [MAP_FILE TILE_FILE"
               " [BUDGET W_MAX W_ACC DWELL_ZENITH IS_DEEPSLOW METHOD_ID"
               " [FORMULATION [SOURCE_TILE [MIP_GAP [TIME_LIMIT"
               " [OUTPUT_CSV]]]]]]]\n"
            << "  METHOD_ID must be 4 for this Gurobi executable.\n"
            << "  FORMULATION is mtz or continuous-time.\n";
        return 2;
    }

    try {
        if (argc > 2) {
            mapfile = argv[1];
            tilefile = argv[2];
        }
        if (argc > 8) {
            budget = std::stod(argv[3]);
            w_max = std::stod(argv[4]);
            w_acc = std::stod(argv[5]);
            dwell_zenith = std::stod(argv[6]);
            is_deepslow = std::stoi(argv[7]) != 0;
            m = std::stoi(argv[8]);
        }
        if (argc > 9) formulation = argv[9];
        if (argc > 10) source_tile = std::stoi(argv[10]);
        if (argc > 11) mipGap = std::stod(argv[11]);
        if (argc > 12) timeLimit = std::stod(argv[12]);
        if (argc > 13) out_file = argv[13];

        if (m != 4) {
            throw std::invalid_argument(
                "This executable implements method 4 (Gurobi) only.");
        }
        if (!std::isfinite(budget) || budget < 0.0 ||
            !std::isfinite(dwell_zenith) || dwell_zenith < 0.0 ||
            !std::isfinite(w_max) || w_max <= 0.0 ||
            !std::isfinite(w_acc) || w_acc <= 0.0 ||
            !std::isfinite(settle_time) || settle_time < 0.0 ||
            !std::isfinite(mipGap) || mipGap < 0.0 ||
            !std::isfinite(timeLimit) || timeLimit <= 0.0) {
            throw std::invalid_argument(
                "Invalid budget, motion, dwell, MIP-gap, or time-limit "
                "parameter.");
        }

        PixelOverlapILPFormulation selected_formulation;
        if (formulation == "mtz" || formulation == "MTZ") {
            selected_formulation = PixelOverlapILPFormulation::MTZ;
        } else if (formulation == "continuous-time" ||
                   formulation == "continuous" ||
                   formulation == "time") {
            selected_formulation =
                PixelOverlapILPFormulation::ContinuousTime;
        } else {
            throw std::invalid_argument(
                "FORMULATION must be mtz or continuous-time.");
        }

        const std::filesystem::path output_path(out_file);
        if (!output_path.parent_path().empty()) {
            std::filesystem::create_directories(output_path.parent_path());
        }

        testPathsCoverSTPixelOverlapConfigured(
            tilefile, mapfile, out_file, "overlap+",
            budget, dwell_zenith, w_max, w_acc, settle_time,
            is_deepslow, source_tile, mipGap, timeLimit,
            selected_formulation);

        const auto total_end = std::chrono::high_resolution_clock::now();
        const std::chrono::duration<double> total_elapsed =
            total_end - total_start;
        std::cout << "Total running time: " << total_elapsed.count()
                  << " seconds" << std::endl;
    } catch (const std::exception& error) {
        std::cerr << "Pixel-overlap path-cover failed: "
                  << error.what() << std::endl;
        return 1;
    }

    return 0;
}