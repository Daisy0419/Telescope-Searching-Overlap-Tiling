#include "BuildGraph.h"
#include "ILP_gurobi_path_cover_op.h"

#include "AntColonyPixelOverlap.h"
#include "GeneticPixelOverlap.h"
#include "GreedyPixelOverlap.h"
#include "SimulatedAnnealingPixelOverlap.h"
#include "SubmodularSelection.h"
#include "GCP_overlap.h"
#include "PixelOverlapPath.h"
#include "recursive_greedy_op.h"


#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

struct BenchInstance {
    std::vector<std::vector<double>> costs;
    std::vector<double> pixel_probabilities;
    std::vector<std::vector<int>> member_pixels;
    std::vector<double> tile_probabilities;
    std::vector<double> dwell_times;
    std::vector<int> ranks;

    int start_node = 0;
    int end_node = 0;
    int tile_count = 0;

    double padding = 0.0;
    double observing_budget = 0.0;
    double solver_budget = 0.0;
};

struct SolverResult {
    std::string method;
    std::vector<int> path;

    double coverage = 0.0;
    double padded_cost = 0.0;
    double duration = 0.0;
    double bound = std::numeric_limits<double>::quiet_NaN();
    double seconds = 0.0;

    std::string status = "ok";
};

// Return true when candidate is a better feasible result than the incumbent.
// Coverage is primary; path cost breaks coverage ties.
static bool resultIsBetter(const SolverResult& candidate,
                           double incumbent_coverage,
                           double incumbent_cost) {
    const double tolerance = 1.0e-12;
    if (candidate.status != "ok") {
        return false;
    }
    if (candidate.coverage > incumbent_coverage + tolerance) {
        return true;
    }
    return std::abs(candidate.coverage - incumbent_coverage) <= tolerance &&
           candidate.padded_cost < incumbent_cost - tolerance;
}

// Validate a returned route and score it with the union-of-pixels objective.
static void validateAndScore(const BenchInstance& instance,
                             const CoverageOracle& oracle,
                             SolverResult& result) {
    result.coverage = 0.0;
    result.padded_cost = 0.0;
    result.duration = 0.0;

    if (result.path.empty()) {
        result.status = "no path";
        return;
    }
    if (result.path.size() < 2 ||
        result.path.front() != instance.start_node ||
        result.path.back() != instance.end_node) {
        result.status = "bad endpoints";
        return;
    }

    std::vector<char> seen(instance.costs.size(), 0);
    for (int node : result.path) {
        if (node < 0 ||
            node >= static_cast<int>(instance.costs.size())) {
            result.status = "node out of range";
            return;
        }
        if (seen[node] != 0) {
            result.status = "revisits a node";
            return;
        }
        seen[node] = 1;
    }

    for (std::size_t position = 1;
         position < result.path.size(); ++position) {
        result.padded_cost +=
            instance.costs[result.path[position - 1]][result.path[position]];
    }
    result.duration = result.padded_cost - instance.padding;
    result.coverage = oracle.rewardOf(result.path);

    if (result.padded_cost > instance.solver_budget + 1.0e-6) {
        result.status = "over budget";
    } else if (result.duration < -1.0e-6) {
        result.status = "cost below padding";
    }
}

// Time one solver, catch its exception, then validate and score its path.
static SolverResult runSolver(
    const BenchInstance& instance,
    const CoverageOracle& oracle,
    const std::string& name,
    const std::function<std::vector<int>()>& solve) {

    SolverResult result;
    result.method = name;

    const std::chrono::high_resolution_clock::time_point started =
        std::chrono::high_resolution_clock::now();
    try {
        result.path = solve();
    } catch (const std::exception& error) {
        result.status = std::string("threw: ") + error.what();
    } catch (...) {
        result.status = "threw an unknown exception";
    }
    result.seconds = std::chrono::duration<double>(
        std::chrono::high_resolution_clock::now() - started).count();

    if (result.status == "ok") {
        validateAndScore(instance, oracle, result);
    }
    return result;
}

// Convert the single-telescope ILP result wrapper to one path.
static std::vector<int> firstPathOf(
    const std::pair<std::vector<std::vector<int>>, double>& ilp_result) {
    if (ilp_result.first.size() > 1) {
        throw std::runtime_error(
            "single-telescope ILP returned more than one path");
    }
    if (ilp_result.first.empty()) {
        return {};
    }
    return ilp_result.first.front();
}

// Escape a string for one RFC-4180-style CSV field.
static std::string csvField(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos) {
        return value;
    }

    std::string escaped = "\"";
    for (char character : value) {
        if (character == '\"') {
            escaped += "\"\"";
        } else {
            escaped += character;
        }
    }
    escaped += '\"';
    return escaped;
}

// Convert a validated path to the original tile ranks used in output files.
// Invalid solver paths are left blank instead of being indexed.
static std::string outputPath(const BenchInstance& instance,
                              const SolverResult& result) {
    if (result.status != "ok" || result.path.size() < 2) {
        return "";
    }

    std::ostringstream path_stream;
    for (std::size_t position = 1;
         position + 1 < result.path.size(); ++position) {
        const int tile = result.path[position];
        int output_tile = tile;
        if (tile >= 0 &&
            tile < static_cast<int>(instance.ranks.size())) {
            output_tile = instance.ranks[tile];
        }

        if (position > 1) {
            path_stream << ' ';
        }
        path_stream << output_tile;
    }
    return path_stream.str();
}

// Print a compact comparison table. Percent-best is relative to the best
// feasible solution found by any solver in this run.
static void printTable(const BenchInstance& instance,
                       const std::vector<SolverResult>& results) {
    double best_coverage = 0.0;
    for (const SolverResult& result : results) {
        if (result.status == "ok") {
            best_coverage = std::max(best_coverage, result.coverage);
        }
    }

    std::cout << "\n"
              << std::left << std::setw(34) << "Method"
              << std::right << std::setw(7) << "Tiles"
              << std::setw(12) << "Coverage"
              << std::setw(9) << "%best"
              << std::setw(11) << "Duration"
              << std::setw(12) << "Bound"
              << std::setw(11) << "Seconds"
              << "  Status\n"
              << std::string(110, '-') << "\n";

    for (const SolverResult& result : results) {
        const std::size_t selected_tiles = result.path.size() >= 2
            ? result.path.size() - 2
            : 0;

        std::cout << std::left << std::setw(34) << result.method
                  << std::right << std::setw(7) << selected_tiles
                  << std::setw(12) << std::fixed << std::setprecision(6)
                  << result.coverage
                  << std::setw(8) << std::setprecision(1)
                  << (best_coverage > 0.0
                          ? 100.0 * result.coverage / best_coverage
                          : 0.0)
                  << "%"
                  << std::setw(11) << result.duration;

        if (std::isfinite(result.bound)) {
            std::cout << std::setw(12) << std::setprecision(6)
                      << result.bound;
        } else {
            std::cout << std::setw(12) << "-";
        }

        std::cout << std::setw(11) << std::setprecision(2)
                  << result.seconds
                  << "  "
                  << (result.status == "ok" ? "" : result.status)
                  << "\n";
    }

    std::cout << std::string(110, '-') << "\n"
              << "observing budget D = " << std::setprecision(1)
              << instance.observing_budget
              << " s, padding = " << instance.padding
              << " s, solver budget = " << instance.solver_budget
              << " s\n";
}

// Append all results to a CSV file. Only validated paths are serialized.
static void appendCsv(const std::string& filename,
                      const BenchInstance& instance,
                      const std::string& tiling_file,
                      const std::string& probability_file,
                      double maximum_slew_rate,
                      const std::vector<SolverResult>& results) {
    const bool file_exists = std::filesystem::exists(filename);
    const bool write_header = !file_exists ||
        std::filesystem::file_size(filename) == 0;

    std::ofstream output(filename, std::ios::app);
    if (!output.is_open()) {
        throw std::runtime_error(
            "could not open results file: " + filename);
    }
    output << std::setprecision(17);

    if (write_header) {
        output
            << "Method,TilingData,LikelihoodData,ObservingBudget,Padding,"
               "SlewRate,NumTilesRetained,NumTilesInPath,Coverage,BestBound,"
               "Duration,PaddedCost,TimeSec,Status,Path\n";
    }

    for (const SolverResult& result : results) {
        const std::size_t selected_tiles = result.path.size() >= 2
            ? result.path.size() - 2
            : 0;

        output
            << csvField(result.method) << ','
            << csvField(std::filesystem::path(tiling_file)
                            .filename().string()) << ','
            << csvField(std::filesystem::path(probability_file)
                            .filename().string()) << ','
            << instance.observing_budget << ','
            << instance.padding << ','
            << maximum_slew_rate << ','
            << instance.tile_count << ','
            << selected_tiles << ','
            << result.coverage << ','
            << result.bound << ','
            << result.duration << ','
            << result.padded_cost << ','
            << result.seconds << ','
            << csvField(result.status) << ','
            << csvField(outputPath(instance, result)) << '\n';
    }
}

int main(int argc, char** argv) {
    const std::chrono::high_resolution_clock::time_point total_started =
        std::chrono::high_resolution_clock::now();

    std::string output_file = "../results/out.csv";
    std::string tiling_file = "../data/tilings_small/6.9x6.9_tiling.csv";
    std::string probability_file =
        "../data/test_maps/GW200316_215756_50.txt";

    double observing_budget = 60.0;
    double maximum_slew_rate = 10.0;
    double maximum_slew_acceleration = 10.0;
    double zenith_dwell = 2.0;
    bool is_deep_slow = false;
    double ilp_time_limit = 1800.0;

    const double settle_time = 0.0;
    // The initial telescope pointing is independent of the zenith direction.
    // A negative value keeps the graph builder's deterministic random choice.
    // If the real initial pointing is a retained tile, pass that tile's index.
    const int source_tile = -1;
    const double mip_gap = 0.001;

    if (argc > 10) {
        std::cerr
            << "Usage: " << argv[0]
            << " [MAP TILES D W_MAX W_ACC DWELL_ZENITH IS_DEEPSLOW "
               "OUT_CSV ILP_TIME_LIMIT]\n";
        return 2;
    }

    try {
        if (argc > 1) probability_file = argv[1];
        if (argc > 2) tiling_file = argv[2];
        if (argc > 3) observing_budget = std::stod(argv[3]);
        if (argc > 4) maximum_slew_rate = std::stod(argv[4]);
        if (argc > 5) maximum_slew_acceleration = std::stod(argv[5]);
        if (argc > 6) zenith_dwell = std::stod(argv[6]);
        if (argc > 7) is_deep_slow = std::stoi(argv[7]) != 0;
        if (argc > 8) output_file = argv[8];
        if (argc > 9) ilp_time_limit = std::stod(argv[9]);

        if (!std::isfinite(observing_budget) || observing_budget < 0.0 ||
            !std::isfinite(maximum_slew_rate) || maximum_slew_rate <= 0.0 ||
            !std::isfinite(maximum_slew_acceleration) ||
                maximum_slew_acceleration <= 0.0 ||
            !std::isfinite(zenith_dwell) || zenith_dwell < 0.0 ||
            !std::isfinite(ilp_time_limit) || ilp_time_limit < 0.0) {
            throw std::invalid_argument(
                "invalid budget, motion, dwell, or time-limit parameter");
        }

        const std::filesystem::path output_path(output_file);
        if (!output_path.parent_path().empty()) {
            std::filesystem::create_directories(
                output_path.parent_path());
        }

        // // Build the instance once. Every solver receives these same arrays.
        // BenchInstance instance;
        // const std::tuple<int, int, double> graph = buildGraphwithPixelsOP(
        //                                                     tiling_file,
        //                                                     probability_file,
        //                                                     instance.costs,
        //                                                     instance.dwell_times,
        //                                                     instance.tile_probabilities,
        //                                                     instance.ranks,
        //                                                     instance.member_pixels,
        //                                                     instance.pixel_probabilities,
        //                                                     zenith_dwell,
        //                                                     maximum_slew_rate,
        //                                                     maximum_slew_acceleration,
        //                                                     settle_time,
        //                                                     is_deep_slow,
        //                                                     instance.tile_count,
        //                                                     source_tile);

        // instance.start_node = std::get<0>(graph);
        // instance.end_node = std::get<1>(graph);
        // instance.padding = std::get<2>(graph);
        // instance.observing_budget = observing_budget;
        // instance.solver_budget = observing_budget + instance.padding;

        // if (instance.costs.size() < 2 ||
        //     instance.member_pixels.size() != instance.costs.size()) {
        //     throw std::invalid_argument(
        //         "graph and member_pixels have inconsistent node counts");
        // }
        // for (const std::vector<double>& row : instance.costs) {
        //     if (row.size() != instance.costs.size()) {
        //         throw std::invalid_argument("cost matrix is not square");
        //     }
        // }
        // if (instance.start_node < 0 ||
        //     instance.start_node >= static_cast<int>(instance.costs.size()) ||
        //     instance.end_node < 0 ||
        //     instance.end_node >= static_cast<int>(instance.costs.size()) ||
        //     instance.start_node == instance.end_node) {
        //     throw std::invalid_argument("graph returned invalid endpoints");
        // }
        // if (!instance.member_pixels[instance.start_node].empty() ||
        //     !instance.member_pixels[instance.end_node].empty()) {
        //     throw std::invalid_argument(
        //         "dummy start and end nodes must cover no pixels");
        // }
        // if (instance.costs[instance.start_node][instance.end_node] >
        //     instance.solver_budget + 1.0e-9) {
        //     throw std::runtime_error(
        //         "budget is smaller than the direct start-to-end cost");
        // }

        // const CoverageOracle oracle(
        //     instance.pixel_probabilities, instance.member_pixels);
        // if (oracle.nodeCount() !=
        //     static_cast<int>(instance.costs.size())) {
        //     throw std::invalid_argument(
        //         "member_pixels must have one list per graph node");
        // }

        // double retained_probability = 0.0;
        // for (double probability : instance.pixel_probabilities) {
        //     retained_probability += probability;
        // }

        // std::cout << "start " << instance.start_node
        //           << ", end " << instance.end_node
        //           << ", retained tiles " << instance.tile_count
        //           << ", retained pixels "
        //           << instance.pixel_probabilities.size()
        //           << ", retained pixel probability "
        //           << std::fixed << std::setprecision(6)
        //           << retained_probability << "\n";

        // double solver_budget = instance.solver_budget;
        // std::vector<SolverResult> results;

        // // The two deterministic greedy routes also seed GA and SA.
        // results.push_back(runSolver(
        //     instance, oracle, "greedy-marginal", [&]() {
        //         return prizeGreedyPathTwoFixedPixelOverlap(
        //             instance.costs,
        //             instance.pixel_probabilities,
        //             instance.member_pixels,
        //             solver_budget,
        //             instance.start_node,
        //             instance.end_node);
        //     }));

        // results.push_back(runSolver(
        //     instance, oracle, "greedy-ratio", [&]() {
        //         return prizeRatioGreedyPathTwoFixedPixelOverlap(
        //             instance.costs,
        //             instance.pixel_probabilities,
        //             instance.member_pixels,
        //             solver_budget,
        //             instance.start_node,
        //             instance.end_node);
        //     }));

        // std::vector<std::vector<int>> genetic_seeds;
        // std::vector<int> annealing_seed;
        // double seed_coverage = -1.0;
        // double seed_cost = std::numeric_limits<double>::infinity();

        // for (std::size_t index = 0;
        //      index < 2 && index < results.size(); ++index) {
        //     if (results[index].status == "ok") {
        //         genetic_seeds.push_back(results[index].path);
        //         if (resultIsBetter(
        //                 results[index], seed_coverage, seed_cost)) {
        //             annealing_seed = results[index].path;
        //             seed_coverage = results[index].coverage;
        //             seed_cost = results[index].padded_cost;
        //         }
        //     }
        // }


        // results.push_back(runSolver(
        //     instance, oracle, "GCP", [&]() {
        //         return GCP_Coverage(
        //             instance.costs,
        //             instance.pixel_probabilities,
        //             instance.member_pixels,
        //             solver_budget,
        //             instance.start_node,
        //             instance.end_node);
        //     }));


        // results.push_back(runSolver(
        //     instance, oracle, "simulated-annealing", [&]() {
        //         return simulated_annealing_optimization_pixel_overlap(
        //             instance.costs,
        //             instance.pixel_probabilities,
        //             instance.member_pixels,
        //             solver_budget,
        //             instance.start_node,
        //             instance.end_node,
        //             annealing_seed);
        //     }));

        // results.push_back(runSolver(
        //     instance, oracle, "genetic", [&]() {
        //         return genetic_optimization_st_pixel_overlap(
        //             instance.costs,
        //             instance.pixel_probabilities,
        //             instance.member_pixels,
        //             solver_budget,
        //             instance.start_node,
        //             instance.end_node,
        //             genetic_seeds);
        //     }));

        // // results.push_back(runSolver(
        // //     instance, oracle, "ant-colony", [&]() {
        // //         AntColonyPixelOverlap ant_colony(
        // //             instance.costs,
        // //             instance.pixel_probabilities,
        // //             instance.member_pixels);
        // //         return ant_colony.ant_colony_optimization_pixel_overlap(
        // //             solver_budget,
        // //             instance.start_node,
        // //             instance.end_node);
        // //     }));


        // // results.push_back(runSolver(
        // //     instance, oracle, "RG-QP", [&]() {
        // //         RecursiveGreedyParams parameters;
        // //         const int relevant_arc_bound =
        // //             std::max(1, instance.tile_count + 1);
        // //         parameters.max_depth = std::min(
        // //             50,
        // //             static_cast<int>(std::ceil(
        // //                 1.0 + std::log2(
        // //                     static_cast<double>(relevant_arc_bound)))));
        // //         parameters.reward_levels = 10000;
        // //         parameters.time_quantum = 0.001;
        // //         parameters.time_limit = 0;
        // //         parameters.iterative_deepening = true;
        // //         parameters.candidate_limit = 0;
        // //         parameters.verbose = false;

        // //         const RecursiveGreedyResult rg_result =
        // //             recursiveGreedyOrienteering(
        // //                 instance.costs,
        // //                 instance.pixel_probabilities,
        // //                 instance.member_pixels,
        // //                 instance.start_node,
        // //                 instance.end_node,
        // //                 solver_budget,
        // //                 instance.padding,
        // //                 parameters);
        // //         return rg_result.path;
        // //     }));


        // // // Select the best feasible heuristic route for both ILP warm starts.
        // // std::vector<int> warm_start;
        // // double warm_coverage = -1.0;
        // // double warm_cost = std::numeric_limits<double>::infinity();
        // // for (const SolverResult& result : results) {
        // //     if (resultIsBetter(result, warm_coverage, warm_cost)) {
        // //         warm_start = result.path;
        // //         warm_coverage = result.coverage;
        // //         warm_cost = result.padded_cost;
        // //     }
        // // }

        // // if (warm_start.empty()) {
        // //     std::cout << "warm start: none\n";
        // // } else {
        // //     std::cout << "warm start: " << warm_start.size() - 2
        // //               << " tiles, coverage " << warm_coverage
        // //               << ", padded cost " << warm_cost << "\n";
        // // }

        // // double mtz_bound = std::numeric_limits<double>::quiet_NaN();
        // // SolverResult mtz = runSolver(
        // //     instance, oracle, "ILP MTZ", [&]() {
        // //         const std::pair<std::vector<std::vector<int>>, double> output =
        // //             gurobiPathCoverPixelOverlapMTZ(
        // //                 instance.costs,
        // //                 instance.pixel_probabilities,
        // //                 instance.member_pixels,
        // //                 instance.start_node,
        // //                 instance.end_node,
        // //                 solver_budget,
        // //                 mip_gap,
        // //                 ilp_time_limit,
        // //                 warm_start);
        // //         mtz_bound = output.second;
        // //         return firstPathOf(output);
        // //     });
        // // mtz.bound = mtz_bound;
        // // results.push_back(std::move(mtz));

        // // double continuous_bound =
        // //     std::numeric_limits<double>::quiet_NaN();
        // // SolverResult continuous_time = runSolver(
        // //     instance, oracle, "ILP continuous-time", [&]() {
        // //         const std::pair<std::vector<std::vector<int>>, double> output =
        // //             gurobiPathCoverPixelOverlapContinuousTime(
        // //                 instance.costs,
        // //                 instance.pixel_probabilities,
        // //                 instance.member_pixels,
        // //                 instance.start_node,
        // //                 instance.end_node,
        // //                 solver_budget,
        // //                 mip_gap,
        // //                 ilp_time_limit,
        // //                 warm_start);
        // //         continuous_bound = output.second;
        // //         return firstPathOf(output);
        // //     });
        // // continuous_time.bound = continuous_bound;
        // // results.push_back(std::move(continuous_time));


        // // Save the full-instance results before constructing the
        // // reduced instance. Otherwise they would later be interpreted using
        // // the reduced graph's node indices and ranks.
        // printTable(instance, results);
        // appendCsv(
        //     output_file,
        //     instance,
        //     tiling_file,
        //     probability_file,
        //     maximum_slew_rate,
        //     results);

        // Build a separate top-N instance. In particular, all vectors passed
        // as output parameters are empty when readTiles() is called again.
        const int maximum_tiles = 50;
        BenchInstance reduced_instance;
        const std::tuple<int, int, double> graph2 = buildGraphwithPixelsOP(
                                                            tiling_file,
                                                            probability_file,
                                                            reduced_instance.costs,
                                                            reduced_instance.dwell_times,
                                                            reduced_instance.tile_probabilities,
                                                            reduced_instance.ranks,
                                                            reduced_instance.member_pixels,
                                                            reduced_instance.pixel_probabilities,
                                                            zenith_dwell,
                                                            maximum_slew_rate,
                                                            maximum_slew_acceleration,
                                                            settle_time,
                                                            is_deep_slow,
                                                            maximum_tiles,
                                                            reduced_instance.tile_count,
                                                            source_tile);

        reduced_instance.start_node = std::get<0>(graph2);
        reduced_instance.end_node = std::get<1>(graph2);
        reduced_instance.padding = std::get<2>(graph2);
        reduced_instance.observing_budget = observing_budget;
        reduced_instance.solver_budget =
            observing_budget + reduced_instance.padding;

        if (reduced_instance.costs.size() < 2 ||
            reduced_instance.member_pixels.size() !=
                reduced_instance.costs.size()) {
            throw std::invalid_argument(
                "graph and member_pixels have inconsistent node counts");
        }
        for (const std::vector<double>& row : reduced_instance.costs) {
            if (row.size() != reduced_instance.costs.size()) {
                throw std::invalid_argument("cost matrix is not square");
            }
        }
        if (reduced_instance.start_node < 0 ||
            reduced_instance.start_node >=
                static_cast<int>(reduced_instance.costs.size()) ||
            reduced_instance.end_node < 0 ||
            reduced_instance.end_node >=
                static_cast<int>(reduced_instance.costs.size()) ||
            reduced_instance.start_node == reduced_instance.end_node) {
            throw std::invalid_argument("graph returned invalid endpoints");
        }
        if (!reduced_instance.member_pixels[reduced_instance.start_node].empty() ||
            !reduced_instance.member_pixels[reduced_instance.end_node].empty()) {
            throw std::invalid_argument(
                "dummy start and end nodes must cover no pixels");
        }
        if (reduced_instance.costs[reduced_instance.start_node]
                                  [reduced_instance.end_node] >
            reduced_instance.solver_budget + 1.0e-9) {
            throw std::runtime_error(
                "budget is smaller than the direct start-to-end cost");
        }

        const CoverageOracle reduced_oracle(
            reduced_instance.pixel_probabilities,
            reduced_instance.member_pixels);
        if (reduced_oracle.nodeCount() !=
            static_cast<int>(reduced_instance.costs.size())) {
            throw std::invalid_argument(
                "member_pixels must have one list per graph node");
        }

        double reduced_retained_probability = 0.0;
        for (double probability : reduced_instance.pixel_probabilities) {
            reduced_retained_probability += probability;
        }

        std::cout << "start " << reduced_instance.start_node
                  << ", end " << reduced_instance.end_node
                  << ", retained tiles " << reduced_instance.tile_count
                  << ", retained pixels "
                  << reduced_instance.pixel_probabilities.size()
                  << ", retained pixel probability "
                  << std::fixed << std::setprecision(6)
                  << reduced_retained_probability << "\n";

        const double reduced_solver_budget = reduced_instance.solver_budget;

        std::vector<SolverResult> reduced_results;


        reduced_results.push_back(runSolver(
            reduced_instance, reduced_oracle, "RG-QP", [&]() {
                RecursiveGreedyParams parameters;
                const int relevant_arc_bound =
                    std::max(1, reduced_instance.tile_count + 1);
                parameters.max_depth = std::min(
                    4,
                    static_cast<int>(std::ceil(
                        1.0 + std::log2(
                            static_cast<double>(relevant_arc_bound)))));
                parameters.reward_levels = 100;
                parameters.time_quantum = 0.01;
                parameters.time_limit = 7200;
                parameters.iterative_deepening = true;
                parameters.candidate_limit = 0;
                parameters.verbose = false;

                const RecursiveGreedyResult rg_result =
                    recursiveGreedyOrienteering(
                        reduced_instance.costs,
                        reduced_instance.pixel_probabilities,
                        reduced_instance.member_pixels,
                        reduced_instance.start_node,
                        reduced_instance.end_node,
                        reduced_solver_budget,
                        reduced_instance.padding,
                        parameters);
                return rg_result.path;
            }));


        // // Select the best feasible heuristic route for both ILP warm starts.
        // std::vector<int> warm_start;
        // double warm_coverage = -1.0;
        // double warm_cost = std::numeric_limits<double>::infinity();
        // for (const SolverResult& result : results) {
        //     if (resultIsBetter(result, warm_coverage, warm_cost)) {
        //         warm_start = result.path;
        //         warm_coverage = result.coverage;
        //         warm_cost = result.padded_cost;
        //     }
        // }

        // if (warm_start.empty()) {
        //     std::cout << "warm start: none\n";
        // } else {
        //     std::cout << "warm start: " << warm_start.size() - 2
        //               << " tiles, coverage " << warm_coverage
        //               << ", padded cost " << warm_cost << "\n";
        // }

        // double mtz_bound = std::numeric_limits<double>::quiet_NaN();
        // SolverResult mtz = runSolver(
        //     instance, oracle, "ILP MTZ", [&]() {
        //         const std::pair<std::vector<std::vector<int>>, double> output =
        //             gurobiPathCoverPixelOverlapMTZ(
        //                 instance.costs,
        //                 instance.pixel_probabilities,
        //                 instance.member_pixels,
        //                 instance.start_node,
        //                 instance.end_node,
        //                 solver_budget,
        //                 mip_gap,
        //                 ilp_time_limit,
        //                 warm_start);
        //         mtz_bound = output.second;
        //         return firstPathOf(output);
        //     });
        // mtz.bound = mtz_bound;
        // results.push_back(std::move(mtz));

        // std::vector<int> warm_start = {};
        // double continuous_bound =
        //     std::numeric_limits<double>::quiet_NaN();
        // SolverResult continuous_time = runSolver(
        //     reduced_instance,
        //     reduced_oracle,
        //     "ILP continuous-time",
        //     [&]() {
        //         const std::pair<std::vector<std::vector<int>>, double> output =
        //             gurobiPathCoverPixelOverlapContinuousTime(
        //                 reduced_instance.costs,
        //                 reduced_instance.pixel_probabilities,
        //                 reduced_instance.member_pixels,
        //                 reduced_instance.start_node,
        //                 reduced_instance.end_node,
        //                 reduced_solver_budget,
        //                 mip_gap,
        //                 ilp_time_limit,
        //                 warm_start);
        //         continuous_bound = output.second;
        //         return firstPathOf(output);
        //     });
        // continuous_time.bound = continuous_bound;
        // reduced_results.push_back(std::move(continuous_time));


        printTable(reduced_instance, reduced_results);
        appendCsv(
            output_file,
            reduced_instance,
            tiling_file,
            probability_file,
            maximum_slew_rate,
            reduced_results);

        const double total_seconds = std::chrono::duration<double>(
            std::chrono::high_resolution_clock::now() - total_started).count();
        std::cout << "\nTotal running time: " << std::setprecision(2)
                  << total_seconds << " seconds\n";
    } catch (const std::exception& error) {
        std::cerr << "Pixel-overlap benchmark failed: "
                  << error.what() << "\n";
        return 1;
    }

    return 0;
}


// int main(int argc, char** argv) {
//     const std::chrono::high_resolution_clock::time_point total_started =
//         std::chrono::high_resolution_clock::now();

//     std::string output_file = "../results/out.csv";
//     std::string tiling_file = "../data/tilings_small/6.9x6.9_tiling.csv";
//     std::string probability_file =
//         "../data/test_maps/GW200316_215756_50.txt";

//     double observing_budget = 60.0;
//     double maximum_slew_rate = 10.0;
//     double maximum_slew_acceleration = 10.0;
//     double zenith_dwell = 2.0;
//     bool is_deep_slow = false;
//     double ilp_time_limit = 14400.0;

//     const double settle_time = 0.0;
//     const int source_tile = -1;
//     const double mip_gap = 0.001;

//     if (argc > 10) {
//         std::cerr
//             << "Usage: " << argv[0]
//             << " [MAP TILES D W_MAX W_ACC DWELL_ZENITH IS_DEEPSLOW "
//                "OUT_CSV ILP_TIME_LIMIT]\n";
//         return 2;
//     }

//     try {
//         if (argc > 1) probability_file = argv[1];
//         if (argc > 2) tiling_file = argv[2];
//         if (argc > 3) observing_budget = std::stod(argv[3]);
//         if (argc > 4) maximum_slew_rate = std::stod(argv[4]);
//         if (argc > 5) maximum_slew_acceleration = std::stod(argv[5]);
//         if (argc > 6) zenith_dwell = std::stod(argv[6]);
//         if (argc > 7) is_deep_slow = std::stoi(argv[7]) != 0;
//         if (argc > 8) output_file = argv[8];
//         if (argc > 9) ilp_time_limit = std::stod(argv[9]);

//         if (!std::isfinite(observing_budget) || observing_budget < 0.0 ||
//             !std::isfinite(maximum_slew_rate) || maximum_slew_rate <= 0.0 ||
//             !std::isfinite(maximum_slew_acceleration) ||
//                 maximum_slew_acceleration <= 0.0 ||
//             !std::isfinite(zenith_dwell) || zenith_dwell < 0.0 ||
//             !std::isfinite(ilp_time_limit) || ilp_time_limit < 0.0) {
//             throw std::invalid_argument(
//                 "invalid budget, motion, dwell, or time-limit parameter");
//         }

//         const std::filesystem::path output_path(output_file);
//         if (!output_path.parent_path().empty()) {
//             std::filesystem::create_directories(
//                 output_path.parent_path());
//         }

//         // Build the instance once. Every solver receives these same arrays.
//         BenchInstance instance;
//         const std::tuple<int, int, double> graph = buildGraphwithPixelsOP(
//                                                             tiling_file,
//                                                             probability_file,
//                                                             instance.costs,
//                                                             instance.dwell_times,
//                                                             instance.tile_probabilities,
//                                                             instance.ranks,
//                                                             instance.member_pixels,
//                                                             instance.pixel_probabilities,
//                                                             zenith_dwell,
//                                                             maximum_slew_rate,
//                                                             maximum_slew_acceleration,
//                                                             settle_time,
//                                                             is_deep_slow,
//                                                             instance.tile_count,
//                                                             source_tile);

//         instance.start_node = std::get<0>(graph);
//         instance.end_node = std::get<1>(graph);
//         instance.padding = std::get<2>(graph);
//         instance.observing_budget = observing_budget;
//         instance.solver_budget = observing_budget + instance.padding;

//         if (instance.costs.size() < 2 ||
//             instance.member_pixels.size() != instance.costs.size()) {
//             throw std::invalid_argument(
//                 "graph and member_pixels have inconsistent node counts");
//         }
//         for (const std::vector<double>& row : instance.costs) {
//             if (row.size() != instance.costs.size()) {
//                 throw std::invalid_argument("cost matrix is not square");
//             }
//         }
//         if (instance.start_node < 0 ||
//             instance.start_node >= static_cast<int>(instance.costs.size()) ||
//             instance.end_node < 0 ||
//             instance.end_node >= static_cast<int>(instance.costs.size()) ||
//             instance.start_node == instance.end_node) {
//             throw std::invalid_argument("graph returned invalid endpoints");
//         }
//         if (!instance.member_pixels[instance.start_node].empty() ||
//             !instance.member_pixels[instance.end_node].empty()) {
//             throw std::invalid_argument(
//                 "dummy start and end nodes must cover no pixels");
//         }
//         if (instance.costs[instance.start_node][instance.end_node] >
//             instance.solver_budget + 1.0e-9) {
//             throw std::runtime_error(
//                 "budget is smaller than the direct start-to-end cost");
//         }

//         const CoverageOracle oracle(
//             instance.pixel_probabilities, instance.member_pixels);
//         if (oracle.nodeCount() !=
//             static_cast<int>(instance.costs.size())) {
//             throw std::invalid_argument(
//                 "member_pixels must have one list per graph node");
//         }

//         double retained_probability = 0.0;
//         for (double probability : instance.pixel_probabilities) {
//             retained_probability += probability;
//         }

//         std::cout << "start " << instance.start_node
//                   << ", end " << instance.end_node
//                   << ", retained tiles " << instance.tile_count
//                   << ", retained pixels "
//                   << instance.pixel_probabilities.size()
//                   << ", retained pixel probability "
//                   << std::fixed << std::setprecision(6)
//                   << retained_probability << "\n";

//         double solver_budget = instance.solver_budget;
//         std::vector<SolverResult> results;

//         // The two deterministic greedy routes also seed GA and SA.
//         results.push_back(runSolver(
//             instance, oracle, "greedy-marginal", [&]() {
//                 return prizeGreedyPathTwoFixedPixelOverlap(
//                     instance.costs,
//                     instance.pixel_probabilities,
//                     instance.member_pixels,
//                     solver_budget,
//                     instance.start_node,
//                     instance.end_node);
//             }));

//         results.push_back(runSolver(
//             instance, oracle, "greedy-ratio", [&]() {
//                 return prizeRatioGreedyPathTwoFixedPixelOverlap(
//                     instance.costs,
//                     instance.pixel_probabilities,
//                     instance.member_pixels,
//                     solver_budget,
//                     instance.start_node,
//                     instance.end_node);
//             }));

//         std::vector<std::vector<int>> genetic_seeds;
//         std::vector<int> annealing_seed;
//         double seed_coverage = -1.0;
//         double seed_cost = std::numeric_limits<double>::infinity();

//         for (std::size_t index = 0;
//              index < 2 && index < results.size(); ++index) {
//             if (results[index].status == "ok") {
//                 genetic_seeds.push_back(results[index].path);
//                 if (resultIsBetter(
//                         results[index], seed_coverage, seed_cost)) {
//                     annealing_seed = results[index].path;
//                     seed_coverage = results[index].coverage;
//                     seed_cost = results[index].padded_cost;
//                 }
//             }
//         }


//         results.push_back(runSolver(
//             instance, oracle, "GCP", [&]() {
//                 return GCP_Coverage(
//                     instance.costs,
//                     instance.pixel_probabilities,
//                     instance.member_pixels,
//                     solver_budget,
//                     instance.start_node,
//                     instance.end_node);
//             }));


//         results.push_back(runSolver(
//             instance, oracle, "simulated-annealing", [&]() {
//                 return simulated_annealing_optimization_pixel_overlap(
//                     instance.costs,
//                     instance.pixel_probabilities,
//                     instance.member_pixels,
//                     solver_budget,
//                     instance.start_node,
//                     instance.end_node,
//                     annealing_seed);
//             }));

//         results.push_back(runSolver(
//             instance, oracle, "genetic", [&]() {
//                 return genetic_optimization_st_pixel_overlap(
//                     instance.costs,
//                     instance.pixel_probabilities,
//                     instance.member_pixels,
//                     solver_budget,
//                     instance.start_node,
//                     instance.end_node,
//                     genetic_seeds);
//             }));

//         // results.push_back(runSolver(
//         //     instance, oracle, "ant-colony", [&]() {
//         //         AntColonyPixelOverlap ant_colony(
//         //             instance.costs,
//         //             instance.pixel_probabilities,
//         //             instance.member_pixels);
//         //         return ant_colony.ant_colony_optimization_pixel_overlap(
//         //             solver_budget,
//         //             instance.start_node,
//         //             instance.end_node);
//         //     }));


//         // results.push_back(runSolver(
//         //     instance, oracle, "RG-QP", [&]() {
//         //         RecursiveGreedyParams parameters;
//         //         const int relevant_arc_bound =
//         //             std::max(1, instance.tile_count + 1);
//         //         parameters.max_depth = std::min(
//         //             50,
//         //             static_cast<int>(std::ceil(
//         //                 1.0 + std::log2(
//         //                     static_cast<double>(relevant_arc_bound)))));
//         //         parameters.reward_levels = 10000;
//         //         parameters.time_quantum = 0.001;
//         //         parameters.time_limit = 0;
//         //         parameters.iterative_deepening = true;
//         //         parameters.candidate_limit = 0;
//         //         parameters.verbose = false;

//         //         const RecursiveGreedyResult rg_result =
//         //             recursiveGreedyOrienteering(
//         //                 instance.costs,
//         //                 instance.pixel_probabilities,
//         //                 instance.member_pixels,
//         //                 instance.start_node,
//         //                 instance.end_node,
//         //                 solver_budget,
//         //                 instance.padding,
//         //                 parameters);
//         //         return rg_result.path;
//         //     }));


//         // // Select the best feasible heuristic route for both ILP warm starts.
//         // std::vector<int> warm_start;
//         // double warm_coverage = -1.0;
//         // double warm_cost = std::numeric_limits<double>::infinity();
//         // for (const SolverResult& result : results) {
//         //     if (resultIsBetter(result, warm_coverage, warm_cost)) {
//         //         warm_start = result.path;
//         //         warm_coverage = result.coverage;
//         //         warm_cost = result.padded_cost;
//         //     }
//         // }

//         // if (warm_start.empty()) {
//         //     std::cout << "warm start: none\n";
//         // } else {
//         //     std::cout << "warm start: " << warm_start.size() - 2
//         //               << " tiles, coverage " << warm_coverage
//         //               << ", padded cost " << warm_cost << "\n";
//         // }

//         // double mtz_bound = std::numeric_limits<double>::quiet_NaN();
//         // SolverResult mtz = runSolver(
//         //     instance, oracle, "ILP MTZ", [&]() {
//         //         const std::pair<std::vector<std::vector<int>>, double> output =
//         //             gurobiPathCoverPixelOverlapMTZ(
//         //                 instance.costs,
//         //                 instance.pixel_probabilities,
//         //                 instance.member_pixels,
//         //                 instance.start_node,
//         //                 instance.end_node,
//         //                 solver_budget,
//         //                 mip_gap,
//         //                 ilp_time_limit,
//         //                 warm_start);
//         //         mtz_bound = output.second;
//         //         return firstPathOf(output);
//         //     });
//         // mtz.bound = mtz_bound;
//         // results.push_back(std::move(mtz));

//         // double continuous_bound =
//         //     std::numeric_limits<double>::quiet_NaN();
//         // SolverResult continuous_time = runSolver(
//         //     instance, oracle, "ILP continuous-time", [&]() {
//         //         const std::pair<std::vector<std::vector<int>>, double> output =
//         //             gurobiPathCoverPixelOverlapContinuousTime(
//         //                 instance.costs,
//         //                 instance.pixel_probabilities,
//         //                 instance.member_pixels,
//         //                 instance.start_node,
//         //                 instance.end_node,
//         //                 solver_budget,
//         //                 mip_gap,
//         //                 ilp_time_limit,
//         //                 warm_start);
//         //         continuous_bound = output.second;
//         //         return firstPathOf(output);
//         //     });
//         // continuous_time.bound = continuous_bound;
//         // results.push_back(std::move(continuous_time));


//         // Modify the instance to retain top N tiles:
//         const int maximum_tiles = 50;
//         const std::tuple<int, int, double> graph2 = buildGraphwithPixelsOP(
//                                                             tiling_file,
//                                                             probability_file,
//                                                             instance.costs,
//                                                             instance.dwell_times,
//                                                             instance.tile_probabilities,
//                                                             instance.ranks,
//                                                             instance.member_pixels,
//                                                             instance.pixel_probabilities,
//                                                             zenith_dwell,
//                                                             maximum_slew_rate,
//                                                             maximum_slew_acceleration,
//                                                             settle_time,
//                                                             is_deep_slow,
//                                                             maximum_tiles,
//                                                             instance.tile_count,
//                                                             source_tile);

//         instance.start_node = std::get<0>(graph2);
//         instance.end_node = std::get<1>(graph2);
//         instance.padding = std::get<2>(graph2);
//         instance.observing_budget = observing_budget;
//         instance.solver_budget = observing_budget + instance.padding;

//         if (instance.costs.size() < 2 ||
//             instance.member_pixels.size() != instance.costs.size()) {
//             throw std::invalid_argument(
//                 "graph and member_pixels have inconsistent node counts");
//         }
//         for (const std::vector<double>& row : instance.costs) {
//             if (row.size() != instance.costs.size()) {
//                 throw std::invalid_argument("cost matrix is not square");
//             }
//         }
//         if (instance.start_node < 0 ||
//             instance.start_node >= static_cast<int>(instance.costs.size()) ||
//             instance.end_node < 0 ||
//             instance.end_node >= static_cast<int>(instance.costs.size()) ||
//             instance.start_node == instance.end_node) {
//             throw std::invalid_argument("graph returned invalid endpoints");
//         }
//         if (!instance.member_pixels[instance.start_node].empty() ||
//             !instance.member_pixels[instance.end_node].empty()) {
//             throw std::invalid_argument(
//                 "dummy start and end nodes must cover no pixels");
//         }
//         if (instance.costs[instance.start_node][instance.end_node] >
//             instance.solver_budget + 1.0e-9) {
//             throw std::runtime_error(
//                 "budget is smaller than the direct start-to-end cost");
//         }

//         // const CoverageOracle oracle(
//         //     instance.pixel_probabilities, instance.member_pixels);
//         // if (oracle.nodeCount() !=
//         //     static_cast<int>(instance.costs.size())) {
//         //     throw std::invalid_argument(
//         //         "member_pixels must have one list per graph node");
//         // }

//         retained_probability = 0.0;
//         for (double probability : instance.pixel_probabilities) {
//             retained_probability += probability;
//         }

//         std::cout << "start " << instance.start_node
//                   << ", end " << instance.end_node
//                   << ", retained tiles " << instance.tile_count
//                   << ", retained pixels "
//                   << instance.pixel_probabilities.size()
//                   << ", retained pixel probability "
//                   << std::fixed << std::setprecision(6)
//                   << retained_probability << "\n";

//         solver_budget = instance.solver_budget;
//         // std::vector<SolverResult> results;


//         results.push_back(runSolver(
//             instance, oracle, "RG-QP", [&]() {
//                 RecursiveGreedyParams parameters;
//                 const int relevant_arc_bound =
//                     std::max(1, instance.tile_count + 1);
//                 parameters.max_depth = std::min(
//                     50,
//                     static_cast<int>(std::ceil(
//                         1.0 + std::log2(
//                             static_cast<double>(relevant_arc_bound)))));
//                 parameters.reward_levels = 100;
//                 parameters.time_quantum = 0.01;
//                 parameters.time_limit = 0;
//                 parameters.iterative_deepening = true;
//                 parameters.candidate_limit = 0;
//                 parameters.verbose = false;

//                 const RecursiveGreedyResult rg_result =
//                     recursiveGreedyOrienteering(
//                         instance.costs,
//                         instance.pixel_probabilities,
//                         instance.member_pixels,
//                         instance.start_node,
//                         instance.end_node,
//                         solver_budget,
//                         instance.padding,
//                         parameters);
//                 return rg_result.path;
//             }));


//         // // Select the best feasible heuristic route for both ILP warm starts.
//         // std::vector<int> warm_start;
//         // double warm_coverage = -1.0;
//         // double warm_cost = std::numeric_limits<double>::infinity();
//         // for (const SolverResult& result : results) {
//         //     if (resultIsBetter(result, warm_coverage, warm_cost)) {
//         //         warm_start = result.path;
//         //         warm_coverage = result.coverage;
//         //         warm_cost = result.padded_cost;
//         //     }
//         // }

//         // if (warm_start.empty()) {
//         //     std::cout << "warm start: none\n";
//         // } else {
//         //     std::cout << "warm start: " << warm_start.size() - 2
//         //               << " tiles, coverage " << warm_coverage
//         //               << ", padded cost " << warm_cost << "\n";
//         // }

//         // double mtz_bound = std::numeric_limits<double>::quiet_NaN();
//         // SolverResult mtz = runSolver(
//         //     instance, oracle, "ILP MTZ", [&]() {
//         //         const std::pair<std::vector<std::vector<int>>, double> output =
//         //             gurobiPathCoverPixelOverlapMTZ(
//         //                 instance.costs,
//         //                 instance.pixel_probabilities,
//         //                 instance.member_pixels,
//         //                 instance.start_node,
//         //                 instance.end_node,
//         //                 solver_budget,
//         //                 mip_gap,
//         //                 ilp_time_limit,
//         //                 warm_start);
//         //         mtz_bound = output.second;
//         //         return firstPathOf(output);
//         //     });
//         // mtz.bound = mtz_bound;
//         // results.push_back(std::move(mtz));

//         std::vector<int> warm_start = {};
//         double continuous_bound =
//             std::numeric_limits<double>::quiet_NaN();
//         SolverResult continuous_time = runSolver(
//             instance, oracle, "ILP continuous-time", [&]() {
//                 const std::pair<std::vector<std::vector<int>>, double> output =
//                     gurobiPathCoverPixelOverlapContinuousTime(
//                         instance.costs,
//                         instance.pixel_probabilities,
//                         instance.member_pixels,
//                         instance.start_node,
//                         instance.end_node,
//                         solver_budget,
//                         mip_gap,
//                         ilp_time_limit,
//                         warm_start);
//                 continuous_bound = output.second;
//                 return firstPathOf(output);
//             });
//         continuous_time.bound = continuous_bound;
//         results.push_back(std::move(continuous_time));


//         printTable(instance, results);
//         appendCsv(
//             output_file,
//             instance,
//             tiling_file,
//             probability_file,
//             maximum_slew_rate,
//             results);

//         const double total_seconds = std::chrono::duration<double>(
//             std::chrono::high_resolution_clock::now() - total_started).count();
//         std::cout << "\nTotal running time: " << std::setprecision(2)
//                   << total_seconds << " seconds\n";
//     } catch (const std::exception& error) {
//         std::cerr << "Pixel-overlap benchmark failed: "
//                   << error.what() << "\n";
//         return 1;
//     }

//     return 0;
// }





// #include "BuildGraph.h"
// #include "ILP_gurobi_path_cover_op.h"


// #include <chrono>
// #include <cmath>
// #include <filesystem>
// #include <fstream>
// #include <iostream>
// #include <stdexcept>
// #include <string>
// #include <tuple>
// #include <utility>
// #include <vector>


// double coveredProbabilityUnionSTPixelOverlap(
//         const std::vector<int>& path_tiles,
//         const std::vector<std::vector<int>>& member_pixels,
//         const std::vector<double>& pixel_probs) {

//     const int nTiles  = static_cast<int>(member_pixels.size());
//     const int nPixels = static_cast<int>(pixel_probs.size());
//     if (nPixels == 0) return 0.0;

//     std::vector<char> seen(nPixels, 0);
//     double total = 0.0;

//     for (int tile : path_tiles) {
//         if (tile < 0 || tile >= nTiles) {
//             throw std::runtime_error(
//                 "coveredProbabilityUnionSTPixelOverlap: "
//                 "tile index out of range.");
//         }
//         for (int p : member_pixels[tile]) {
//             if (p < 0 || p >= nPixels) {
//                 throw std::runtime_error(
//                     "coveredProbabilityUnionSTPixelOverlap: "
//                     "pixel index out of range.");
//             }
//             if (!seen[p]) {
//                 seen[p] = 1;
//                 total += pixel_probs[p];
//             }
//         }
//     }

//     return total;
// }


// // costs is the all-in edge matrix returned by buildGraph. Do not add dwell
// // time again: each selected tile's dwell and the end padding are already on
// // the corresponding edges.
// double printPathSTPixelOverlap(
//         const std::vector<std::vector<double>>& costs,
//         const std::vector<std::vector<int>>& member_pixels,
//         const std::vector<double>& pixel_probs,
//         const std::vector<int>& rank,
//         const std::vector<int>& best_path,
//         int start_idx,
//         int end_idx) {

//     if (best_path.empty()) {
//         std::cout << "No incumbent path was returned.\n";
//         return 0.0;
//     }
//     if (best_path.front() != start_idx || best_path.back() != end_idx) {
//         throw std::runtime_error(
//             "printPathSTPixelOverlap: invalid path endpoints.");
//     }

//     double sum_probability = coveredProbabilityUnionSTPixelOverlap(
//         best_path, member_pixels, pixel_probs);

//     double total_cost = 0.0;
//     for (std::size_t i = 1; i < best_path.size(); ++i) {
//         total_cost += costs[best_path[i - 1]][best_path[i]];
//     }

//     std::cout << "Path(tile rank): ";
//     for (std::size_t i = 1; i + 1 < best_path.size(); ++i) {
//         const int tile = best_path[i];
//         std::cout << (rank.empty() ? tile : rank.at(tile)) << " ";
//     }
//     std::cout << std::endl;

//     const std::size_t num_selected_tiles =
//         best_path.size() >= 2 ? best_path.size() - 2 : 0;
//     std::cout << "Num Tiles in Path: " << num_selected_tiles << std::endl;
//     std::cout << "Sum Probability: " << sum_probability << std::endl;
//     std::cout << "Total Cost: " << total_cost << std::endl;

//     return sum_probability;
// }

// void saveResultSTPixelOverlap(
//         const std::string& filename,
//         const std::string& method,
//         const std::string& tiling_data,
//         const std::string& likelihood_data,
//         double budget,
//         double slew_rate,
//         const std::vector<std::vector<double>>& costs,
//         const std::vector<std::vector<int>>& member_pixels,
//         const std::vector<double>& pixel_probs,
//         const std::vector<int>& rank,
//         const std::pair<std::vector<std::vector<int>>, double>& result,
//         int start_idx,
//         int end_idx,
//         double elapsed_time) {

//     bool file_exists = std::filesystem::exists(filename);
//     std::ofstream outfile(filename, std::ios::app);

//     if (!outfile.is_open()) {
//         std::cerr << "Error: Could not open results file: "
//                   << filename << std::endl;
//         return;
//     }

//     if (!file_exists) {
//         outfile
//             << "Method,TilingData,LikelihoodData,Budget,SlewRate,NumTiles,"
//                "SumProb,BestBound,TotalCost,TimeSec,Path\n";
//     }

//     if (result.first.size() > 1) {
//         throw std::runtime_error(
//             "saveResultSTPixelOverlap: the single-telescope solver "
//             "returned more than one path.");
//     }

//     const std::vector<int> best_path =
//         result.first.empty() ? std::vector<int>{} : result.first.front();
//     double sum_probability = 0.0;
//     double total_cost = 0.0;
//     std::size_t num_selected_tiles = 0;

//     if (!best_path.empty()) {
//         if (best_path.front() != start_idx || best_path.back() != end_idx) {
//             throw std::runtime_error(
//                 "saveResultSTPixelOverlap: invalid path endpoints.");
//         }
//         sum_probability = coveredProbabilityUnionSTPixelOverlap(
//             best_path, member_pixels, pixel_probs);
//         for (std::size_t i = 1; i < best_path.size(); ++i) {
//             total_cost += costs[best_path[i - 1]][best_path[i]];
//         }
//         num_selected_tiles =
//             best_path.size() >= 2 ? best_path.size() - 2 : 0;
//     }

//     outfile << method << ","
//             << std::filesystem::path(tiling_data).filename().string() << ","
//             << std::filesystem::path(likelihood_data).filename().string()
//             << ","
//             << budget << ","
//             << slew_rate << ","
//             << num_selected_tiles << ","
//             << sum_probability << ","
//             << result.second << ","
//             << total_cost << ","
//             << elapsed_time << ",";

//     // Save only real tile ranks; omit the dummy start and end nodes.
//     for (std::size_t i = 1; i + 1 < best_path.size(); ++i) {
//         const int tile = best_path[i];
//         outfile << (rank.empty() ? tile : rank.at(tile));
//         if (i + 2 < best_path.size()) {
//             outfile << " ";
//         }
//     }
//     outfile << "\n";
// }


// void testPathsCoverSTPixelOverlapConfigured(
//                                             std::string tilefile,
//                                             std::string mapfile,
//                                             std::string out_file,
//                                             std::string overlap,
//                                             double budget,
//                                             double dwell_time,
//                                             double w_max,
//                                             double w_acc,
//                                             double settle_time,
//                                             bool is_deepslow,
//                                             int init_pos_idx,
//                                             double accu_thr,
//                                             double time_limit) {

//     std::vector<std::vector<double>> costs;
//     std::vector<double> tile_probs;
//     std::vector<int> ranks;
//     std::vector<double> dwell_times;
//     std::vector<std::vector<int>> member_pixels;
//     std::vector<double> pixel_probs;

//     int num_tiles = 0;

//     const auto [start_idx, end_idx, padding] = buildGraphwithPixelsOP(
//         tilefile, mapfile,
//         costs, dwell_times, tile_probs, ranks,
//         member_pixels, pixel_probs,
//         dwell_time, w_max, w_acc, settle_time,
//         is_deepslow, num_tiles, init_pos_idx);

//     if (costs[start_idx][end_idx] > budget + 1e-9) {
//         throw std::runtime_error(
//             "testPathsCoverSTPixelOverlapConfigured: budget is smaller "
//             "than the direct start-to-end padding cost.");
//     }

//     std::cout << "start: " << start_idx
//               << ", end: " << end_idx
//               << ", retained tiles: " << num_tiles
//               << ", padding: " << padding << std::endl;

//     std::cout
//         << "*********ILP Path Cover (Continuous Time, Pixel Overlap)*********" << std::endl;

//     const auto solve_start = std::chrono::high_resolution_clock::now();

//     const std::pair<std::vector<std::vector<int>>, double> result =
//         gurobiPathCoverPixelOverlapContinuousTime(costs, pixel_probs, member_pixels,
//                                                 start_idx, end_idx, budget,
//                                                 accu_thr, time_limit);


//     // result = gurobiPathCoverPixelOverlapMTZ(
//     //               costs, pixel_probs, member_pixels,
//     //               start_idx, end_idx, budget,
//     //               accu_thr, time_limit);

//     const auto solve_end = std::chrono::high_resolution_clock::now();

//     const std::chrono::duration<double> elapsed_seconds = solve_end - solve_start;

//     if (result.first.size() > 1) {
//         throw std::runtime_error(
//             "The single-telescope ILP returned more than one path.");
//     }

//     const std::vector<int> best_path = result.first.empty() ? std::vector<int>{} : result.first.front();

//     std::cout << "Best Bound: " << result.second << std::endl;

//     printPathSTPixelOverlap( costs, member_pixels, pixel_probs, ranks, best_path, start_idx, end_idx);

//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

//     saveResultSTPixelOverlap(out_file,
//         overlap + "path+cover+continuous-time",
//         tilefile, mapfile,
//         budget, w_max,
//         costs, member_pixels, pixel_probs,
//         ranks, result,
//         start_idx, end_idx,
//         elapsed_seconds.count());
// }


// int main(int argc, char** argv) {
//     const auto total_start =
//         std::chrono::high_resolution_clock::now();

//     std::string out_file = "../results/out.csv";
//     std::string tilefile =
//         "../data/tilings_small/6.9x6.9_tiling.csv";
//     std::string mapfile =
//         "../data/test_maps/GW200316_215756_50.txt";

//     double budget = 60.0;
//     double w_max = 10.0;
//     double w_acc = 10.0;
//     double dwell_zenith = 2.0;
//     bool is_deepslow = false;

//     // Fixed solver/graph parameters.
//     const double settle_time = 0.0;
//     const int source_tile = -2;
//     const double mip_gap = 0.001;
//     const double time_limit = 28800.0;

//     if (argc > 9) {
//         std::cerr
//             << "Usage: " << argv[0]
//             << " [MAP_FILE TILE_FILE BUDGET W_MAX W_ACC "
//                "DWELL_ZENITH IS_DEEPSLOW OUTPUT_CSV]\n";
//         return 2;
//     }

//     try {
//         if (argc > 1) mapfile = argv[1];
//         if (argc > 2) tilefile = argv[2];
//         if (argc > 3) budget = std::stod(argv[3]);
//         if (argc > 4) w_max = std::stod(argv[4]);
//         if (argc > 5) w_acc = std::stod(argv[5]);
//         if (argc > 6) dwell_zenith = std::stod(argv[6]);
//         if (argc > 7) {
//             is_deepslow = std::stoi(argv[7]) != 0;
//         }
//         if (argc > 8) out_file = argv[8];

//         if (!std::isfinite(budget) || budget < 0.0 ||
//             !std::isfinite(w_max) || w_max <= 0.0 ||
//             !std::isfinite(w_acc) || w_acc <= 0.0 ||
//             !std::isfinite(dwell_zenith) || dwell_zenith < 0.0) {
//             throw std::invalid_argument(
//                 "Invalid budget, motion, or dwell parameter.");
//         }

//         const std::filesystem::path output_path(out_file);
//         if (!output_path.parent_path().empty()) {
//             std::filesystem::create_directories(
//                 output_path.parent_path());
//         }

//         testPathsCoverSTPixelOverlapConfigured(
//             tilefile,
//             mapfile,
//             out_file,
//             "overlap+",
//             budget,
//             dwell_zenith,
//             w_max,
//             w_acc,
//             settle_time,
//             is_deepslow,
//             source_tile,
//             mip_gap,
//             time_limit);

//         const auto total_end =
//             std::chrono::high_resolution_clock::now();

//         const std::chrono::duration<double> total_elapsed =
//             total_end - total_start;

//         std::cout << "Total running time: "
//                   << total_elapsed.count()
//                   << " seconds" << std::endl;

//     } catch (const std::exception& error) {
//         std::cerr << "Pixel-overlap path-cover failed: "
//                   << error.what() << std::endl;
//         return 1;
//     }

//     return 0;
// }
