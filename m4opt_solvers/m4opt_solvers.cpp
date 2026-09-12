#include "m4opt_solvers.h"

#include "AntColonyPixelOverlap.h"
#include "GCP_overlap.h"
#include "GeneticPixelOverlap.h"
#include "GreedyPixelOverlap.h"
#include "ILP_gurobi_path_cover_op.h"
#include "SimulatedAnnealingPixelOverlap.h"
#include "recursive_greedy_op.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

static const double TOLERANCE = 1.0e-9;

// Check the physical slew matrix and scalar time inputs. The caller has
// already confirmed the triangle inequality, so this deliberately avoids an
// expensive all-triples metric test.
static void validatePhysicalGraph(
    const std::vector<std::vector<double>>& slew_seconds,
    double physical_budget_seconds) {

    if (!std::isfinite(physical_budget_seconds) ||
        physical_budget_seconds < 0.0) {
        throw std::invalid_argument(
            "physical budget must be finite and nonnegative");
    }

    const int tile_count = static_cast<int>(slew_seconds.size());
    if (tile_count == 0) {
        throw std::invalid_argument("there must be at least one tile");
    }

    // Check every row before indexing slew_seconds[j][i] below.
    for (int i = 0; i < tile_count; ++i) {
        if (static_cast<int>(slew_seconds[i].size()) != tile_count) {
            throw std::invalid_argument("slew matrix must be square");
        }
    }

    for (int i = 0; i < tile_count; ++i) {
        for (int j = 0; j < tile_count; ++j) {
            const double value = slew_seconds[i][j];
            if (!std::isfinite(value) || value < 0.0) {
                throw std::invalid_argument(
                    "slew times must be finite and nonnegative");
            }
            if (i == j && std::abs(value) > TOLERANCE) {
                throw std::invalid_argument(
                    "the diagonal of the slew matrix must be zero");
            }
            if (std::abs(value - slew_seconds[j][i]) > TOLERANCE) {
                throw std::invalid_argument(
                    "this simple bridge requires symmetric slew times");
            }
        }
    }
}

// Sort and deduplicate every footprint. Pixel bounds are checked later when
// the pixel-probability vector is available.
static std::vector<std::vector<int>> normalizeFootprints(
    const std::vector<std::vector<int>>& member_pixels) {

    std::vector<std::vector<int>> normalized = member_pixels;
    for (std::vector<int>& pixels : normalized) {
        std::sort(pixels.begin(), pixels.end());
        pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
    }
    return normalized;
}

// Shared half-dwell transformation used by both P1 and the P2 mode graph.
// M4OPT does not specify a fixed initial pointing, so the first and last
// physical pointings are free. Equal padding is placed at both dummy ends:
//
//   transformed path cost = physical route duration + 2 * padding.
//
// This keeps the complete graph symmetric and metric.
static OrienteeringGraph addDummyEndpoints(
    const std::vector<std::vector<double>>& physical_slew_seconds,
    const std::vector<std::vector<int>>& mode_pixels,
    const std::vector<int>& tile_of_mode,
    const std::vector<double>& exposure_of_mode,
    double physical_budget_seconds) {

    validatePhysicalGraph(physical_slew_seconds, physical_budget_seconds);

    if (mode_pixels.size() != tile_of_mode.size() ||
        mode_pixels.size() != exposure_of_mode.size()) {
        throw std::invalid_argument("inconsistent mode arrays");
    }

    const int tile_count =
        static_cast<int>(physical_slew_seconds.size());
    const int mode_count = static_cast<int>(mode_pixels.size());
    const int start_node = mode_count;
    const int end_node = mode_count + 1;

    double maximum_slew = 0.0;
    for (const std::vector<double>& row : physical_slew_seconds) {
        for (double value : row) {
            maximum_slew = std::max(maximum_slew, value);
        }
    }
    const double end_padding = std::ceil(maximum_slew);

    OrienteeringGraph graph;
    graph.start_node = start_node;
    graph.end_node = end_node;
    graph.physical_tile_count = tile_count;
    graph.physical_budget = physical_budget_seconds;
    graph.padding = 2.0 * end_padding;
    graph.solver_budget = physical_budget_seconds + graph.padding;

    graph.member_pixels = normalizeFootprints(mode_pixels);
    graph.member_pixels.resize(mode_count + 2);

    graph.tile_of_node = tile_of_mode;
    graph.tile_of_node.resize(mode_count + 2, -1);

    graph.exposure_of_node = exposure_of_mode;
    graph.exposure_of_node.resize(mode_count + 2, 0.0);

    graph.costs.assign(
        mode_count + 2, std::vector<double>(mode_count + 2, 0.0));

    for (int first = 0; first < mode_count; ++first) {
        const int first_tile = graph.tile_of_node[first];
        if (first_tile < 0 || first_tile >= tile_count) {
            throw std::out_of_range("invalid physical tile index");
        }
        if (!std::isfinite(graph.exposure_of_node[first]) ||
            graph.exposure_of_node[first] < 0.0) {
            throw std::invalid_argument(
                "exposure times must be finite and nonnegative");
        }

        for (int second = first + 1; second < mode_count; ++second) {
            const int second_tile = graph.tile_of_node[second];
            if (second_tile < 0 || second_tile >= tile_count) {
                throw std::out_of_range("invalid physical tile index");
            }

            // Different exposure modes of one physical field have zero slew.
            const double slew = first_tile == second_tile
                ? 0.0
                : physical_slew_seconds[first_tile][second_tile];
            const double cost =
                slew +
                0.5 * graph.exposure_of_node[first] +
                0.5 * graph.exposure_of_node[second];

            graph.costs[first][second] = cost;
            graph.costs[second][first] = cost;
        }
    }

    for (int node = 0; node < mode_count; ++node) {
        const double dummy_cost = end_padding + 0.5 * graph.exposure_of_node[node];
        graph.costs[start_node][node] = dummy_cost;
        graph.costs[node][start_node] = dummy_cost;
        graph.costs[node][end_node] = dummy_cost;
        graph.costs[end_node][node] = dummy_cost;
    }

    graph.costs[start_node][end_node] = graph.padding;
    graph.costs[end_node][start_node] = graph.padding;
    return graph;
}

OrienteeringGraph buildP1Graph(
                    const std::vector<std::vector<double>>& slew_seconds,
                    const std::vector<std::vector<int>>& member_pixels,
                    const std::vector<double>& dwell_seconds,
                    double physical_budget_seconds) {

    const int tile_count = static_cast<int>(slew_seconds.size());
    if (static_cast<int>(member_pixels.size()) != tile_count ||
        static_cast<int>(dwell_seconds.size()) != tile_count) {
        throw std::invalid_argument(
            "P1 needs one footprint and one dwell time per tile");
    }

    std::vector<int> tile_of_node(tile_count);
    for (int tile = 0; tile < tile_count; ++tile) {
        tile_of_node[tile] = tile;
    }

    return addDummyEndpoints(
        slew_seconds,
        member_pixels,
        tile_of_node,
        dwell_seconds,
        physical_budget_seconds);
}

OrienteeringGraph buildP2Graph(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& required_exposure_seconds,
    double minimum_exposure_seconds,
    double maximum_exposure_seconds,
    double physical_budget_seconds) {

    const int tile_count = static_cast<int>(slew_seconds.size());
    const int pixel_count = static_cast<int>(pixel_probabilities.size());

    if (static_cast<int>(member_pixels.size()) != tile_count) {
        throw std::invalid_argument("P2 needs one footprint per tile");
    }
    if (static_cast<int>(required_exposure_seconds.size()) != pixel_count) {
        throw std::invalid_argument(
            "P2 needs one exposure threshold per compact pixel");
    }
    if (!std::isfinite(minimum_exposure_seconds) ||
        !std::isfinite(maximum_exposure_seconds) ||
        minimum_exposure_seconds < 0.0 ||
        maximum_exposure_seconds < minimum_exposure_seconds) {
        throw std::invalid_argument("invalid P2 exposure limits");
    }

    for (int pixel = 0; pixel < pixel_count; ++pixel) {
        if (!std::isfinite(pixel_probabilities[pixel]) ||
            pixel_probabilities[pixel] < 0.0) {
            throw std::invalid_argument(
                "pixel probabilities must be finite and nonnegative");
        }
        const double threshold = required_exposure_seconds[pixel];
        if (std::isnan(threshold) || threshold < 0.0) {
            throw std::invalid_argument(
                "pixel exposure thresholds must be nonnegative or infinity");
        }
    }

    const std::vector<std::vector<int>> footprints =
        normalizeFootprints(member_pixels);

    std::vector<std::vector<int>> mode_pixels;
    std::vector<int> tile_of_mode;
    std::vector<double> exposure_of_mode;

    for (int tile = 0; tile < tile_count; ++tile) {
        std::vector<double> levels{minimum_exposure_seconds};

        for (int pixel : footprints[tile]) {
            if (pixel < 0 || pixel >= pixel_count) {
                throw std::out_of_range(
                    "a tile contains an invalid compact pixel index");
            }

            const double threshold = required_exposure_seconds[pixel];
            if (pixel_probabilities[pixel] > 0.0 &&
                std::isfinite(threshold) &&
                threshold > minimum_exposure_seconds &&
                threshold <= maximum_exposure_seconds) {
                levels.push_back(threshold);
            }
        }

        std::sort(levels.begin(), levels.end());
        levels.erase(std::unique(levels.begin(), levels.end()), levels.end());

        // A higher exposure with the same detected-pixel set is dominated,
        // so only keep levels at which positive-probability coverage changes.
        std::vector<int> previous_pixels;
        for (double exposure : levels) {
            std::vector<int> detected_pixels;
            for (int pixel : footprints[tile]) {
                if (pixel_probabilities[pixel] > 0.0 &&
                    required_exposure_seconds[pixel] <= exposure) {
                    detected_pixels.push_back(pixel);
                }
            }

            if (detected_pixels.empty() ||
                detected_pixels == previous_pixels) {
                continue;
            }

            mode_pixels.push_back(detected_pixels);
            tile_of_mode.push_back(tile);
            exposure_of_mode.push_back(exposure);
            previous_pixels = detected_pixels;
        }
    }

    return addDummyEndpoints(
        slew_seconds,
        mode_pixels,
        tile_of_mode,
        exposure_of_mode,
        physical_budget_seconds);
}

// The mode-expanded P2 graph may contain several selected modes of one tile.
// Keep only that tile's largest selected exposure. Because its cumulative
// footprint contains all smaller-mode footprints and the graph is metric,
// this keeps the reward and cannot increase the route cost.
static std::vector<int> removeRepeatedP2Modes(
    const OrienteeringGraph& graph,
    const std::vector<int>& path) {

    std::vector<int> largest_mode(graph.physical_tile_count, -1);
    for (std::size_t position = 1; position + 1 < path.size(); ++position) {
        const int node = path[position];
        if (node < 0 || node >= graph.start_node) {
            throw std::runtime_error("solver returned an invalid mode node");
        }
        const int tile = graph.tile_of_node[node];
        const int old_node = largest_mode[tile];
        if (old_node < 0 ||
            graph.exposure_of_node[node] >
                graph.exposure_of_node[old_node]) {
            largest_mode[tile] = node;
        }
    }

    std::vector<int> clean_path{graph.start_node};
    std::vector<char> emitted(graph.physical_tile_count, 0);
    for (std::size_t position = 1; position + 1 < path.size(); ++position) {
        const int node = path[position];
        const int tile = graph.tile_of_node[node];
        if (node == largest_mode[tile] && emitted[tile] == 0) {
            clean_path.push_back(node);
            emitted[tile] = 1;
        }
    }
    clean_path.push_back(graph.end_node);
    return clean_path;
}

static void validatePixelIndices(
    const std::vector<double>& pixel_probabilities,
    const OrienteeringGraph& graph) {

    for (double probability : pixel_probabilities) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "pixel probabilities must be finite and nonnegative");
        }
    }

    for (const std::vector<int>& pixels : graph.member_pixels) {
        for (int pixel : pixels) {
            if (pixel < 0 ||
                pixel >= static_cast<int>(pixel_probabilities.size())) {
                throw std::out_of_range(
                    "a tile contains an invalid compact pixel index");
            }
        }
    }
}

// Convert names such as "simulated-annealing" to "simulated_annealing".
static std::string normalizeAlgorithmName(std::string name) {
    for (char& character : name) {
        if (character == '-' || character == ' ') {
            character = '_';
        } else {
            character = static_cast<char>(
                std::tolower(static_cast<unsigned char>(character)));
        }
    }
    return name;
}

std::vector<std::string> availableAlgorithms() {
    return {
        "greedy_marginal",
        "greedy_ratio",
        "gcp",
        "simulated_annealing",
        "genetic",
        "ant_colony",
        "rg_qp",
        "ilp_mtz",
        "ilp_continuous",
    };
}

static double graphPathCost(
    const OrienteeringGraph& graph,
    const std::vector<int>& path) {

    double total = 0.0;
    for (std::size_t i = 1; i < path.size(); ++i) {
        total += graph.costs[path[i - 1]][path[i]];
    }
    return total;
}

static double graphPathCoverage(
    const OrienteeringGraph& graph,
    const std::vector<double>& pixel_probabilities,
    const std::vector<int>& path) {

    std::vector<char> covered(pixel_probabilities.size(), 0);
    for (int node : path) {
        for (int pixel : graph.member_pixels[node]) {
            covered[pixel] = 1;
        }
    }

    double total = 0.0;
    for (std::size_t pixel = 0; pixel < covered.size(); ++pixel) {
        if (covered[pixel] != 0) {
            total += pixel_probabilities[pixel];
        }
    }
    return total;
}

// Build the same two greedy seeds used by the shared main function.
static std::vector<std::vector<int>> makeGreedySeeds(
    const OrienteeringGraph& graph,
    const std::vector<double>& pixel_probabilities) {

    std::vector<std::vector<int>> seeds;
    seeds.push_back(prizeGreedyPathTwoFixedPixelOverlap(
        graph.costs,
        pixel_probabilities,
        graph.member_pixels,
        graph.solver_budget,
        graph.start_node,
        graph.end_node));
    seeds.push_back(prizeRatioGreedyPathTwoFixedPixelOverlap(
        graph.costs,
        pixel_probabilities,
        graph.member_pixels,
        graph.solver_budget,
        graph.start_node,
        graph.end_node));
    return seeds;
}

static std::vector<int> bestGreedySeed(
    const OrienteeringGraph& graph,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& seeds) {

    std::vector<int> best;
    double best_reward = -1.0;
    double best_cost = std::numeric_limits<double>::infinity();

    for (const std::vector<int>& path : seeds) {
        if (path.empty()) {
            continue;
        }
        const double reward = graphPathCoverage(
            graph, pixel_probabilities, path);
        const double cost = graphPathCost(graph, path);
        if (reward > best_reward + TOLERANCE ||
            (std::abs(reward - best_reward) <= TOLERANCE &&
             cost < best_cost)) {
            best = path;
            best_reward = reward;
            best_cost = cost;
        }
    }
    return best;
}

struct SolverOutput {
    std::string algorithm;
    std::vector<int> path;
    double best_bound = std::numeric_limits<double>::quiet_NaN();
    int completed_depth = 0;
    bool timed_out = false;
};

static RecursiveGreedyParams makeParameters(
    int physical_tile_count,
    int max_depth,
    long long reward_levels,
    double time_quantum,
    double time_limit_seconds) {

    RecursiveGreedyParams parameters;
    if (max_depth <= 0) {
        max_depth = static_cast<int>(std::ceil(
            1.0 + std::log2(
                static_cast<double>(physical_tile_count + 1))));
    }

    parameters.max_depth = std::max(1, std::min(30, max_depth));
    parameters.reward_levels = reward_levels;
    parameters.time_quantum = time_quantum;
    parameters.time_limit = time_limit_seconds;
    parameters.iterative_deepening = true;
    parameters.candidate_limit = 0;
    parameters.verbose = false;

    // The bridge constructs a metric graph from the caller-certified metric
    // slew matrix. This skips RG-QP's otherwise cubic diagnostic check.
    // parameters.trust_metric_construction = true;
    return parameters;
}

// This is the only dispatch function. It directly mirrors the calls in the
// shared benchmark main, including greedy seeds for SA, GA, and the ILPs.
static SolverOutput runAlgorithm(
    const OrienteeringGraph& graph,
    const std::vector<double>& pixel_probabilities,
    const std::string& requested_algorithm,
    int max_depth,
    long long reward_levels,
    double time_quantum,
    double time_limit_seconds,
    int thread_count,
    double mip_gap) {

    SolverOutput output;
    output.algorithm = normalizeAlgorithmName(requested_algorithm);
    if (output.algorithm == "ilp_continuous_time") {
        output.algorithm = "ilp_continuous";
    }

    if (thread_count < 0) {
        throw std::invalid_argument("thread_count must be nonnegative");
    }
#ifdef _OPENMP
    if (thread_count > 0) {
        omp_set_num_threads(thread_count);
    }
#else
    (void)thread_count;
#endif

    if (output.algorithm == "greedy_marginal") {
        output.path = prizeGreedyPathTwoFixedPixelOverlap(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.solver_budget,
            graph.start_node,
            graph.end_node);
    } else if (output.algorithm == "greedy_ratio") {
        output.path = prizeRatioGreedyPathTwoFixedPixelOverlap(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.solver_budget,
            graph.start_node,
            graph.end_node);
    } else if (output.algorithm == "gcp") {
        output.path = GCP_Coverage(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.solver_budget,
            graph.start_node,
            graph.end_node);
    } else if (output.algorithm == "simulated_annealing") {
        const std::vector<std::vector<int>> seeds = makeGreedySeeds(
            graph, pixel_probabilities);
        output.path = simulated_annealing_optimization_pixel_overlap(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.solver_budget,
            graph.start_node,
            graph.end_node,
            bestGreedySeed(graph, pixel_probabilities, seeds));
    } else if (output.algorithm == "genetic") {
        output.path = genetic_optimization_st_pixel_overlap(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.solver_budget,
            graph.start_node,
            graph.end_node,
            makeGreedySeeds(graph, pixel_probabilities));
    } else if (output.algorithm == "ant_colony") {
        AntColonyPixelOverlap ant_colony(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels);
        output.path = ant_colony.ant_colony_optimization_pixel_overlap(
            graph.solver_budget,
            graph.start_node,
            graph.end_node);
    } else if (output.algorithm == "rg_qp") {
        const RecursiveGreedyParams parameters = makeParameters(
            graph.physical_tile_count,
            max_depth,
            reward_levels,
            time_quantum,
            time_limit_seconds);
        const RecursiveGreedyResult result = recursiveGreedyOrienteering(
            graph.costs,
            pixel_probabilities,
            graph.member_pixels,
            graph.start_node,
            graph.end_node,
            graph.solver_budget,
            graph.padding,
            parameters);
        output.path = result.path;
        output.completed_depth = result.completed_depth;
        output.timed_out = result.timed_out;
    } else if (output.algorithm == "ilp_mtz" ||
               output.algorithm == "ilp_continuous") {
        const std::vector<std::vector<int>> seeds = makeGreedySeeds(
            graph, pixel_probabilities);
        const std::vector<int> warm_start = bestGreedySeed(
            graph, pixel_probabilities, seeds);
        const double ilp_time_limit =
            time_limit_seconds > 0.0 ? time_limit_seconds : -1.0;

        std::pair<std::vector<std::vector<int>>, double> result;
        if (output.algorithm == "ilp_mtz") {
            result = gurobiPathCoverPixelOverlapMTZ(
                graph.costs,
                pixel_probabilities,
                graph.member_pixels,
                graph.start_node,
                graph.end_node,
                graph.solver_budget,
                mip_gap,
                ilp_time_limit,
                warm_start);
        } else {
            result = gurobiPathCoverPixelOverlapContinuousTime(
                graph.costs,
                pixel_probabilities,
                graph.member_pixels,
                graph.start_node,
                graph.end_node,
                graph.solver_budget,
                mip_gap,
                ilp_time_limit,
                warm_start);
        }

        if (result.first.size() > 1) {
            throw std::runtime_error(
                "a single-telescope ILP returned more than one path");
        }
        if (!result.first.empty()) {
            output.path = result.first.front();
        }
        output.best_bound = result.second;
    } else {
        throw std::invalid_argument(
            "unknown algorithm: " + requested_algorithm);
    }

    return output;
}

// Convert the graph path back to an M4OPT tile sequence, recompute physical
// timing, and score the union of covered pixels.
static M4OptResult makeM4OptResult(
    const OrienteeringGraph& graph,
    const std::vector<std::vector<double>>& physical_slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<int>& original_path,
    bool p2_mode_graph,
    const SolverOutput& solver_output,
    double runtime_seconds) {

    M4OptResult result;
    result.algorithm = solver_output.algorithm;
    result.runtime_seconds = runtime_seconds;
    result.best_bound = solver_output.best_bound;
    result.completed_depth = solver_output.completed_depth;
    result.timed_out = solver_output.timed_out;

    if (original_path.empty()) {
        return result;
    }
    if (original_path.size() < 2 ||
        original_path.front() != graph.start_node ||
        original_path.back() != graph.end_node) {
        throw std::runtime_error("solver did not return a complete s-t path");
    }

    const std::vector<int> path = p2_mode_graph
        ? removeRepeatedP2Modes(graph, original_path)
        : original_path;

    result.has_path = true;

    std::vector<char> used_tile(graph.physical_tile_count, 0);
    std::vector<char> covered(pixel_probabilities.size(), 0);

    double current_time = 0.0;
    int previous_tile = -1;

    for (std::size_t position = 1; position + 1 < path.size(); ++position) {
        const int node = path[position];
        const int tile = graph.tile_of_node[node];

        if (tile < 0 || tile >= graph.physical_tile_count) {
            throw std::runtime_error("solver returned an invalid graph node");
        }
        if (used_tile[tile] != 0) {
            throw std::runtime_error("result contains a repeated physical tile");
        }
        used_tile[tile] = 1;

        if (previous_tile >= 0) {
            current_time += physical_slew_seconds[previous_tile][tile];
        }

        result.tile_indices.push_back(tile);
        result.exposure_seconds.push_back(graph.exposure_of_node[node]);
        result.start_seconds.push_back(current_time);

        current_time += graph.exposure_of_node[node];
        previous_tile = tile;

        for (int pixel : graph.member_pixels[node]) {
            covered[pixel] = 1;
        }
    }

    result.duration_seconds = current_time;
    if (result.duration_seconds > graph.physical_budget + TOLERANCE) {
        throw std::runtime_error(
            "the converted physical route exceeds the observing budget");
    }

    for (std::size_t pixel = 0; pixel < covered.size(); ++pixel) {
        if (covered[pixel] != 0) {
            result.coverage += pixel_probabilities[pixel];
        }
    }
    return result;
}

static M4OptResult solveGraph(
    const OrienteeringGraph& graph,
    const std::vector<std::vector<double>>& physical_slew_seconds,
    const std::vector<double>& pixel_probabilities,
    bool p2_mode_graph,
    const std::string& algorithm,
    int max_depth,
    long long reward_levels,
    double time_quantum,
    double time_limit_seconds,
    int thread_count,
    double mip_gap) {

    validatePixelIndices(pixel_probabilities, graph);

    const std::chrono::steady_clock::time_point begin =
        std::chrono::steady_clock::now();

    const SolverOutput solver_output = runAlgorithm(
        graph,
        pixel_probabilities,
        algorithm,
        max_depth,
        reward_levels,
        time_quantum,
        time_limit_seconds,
        thread_count,
        mip_gap);

    const double runtime_seconds = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - begin).count();

    return makeM4OptResult(
        graph,
        physical_slew_seconds,
        pixel_probabilities,
        solver_output.path,
        p2_mode_graph,
        solver_output,
        runtime_seconds);
}

M4OptResult solveP1(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& dwell_seconds,
    double physical_budget_seconds,
    const std::string& algorithm,
    int max_depth,
    long long reward_levels,
    double time_quantum,
    double time_limit_seconds,
    int thread_count,
    double mip_gap) {

    const OrienteeringGraph graph = buildP1Graph(
        slew_seconds,
        member_pixels,
        dwell_seconds,
        physical_budget_seconds);

    return solveGraph(
        graph,
        slew_seconds,
        pixel_probabilities,
        false,
        algorithm,
        max_depth,
        reward_levels,
        time_quantum,
        time_limit_seconds,
        thread_count,
        mip_gap);
}

M4OptResult solveP2(
    const std::vector<std::vector<double>>& slew_seconds,
    const std::vector<double>& pixel_probabilities,
    const std::vector<std::vector<int>>& member_pixels,
    const std::vector<double>& required_exposure_seconds,
    double minimum_exposure_seconds,
    double maximum_exposure_seconds,
    double physical_budget_seconds,
    const std::string& algorithm,
    int max_depth,
    long long reward_levels,
    double time_quantum,
    double time_limit_seconds,
    int thread_count,
    double mip_gap) {

    const OrienteeringGraph graph = buildP2Graph(
        slew_seconds,
        pixel_probabilities,
        member_pixels,
        required_exposure_seconds,
        minimum_exposure_seconds,
        maximum_exposure_seconds,
        physical_budget_seconds);

    return solveGraph(
        graph,
        slew_seconds,
        pixel_probabilities,
        true,
        algorithm,
        max_depth,
        reward_levels,
        time_quantum,
        time_limit_seconds,
        thread_count,
        mip_gap);
}