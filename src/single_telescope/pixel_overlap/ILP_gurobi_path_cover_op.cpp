#include "ILP_gurobi_path_cover_op.h"

#include "gurobi_c++.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>


// --------------------------------------------------------------------------
// MTZ formulation: subtours are eliminated only by ordering variables.
// --------------------------------------------------------------------------
std::pair<std::vector<std::vector<int>>, double>
gurobiPathCoverPixelOverlapMTZ(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap,
    double time_limit,
    const std::vector<int>& initial_route) {

    double best_bound = std::numeric_limits<double>::quiet_NaN();
    const int node_count = static_cast<int>(costs.size());
    const int pixel_count = static_cast<int>(pixel_probs.size());

    // 1. Validate the all-in s-t graph and compact pixel coverage.
    if (node_count < 2) {
        throw std::invalid_argument(
            "MTZ path cover: graph must contain at least two nodes.");
    }
    if (start < 0 || start >= node_count ||
        end < 0 || end >= node_count || start == end) {
        throw std::invalid_argument(
            "MTZ path cover: invalid or identical start/end nodes.");
    }
    if (!std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "MTZ path cover: budget must be finite and nonnegative.");
    }
    if (!std::isfinite(mip_gap) || !std::isfinite(time_limit)) {
        throw std::invalid_argument(
            "MTZ path cover: MIP gap and time limit must be finite.");
    }
    if (static_cast<int>(member_pixels.size()) != node_count) {
        throw std::invalid_argument(
            "MTZ path cover: member_pixels must have one list per node.");
    }

    for (int i = 0; i < node_count; ++i) {
        if (static_cast<int>(costs[i].size()) != node_count) {
            throw std::invalid_argument(
                "MTZ path cover: cost matrix must be square.");
        }
        for (int j = 0; j < node_count; ++j) {
            if (!std::isfinite(costs[i][j]) || costs[i][j] < 0.0) {
                throw std::invalid_argument(
                    "MTZ path cover: costs must be finite and nonnegative.");
            }
        }
        for (int pixel : member_pixels[i]) {
            if (pixel < 0 || pixel >= pixel_count) {
                throw std::out_of_range(
                    "MTZ path cover: invalid compact HEALPixel id.");
            }
        }
    }
    for (double probability : pixel_probs) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "MTZ path cover: pixel probabilities must be finite and "
                "nonnegative.");
        }
    }
    if (!member_pixels[start].empty() || !member_pixels[end].empty()) {
        throw std::invalid_argument(
            "MTZ path cover: dummy start/end pixel lists must be empty.");
    }

    const auto usable_arc = [&](int from, int to) {
        return from != to && from != end && to != start &&
               costs[from][to] <= budget + 1e-12;
    };

    // Validate the optional warm-start route.
    bool use_warm_start = !initial_route.empty();
    std::string warm_start_error;
    if (use_warm_start) {
        if (initial_route.front() != start || initial_route.back() != end) {
            use_warm_start = false;
            warm_start_error =
                "route must begin at start and end at the dummy end";
        } else {
            std::vector<char> seen(node_count, 0);
            double route_cost = 0.0;
            for (std::size_t k = 0;
                 k < initial_route.size() && use_warm_start; ++k) {
                const int node = initial_route[k];
                if (node < 0 || node >= node_count || seen[node]) {
                    use_warm_start = false;
                    warm_start_error =
                        "route contains an invalid or repeated node";
                    break;
                }
                seen[node] = 1;
                if (k + 1 < initial_route.size()) {
                    const int next = initial_route[k + 1];
                    if (next < 0 || next >= node_count ||
                        !usable_arc(node, next)) {
                        use_warm_start = false;
                        warm_start_error =
                            "route contains an invalid or over-budget arc";
                        break;
                    }
                    route_cost += costs[node][next];
                }
            }
            if (use_warm_start && route_cost > budget + 1e-9) {
                use_warm_start = false;
                warm_start_error =
                    "route exceeds the already padded budget";
            }
        }
    }

    try {
        // 2. Create variables.
        GRBEnv environment(true);
        environment.start();
        GRBModel model(environment);
        model.set(GRB_StringAttr_ModelName, "ST_Pixel_Path_Cover_MTZ");

        std::vector<GRBVar> visit(node_count);
        for (int node = 0; node < node_count; ++node) {
            visit[node] = model.addVar(
                0.0, 1.0, 0.0, GRB_CONTINUOUS,
                "visit_" + std::to_string(node));
        }

        std::vector<GRBVar> pixel_covered(pixel_count);
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            pixel_covered[pixel] = model.addVar(
                0.0, 1.0, 0.0, GRB_CONTINUOUS,
                "pixel_" + std::to_string(pixel));
        }

        std::vector<std::vector<GRBVar>> arc(
            node_count, std::vector<GRBVar>(node_count));
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (i == j) {
                    continue;
                }
                arc[i][j] = model.addVar(
                    0.0, usable_arc(i, j) ? 1.0 : 0.0,
                    0.0, GRB_BINARY,
                    "arc_" + std::to_string(i) + "_" +
                        std::to_string(j));
            }
        }

        std::vector<GRBVar> order(node_count);
        for (int node = 0; node < node_count; ++node) {
            order[node] = model.addVar(
                0.0, static_cast<double>(node_count - 1),
                0.0, GRB_CONTINUOUS,
                "order_" + std::to_string(node));
        }
        model.update();

        // 3. Maximize the union probability of covered HEALPixels.
        GRBLinExpr objective = 0.0;
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            objective += pixel_probs[pixel] * pixel_covered[pixel];
        }
        model.setObjective(objective, GRB_MAXIMIZE);

        std::vector<GRBLinExpr> cover_expression(pixel_count, 0.0);
        for (int node = 0; node < node_count; ++node) {
            std::unordered_set<int> pixels_seen_in_node;
            for (int pixel : member_pixels[node]) {
                if (pixels_seen_in_node.insert(pixel).second) {
                    cover_expression[pixel] += visit[node];
                }
            }
        }
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            model.addConstr(
                pixel_covered[pixel] <= cover_expression[pixel],
                "pixel_covered_" + std::to_string(pixel));
        }

        // 4. Enforce one dummy-start to dummy-end path.
        model.addConstr(visit[start] == 1.0, "start_active");
        model.addConstr(visit[end] == 1.0, "end_active");

        GRBLinExpr start_in = 0.0;
        GRBLinExpr start_out = 0.0;
        GRBLinExpr end_in = 0.0;
        GRBLinExpr end_out = 0.0;
        for (int node = 0; node < node_count; ++node) {
            if (node != start) {
                start_in += arc[node][start];
                start_out += arc[start][node];
            }
            if (node != end) {
                end_in += arc[node][end];
                end_out += arc[end][node];
            }
        }
        model.addConstr(start_in == 0.0, "start_in_zero");
        model.addConstr(start_out == 1.0, "start_out_one");
        model.addConstr(end_in == 1.0, "end_in_one");
        model.addConstr(end_out == 0.0, "end_out_zero");

        for (int node = 0; node < node_count; ++node) {
            if (node == start || node == end) {
                continue;
            }
            GRBLinExpr incoming = 0.0;
            GRBLinExpr outgoing = 0.0;
            for (int other = 0; other < node_count; ++other) {
                if (other != node) {
                    incoming += arc[other][node];
                    outgoing += arc[node][other];
                }
            }
            model.addConstr(
                incoming == visit[node],
                "node_in_" + std::to_string(node));
            model.addConstr(
                outgoing == visit[node],
                "node_out_" + std::to_string(node));
        }

        // 5. MTZ ordering constraints eliminate all subtours.
        model.addConstr(order[start] == 0.0, "order_start");
        for (int node = 0; node < node_count; ++node) {
            if (node == start) {
                continue;
            }
            model.addConstr(
                order[node] <= (node_count - 1) * visit[node],
                "order_upper_" + std::to_string(node));
            model.addConstr(
                order[node] >= visit[node],
                "order_lower_" + std::to_string(node));
        }

        const double order_big_m = static_cast<double>(node_count);
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (!usable_arc(i, j)) {
                    continue;
                }
                model.addConstr(
                    order[j] >= order[i] + 1.0 -
                        order_big_m * (1.0 - arc[i][j]),
                    "order_arc_" + std::to_string(i) + "_" +
                        std::to_string(j));
            }
        }

        // 6. The budget is already padded and costs are already all-in.
        GRBLinExpr total_cost = 0.0;
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (usable_arc(i, j)) {
                    total_cost += costs[i][j] * arc[i][j];
                }
            }
        }
        model.addConstr(total_cost <= budget, "budget");

        if (mip_gap >= 0.0) {
            model.set(GRB_DoubleParam_MIPGap, mip_gap);
        }
        if (time_limit > 0.0) {
            model.set(GRB_DoubleParam_TimeLimit, time_limit);
        }

        // 7. Apply the optional warm start.
        if (use_warm_start) {
            std::vector<char> in_route(node_count, 0);
            for (int node : initial_route) {
                in_route[node] = 1;
            }
            for (int i = 0; i < node_count; ++i) {
                visit[i].set(
                    GRB_DoubleAttr_Start, in_route[i] ? 1.0 : 0.0);
                order[i].set(GRB_DoubleAttr_Start, 0.0);
                for (int j = 0; j < node_count; ++j) {
                    if (i != j) {
                        arc[i][j].set(GRB_DoubleAttr_Start, 0.0);
                    }
                }
            }
            for (std::size_t k = 0; k < initial_route.size(); ++k) {
                const int node = initial_route[k];
                order[node].set(
                    GRB_DoubleAttr_Start, static_cast<double>(k));
                if (k + 1 < initial_route.size()) {
                    arc[node][initial_route[k + 1]].set(
                        GRB_DoubleAttr_Start, 1.0);
                }
            }
        } else if (!initial_route.empty()) {
            std::cerr << "MTZ warm start skipped: "
                      << warm_start_error << ".\n";
        }

        // 8. Optimize, retain ObjBound, and extract the incumbent path.
        model.optimize();
        const int status = model.get(GRB_IntAttr_Status);
        const int solution_count = model.get(GRB_IntAttr_SolCount);
        try {
            best_bound = model.get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException&) {
            // Leave the bound as NaN when unavailable.
        }

        if (solution_count <= 0) {
            std::cerr << "MTZ path-cover model: no incumbent, best_bound="
                      << best_bound << ", status=" << status << ".\n";
            return {std::vector<std::vector<int>>{}, best_bound};
        }

        const double incumbent = model.get(GRB_DoubleAttr_ObjVal);
        const double gap = model.get(GRB_DoubleAttr_MIPGap);
        std::cout << "MTZ path-cover model: incumbent=" << incumbent
                  << ", best_bound=" << best_bound
                  << ", gap=" << gap
                  << ", status=" << status << "\n";

        std::vector<int> path;
        std::vector<char> seen(node_count, 0);
        int current = start;
        while (current >= 0 && current < node_count && !seen[current]) {
            path.push_back(current);
            seen[current] = 1;
            if (current == end) {
                return {{std::move(path)}, best_bound};
            }

            int successor = -1;
            double largest_value = 0.5;
            for (int next = 0; next < node_count; ++next) {
                if (next == current) {
                    continue;
                }
                const double value =
                    arc[current][next].get(GRB_DoubleAttr_X);
                if (value > largest_value) {
                    largest_value = value;
                    successor = next;
                }
            }
            if (successor < 0) {
                break;
            }
            current = successor;
        }

        std::cerr << "MTZ path-cover model: failed to extract a complete "
                     "s-t path.\n";
        return {std::vector<std::vector<int>>{}, best_bound};

    } catch (const GRBException& error) {
        std::cerr << "MTZ path-cover model: Gurobi error "
                  << error.getErrorCode() << ": "
                  << error.getMessage() << "\n";
        return {std::vector<std::vector<int>>{}, best_bound};
    }
}


// --------------------------------------------------------------------------
// Continuous-time formulation: positive cumulative edge time alone eliminates
// subtours. No MTZ/order variables or MTZ constraints are present.
// --------------------------------------------------------------------------
std::pair<std::vector<std::vector<int>>, double>
gurobiPathCoverPixelOverlapContinuousTime(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start,
    int end,
    double budget,
    double mip_gap,
    double time_limit,
    const std::vector<int>& initial_route) {

    double best_bound = std::numeric_limits<double>::quiet_NaN();
    const int node_count = static_cast<int>(costs.size());
    const int pixel_count = static_cast<int>(pixel_probs.size());

    // 1. Validate the all-in s-t graph and compact pixel coverage.
    if (node_count < 2) {
        throw std::invalid_argument(
            "Continuous-time path cover: graph must contain at least two "
            "nodes.");
    }
    if (start < 0 || start >= node_count ||
        end < 0 || end >= node_count || start == end) {
        throw std::invalid_argument(
            "Continuous-time path cover: invalid or identical start/end "
            "nodes.");
    }
    if (!std::isfinite(budget) || budget < 0.0) {
        throw std::invalid_argument(
            "Continuous-time path cover: budget must be finite and "
            "nonnegative.");
    }
    if (!std::isfinite(mip_gap) || !std::isfinite(time_limit)) {
        throw std::invalid_argument(
            "Continuous-time path cover: MIP gap and time limit must be "
            "finite.");
    }
    if (static_cast<int>(member_pixels.size()) != node_count) {
        throw std::invalid_argument(
            "Continuous-time path cover: member_pixels must have one list "
            "per node.");
    }

    for (int i = 0; i < node_count; ++i) {
        if (static_cast<int>(costs[i].size()) != node_count) {
            throw std::invalid_argument(
                "Continuous-time path cover: cost matrix must be square.");
        }
        for (int j = 0; j < node_count; ++j) {
            if (!std::isfinite(costs[i][j]) || costs[i][j] < 0.0) {
                throw std::invalid_argument(
                    "Continuous-time path cover: costs must be finite and "
                    "nonnegative.");
            }
        }
        for (int pixel : member_pixels[i]) {
            if (pixel < 0 || pixel >= pixel_count) {
                throw std::out_of_range(
                    "Continuous-time path cover: invalid compact HEALPixel "
                    "id.");
            }
        }
    }
    for (double probability : pixel_probs) {
        if (!std::isfinite(probability) || probability < 0.0) {
            throw std::invalid_argument(
                "Continuous-time path cover: pixel probabilities must be "
                "finite and nonnegative.");
        }
    }
    if (!member_pixels[start].empty() || !member_pixels[end].empty()) {
        throw std::invalid_argument(
            "Continuous-time path cover: dummy start/end pixel lists must "
            "be empty.");
    }

    const auto usable_arc = [&](int from, int to) {
        return from != to && from != end && to != start &&
               costs[from][to] <= budget + 1e-12;
    };

    // A disconnected flow component must contain a directed cycle. Summing
    // the selected time-propagation constraints around any positive-time cycle
    // gives the contradiction t >= t + positive_value. The only unsafe case
    // is therefore an all-zero-time directed cycle among non-dummy nodes.
    constexpr double zero_time_tolerance = 1e-9;
    std::vector<std::vector<int>> zero_time_adjacency(node_count);
    for (int i = 0; i < node_count; ++i) {
        if (i == start || i == end) {
            continue;
        }
        for (int j = 0; j < node_count; ++j) {
            if (j == start || j == end) {
                continue;
            }
            if (usable_arc(i, j) &&
                costs[i][j] <= zero_time_tolerance) {
                zero_time_adjacency[i].push_back(j);
            }
        }
    }

    std::vector<char> zero_time_state(node_count, 0);
    std::function<bool(int)> has_zero_time_cycle = [&](int node) {
        zero_time_state[node] = 1;
        for (int next : zero_time_adjacency[node]) {
            if (zero_time_state[next] == 1) {
                return true;
            }
            if (zero_time_state[next] == 0 && has_zero_time_cycle(next)) {
                return true;
            }
        }
        zero_time_state[node] = 2;
        return false;
    };

    for (int node = 0; node < node_count; ++node) {
        if (node != start && node != end &&
            zero_time_state[node] == 0 && has_zero_time_cycle(node)) {
            throw std::invalid_argument(
                "Continuous-time path cover: the graph contains a zero-time "
                "directed cycle among real tiles. Actual-time propagation "
                "cannot eliminate that subtour; use MTZ or ensure positive "
                "tile dwell/edge times.");
        }
    }

    // Validate the optional warm-start route.
    bool use_warm_start = !initial_route.empty();
    std::string warm_start_error;
    if (use_warm_start) {
        if (initial_route.front() != start || initial_route.back() != end) {
            use_warm_start = false;
            warm_start_error =
                "route must begin at start and end at the dummy end";
        } else {
            std::vector<char> seen(node_count, 0);
            double route_cost = 0.0;
            for (std::size_t k = 0;
                 k < initial_route.size() && use_warm_start; ++k) {
                const int node = initial_route[k];
                if (node < 0 || node >= node_count || seen[node]) {
                    use_warm_start = false;
                    warm_start_error =
                        "route contains an invalid or repeated node";
                    break;
                }
                seen[node] = 1;
                if (k + 1 < initial_route.size()) {
                    const int next = initial_route[k + 1];
                    if (next < 0 || next >= node_count ||
                        !usable_arc(node, next)) {
                        use_warm_start = false;
                        warm_start_error =
                            "route contains an invalid or over-budget arc";
                        break;
                    }
                    route_cost += costs[node][next];
                }
            }
            if (use_warm_start && route_cost > budget + 1e-9) {
                use_warm_start = false;
                warm_start_error =
                    "route exceeds the already padded budget";
            }
        }
    }

    try {
        // 2. Create routing, coverage, and actual-time variables.
        GRBEnv environment(true);
        environment.start();
        GRBModel model(environment);
        model.set(
            GRB_StringAttr_ModelName,
            "ST_Pixel_Path_Cover_Continuous_Time");

        std::vector<GRBVar> visit(node_count);
        for (int node = 0; node < node_count; ++node) {
            visit[node] = model.addVar(
                0.0, 1.0, 0.0, GRB_CONTINUOUS,
                "visit_" + std::to_string(node));
        }

        std::vector<GRBVar> pixel_covered(pixel_count);
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            pixel_covered[pixel] = model.addVar(
                0.0, 1.0, 0.0, GRB_CONTINUOUS,
                "pixel_" + std::to_string(pixel));
        }

        std::vector<std::vector<GRBVar>> arc(
            node_count, std::vector<GRBVar>(node_count));
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (i == j) {
                    continue;
                }
                arc[i][j] = model.addVar(
                    0.0, usable_arc(i, j) ? 1.0 : 0.0,
                    0.0, GRB_BINARY,
                    "arc_" + std::to_string(i) + "_" +
                        std::to_string(j));
            }
        }

        std::vector<GRBVar> completion(node_count);
        for (int node = 0; node < node_count; ++node) {
            completion[node] = model.addVar(
                0.0, budget, 0.0, GRB_CONTINUOUS,
                "completion_" + std::to_string(node));
        }
        model.update();

        // 3. Maximize the union probability of covered HEALPixels.
        GRBLinExpr objective = 0.0;
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            objective += pixel_probs[pixel] * pixel_covered[pixel];
        }
        model.setObjective(objective, GRB_MAXIMIZE);

        std::vector<GRBLinExpr> cover_expression(pixel_count, 0.0);
        for (int node = 0; node < node_count; ++node) {
            std::unordered_set<int> pixels_seen_in_node;
            for (int pixel : member_pixels[node]) {
                if (pixels_seen_in_node.insert(pixel).second) {
                    cover_expression[pixel] += visit[node];
                }
            }
        }
        for (int pixel = 0; pixel < pixel_count; ++pixel) {
            model.addConstr(
                pixel_covered[pixel] <= cover_expression[pixel],
                "pixel_covered_" + std::to_string(pixel));
        }

        // 4. Enforce one dummy-start to dummy-end flow path.
        model.addConstr(visit[start] == 1.0, "start_active");
        model.addConstr(visit[end] == 1.0, "end_active");

        GRBLinExpr start_in = 0.0;
        GRBLinExpr start_out = 0.0;
        GRBLinExpr end_in = 0.0;
        GRBLinExpr end_out = 0.0;
        for (int node = 0; node < node_count; ++node) {
            if (node != start) {
                start_in += arc[node][start];
                start_out += arc[start][node];
            }
            if (node != end) {
                end_in += arc[node][end];
                end_out += arc[end][node];
            }
        }
        model.addConstr(start_in == 0.0, "start_in_zero");
        model.addConstr(start_out == 1.0, "start_out_one");
        model.addConstr(end_in == 1.0, "end_in_one");
        model.addConstr(end_out == 0.0, "end_out_zero");

        for (int node = 0; node < node_count; ++node) {
            if (node == start || node == end) {
                continue;
            }
            GRBLinExpr incoming = 0.0;
            GRBLinExpr outgoing = 0.0;
            for (int other = 0; other < node_count; ++other) {
                if (other != node) {
                    incoming += arc[other][node];
                    outgoing += arc[node][other];
                }
            }
            model.addConstr(
                incoming == visit[node],
                "node_in_" + std::to_string(node));
            model.addConstr(
                outgoing == visit[node],
                "node_out_" + std::to_string(node));
        }

        // 5. Actual cumulative time eliminates every positive-time subtour.
        model.addConstr(completion[start] == 0.0, "completion_start");
        for (int node = 0; node < node_count; ++node) {
            model.addConstr(
                completion[node] <= budget * visit[node],
                "completion_active_" + std::to_string(node));
        }
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (!usable_arc(i, j)) {
                    continue;
                }
                const double edge_time = costs[i][j];
                const double time_big_m = budget + edge_time;
                model.addConstr(
                    completion[j] >= completion[i] + edge_time -
                        time_big_m * (1.0 - arc[i][j]),
                    "time_arc_" + std::to_string(i) + "_" +
                        std::to_string(j));
            }
        }

        // Keep the explicit budget inequality as a useful strengthening of
        // the pure time-flow formulation.
        GRBLinExpr total_cost = 0.0;
        for (int i = 0; i < node_count; ++i) {
            for (int j = 0; j < node_count; ++j) {
                if (usable_arc(i, j)) {
                    total_cost += costs[i][j] * arc[i][j];
                }
            }
        }
        model.addConstr(total_cost <= budget, "budget");

        if (mip_gap >= 0.0) {
            model.set(GRB_DoubleParam_MIPGap, mip_gap);
        }
        if (time_limit > 0.0) {
            model.set(GRB_DoubleParam_TimeLimit, time_limit);
        }

        // 6. Apply the optional warm start using true accumulated edge time.
        if (use_warm_start) {
            std::vector<char> in_route(node_count, 0);
            for (int node : initial_route) {
                in_route[node] = 1;
            }
            for (int i = 0; i < node_count; ++i) {
                visit[i].set(
                    GRB_DoubleAttr_Start, in_route[i] ? 1.0 : 0.0);
                completion[i].set(GRB_DoubleAttr_Start, 0.0);
                for (int j = 0; j < node_count; ++j) {
                    if (i != j) {
                        arc[i][j].set(GRB_DoubleAttr_Start, 0.0);
                    }
                }
            }

            double elapsed = 0.0;
            for (std::size_t k = 0; k < initial_route.size(); ++k) {
                const int node = initial_route[k];
                completion[node].set(GRB_DoubleAttr_Start, elapsed);
                if (k + 1 < initial_route.size()) {
                    const int next = initial_route[k + 1];
                    arc[node][next].set(GRB_DoubleAttr_Start, 1.0);
                    elapsed += costs[node][next];
                }
            }
        } else if (!initial_route.empty()) {
            std::cerr << "Continuous-time warm start skipped: "
                      << warm_start_error << ".\n";
        }

        // 7. Optimize, retain ObjBound, and extract the incumbent path.
        model.optimize();
        const int status = model.get(GRB_IntAttr_Status);
        const int solution_count = model.get(GRB_IntAttr_SolCount);
        try {
            best_bound = model.get(GRB_DoubleAttr_ObjBound);
        } catch (const GRBException&) {
            // Leave the bound as NaN when unavailable.
        }

        if (solution_count <= 0) {
            std::cerr << "Continuous-time path-cover model: no incumbent, "
                         "best_bound="
                      << best_bound << ", status=" << status << ".\n";
            return {std::vector<std::vector<int>>{}, best_bound};
        }

        const double incumbent = model.get(GRB_DoubleAttr_ObjVal);
        const double gap = model.get(GRB_DoubleAttr_MIPGap);
        std::cout << "Continuous-time path-cover model: incumbent="
                  << incumbent
                  << ", best_bound=" << best_bound
                  << ", gap=" << gap
                  << ", status=" << status << "\n";

        std::vector<int> path;
        std::vector<char> seen(node_count, 0);
        int current = start;
        while (current >= 0 && current < node_count && !seen[current]) {
            path.push_back(current);
            seen[current] = 1;
            if (current == end) {
                return {{std::move(path)}, best_bound};
            }

            int successor = -1;
            double largest_value = 0.5;
            for (int next = 0; next < node_count; ++next) {
                if (next == current) {
                    continue;
                }
                const double value =
                    arc[current][next].get(GRB_DoubleAttr_X);
                if (value > largest_value) {
                    largest_value = value;
                    successor = next;
                }
            }
            if (successor < 0) {
                break;
            }
            current = successor;
        }

        std::cerr << "Continuous-time path-cover model: failed to extract "
                     "a complete s-t path.\n";
        return {std::vector<std::vector<int>>{}, best_bound};

    } catch (const GRBException& error) {
        std::cerr << "Continuous-time path-cover model: Gurobi error "
                  << error.getErrorCode() << ": "
                  << error.getMessage() << "\n";
        return {std::vector<std::vector<int>>{}, best_bound};
    }
}