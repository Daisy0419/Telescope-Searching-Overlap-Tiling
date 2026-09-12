#include "m4opt_solvers.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Convert the small C++ result into an ordinary Python dictionary.
static pybind11::dict toDictionary(const M4OptResult& result) {
    pybind11::dict output;
    output["algorithm"] = result.algorithm;
    output["tile_indices"] = result.tile_indices;
    output["exposure_seconds"] = result.exposure_seconds;
    output["start_seconds"] = result.start_seconds;
    output["coverage"] = result.coverage;
    output["duration_seconds"] = result.duration_seconds;
    output["runtime_seconds"] = result.runtime_seconds;
    output["best_bound"] = result.best_bound;
    output["completed_depth"] = result.completed_depth;
    output["has_path"] = result.has_path;
    output["timed_out"] = result.timed_out;
    return output;
}

PYBIND11_MODULE(m4opt_solvers, module) {
    module.doc() = "Simple no-window M4OPT interface for all C++ solvers";

    module.def("available_algorithms", &availableAlgorithms);

    module.def(
        "solve_p1",
        [](const std::vector<std::vector<double>>& slew_seconds,
           const std::vector<double>& pixel_probabilities,
           const std::vector<std::vector<int>>& member_pixels,
           const std::vector<double>& dwell_seconds,
           double budget_seconds,
           const std::string& algorithm,
           int max_depth,
           long long reward_levels,
           double time_quantum,
           double time_limit_seconds,
           int thread_count,
           double mip_gap) {
            M4OptResult result;
            {
                pybind11::gil_scoped_release release;
                result = solveP1(
                    slew_seconds,
                    pixel_probabilities,
                    member_pixels,
                    dwell_seconds,
                    budget_seconds,
                    algorithm,
                    max_depth,
                    reward_levels,
                    time_quantum,
                    time_limit_seconds,
                    thread_count,
                    mip_gap);
            }
            return toDictionary(result);
        },
        pybind11::arg("slew_seconds"),
        pybind11::arg("pixel_probabilities"),
        pybind11::arg("member_pixels"),
        pybind11::arg("dwell_seconds"),
        pybind11::arg("budget_seconds"),
        pybind11::arg("algorithm"),
        pybind11::arg("max_depth") = 0,
        pybind11::arg("reward_levels") = 100,
        pybind11::arg("time_quantum") = 0.01,
        pybind11::arg("time_limit_seconds") = 0.0,
        pybind11::arg("thread_count") = 0,
        pybind11::arg("mip_gap") = 0.001);

    module.def(
        "solve_p2",
        [](const std::vector<std::vector<double>>& slew_seconds,
           const std::vector<double>& pixel_probabilities,
           const std::vector<std::vector<int>>& member_pixels,
           const std::vector<double>& required_exposure_seconds,
           double minimum_exposure_seconds,
           double maximum_exposure_seconds,
           double budget_seconds,
           const std::string& algorithm,
           int max_depth,
           long long reward_levels,
           double time_quantum,
           double time_limit_seconds,
           int thread_count,
           double mip_gap) {
            M4OptResult result;
            {
                pybind11::gil_scoped_release release;
                result = solveP2(
                    slew_seconds,
                    pixel_probabilities,
                    member_pixels,
                    required_exposure_seconds,
                    minimum_exposure_seconds,
                    maximum_exposure_seconds,
                    budget_seconds,
                    algorithm,
                    max_depth,
                    reward_levels,
                    time_quantum,
                    time_limit_seconds,
                    thread_count,
                    mip_gap);
            }
            return toDictionary(result);
        },
        pybind11::arg("slew_seconds"),
        pybind11::arg("pixel_probabilities"),
        pybind11::arg("member_pixels"),
        pybind11::arg("required_exposure_seconds"),
        pybind11::arg("minimum_exposure_seconds"),
        pybind11::arg("maximum_exposure_seconds"),
        pybind11::arg("budget_seconds"),
        pybind11::arg("algorithm"),
        pybind11::arg("max_depth") = 0,
        pybind11::arg("reward_levels") = 100,
        pybind11::arg("time_quantum") = 0.01,
        pybind11::arg("time_limit_seconds") = 0.0,
        pybind11::arg("thread_count") = 0,
        pybind11::arg("mip_gap") = 0.001);
}