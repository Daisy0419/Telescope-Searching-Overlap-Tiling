#pragma once
#include <vector>
#include <stdexcept>
#include <string>
#include "gurobi_c++.h"


struct ProblemInstance {
    int nTiles;        // |V|
    int nPixels;       // |P|
    int start;         // index of start tile
    double Deadline;   // global deadline

    std::vector<double> prizes;  //reward for pixel p (0..nPixels-1)
    std::vector<double> dwell_times; //required dwell time at tile i (0..nTiles-1)
    std::vector<std::vector<double>> slew_times; //slew time from tile i to tile j
    std::vector<std::vector<int>> coverage; //list of pixels covered when observing tile i

    ProblemInstance() = default;

    ProblemInstance(int start_,
                    double deadline_,
                    std::vector<double> prizes_,
                    std::vector<double> dwell_times_,
                    std::vector<std::vector<double>> slew_times_,
                    std::vector<std::vector<int>> coverage_
                    ): start(start_),
                    Deadline(deadline_),
                    prizes(std::move(prizes_)),
                    dwell_times(std::move(dwell_times_)),
                    slew_times(std::move(slew_times_)),
                    coverage(std::move(coverage_)) {
                        
        nPixels = static_cast<int>(prizes.size());
        nTiles  = static_cast<int>(dwell_times.size());

        if (nTiles <= 0 || nPixels <= 0) {
            throw std::runtime_error("ProblemInstance: nTiles or nPixels is zero.");
        }
        if (static_cast<int>(slew_times.size()) != nTiles) {
            throw std::runtime_error("ProblemInstance: slew_times outer size != nTiles.");
        }
        for (int i = 0; i < nTiles; ++i) {
            if (static_cast<int>(slew_times[i].size()) != nTiles) {
                throw std::runtime_error("ProblemInstance: slew_times[" + std::to_string(i) +
                                         "] size != nTiles.");
            }
        }
        if (static_cast<int>(coverage.size()) != nTiles) {
            throw std::runtime_error("ProblemInstance: coverage size != nTiles.");
        }
        if (start < 0 || start >= nTiles) {
            throw std::runtime_error("ProblemInstance: start index out of range.");
        }
        if (Deadline <= 0.0) {
            throw std::runtime_error("ProblemInstance: Deadline must be positive.");
        }
    }
};


GRBModel buildContinuousTimeModel(const ProblemInstance& inst);
GRBModel buildMTZModel(const ProblemInstance& inst);

std::vector<int> gurobiPathCover(
    const std::vector<std::vector<double>>& slew_times,
    const std::vector<double>& dwell_times,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    int start, double Budget,
    double mipGap, double timeLimit,
    const std::vector<int>& initialRoute={});

    