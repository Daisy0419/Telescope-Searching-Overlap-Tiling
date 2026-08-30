#include "SubmodularSelection.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>


CoverageOracle::CoverageOracle(const std::vector<double>& pixel_probs,
                               const std::vector<std::vector<int>>& member_pixels)
                            : pixel_probs_(pixel_probs),
                            member_(member_pixels),
                            cover_count_(pixel_probs.size(), 0),
                            stamp_(pixel_probs.size(), 0) {

    for (double p : pixel_probs_)
        if (!std::isfinite(p) || p < 0.0)
            throw std::invalid_argument(
                "CoverageOracle: pixel probabilities must be finite and nonnegative.");

    const int pixel_count = static_cast<int>(pixel_probs_.size());
    for (std::size_t node = 0; node < member_.size(); ++node) {
        auto& pixels = member_[node];
        std::sort(pixels.begin(), pixels.end());
        pixels.erase(std::unique(pixels.begin(), pixels.end()), pixels.end());
        for (int p : pixels)
            if (p < 0 || p >= pixel_count)
                throw std::out_of_range(
                    "CoverageOracle: node " + std::to_string(node) +
                    " has out-of-range pixel " + std::to_string(p) + ".");
    }
}


void CoverageOracle::reset() {
    std::fill(cover_count_.begin(), cover_count_.end(), 0);
    committed_ = 0.0;
}


void CoverageOracle::commit(int node) {
    for (int p : member_[node])
        if (cover_count_[p]++ == 0) committed_ += pixel_probs_[p];
}

void CoverageOracle::commitAll(const std::vector<int>& nodes) {
    for (int v : nodes) commit(v);
}

//not change the committed state
double CoverageOracle::marginal(int node) const {
    double gain = 0.0;
    for (int p : member_[node])
        if (cover_count_[p] == 0) gain += pixel_probs_[p];
    return gain;
}


double CoverageOracle::rewardOf(const std::vector<int>& nodes) const {
    ++epoch_;
    double total = 0.0;
    for (int v : nodes)
        for (int p : member_[v])
            if (stamp_[p] != epoch_) {
                stamp_[p] = epoch_;
                total += pixel_probs_[p];
            }
    return total;
}


double CoverageOracle::marginalOver(int node, const std::vector<int>& nodes) const {
    ++epoch_;
    for (int v : nodes)
        for (int p : member_[v]) stamp_[p] = epoch_;
    double gain = 0.0;
    for (int p : member_[node])
        if (stamp_[p] != epoch_) {
            stamp_[p] = epoch_;
            gain += pixel_probs_[p];
        }
    return gain;
}

