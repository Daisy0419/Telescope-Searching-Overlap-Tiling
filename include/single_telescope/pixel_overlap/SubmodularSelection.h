#pragma once

#include <utility>
#include <vector>

// Oracle for Margional coverage
// f(S) = sum of pixel_probs over the UNION of member_pixels[v], v in S

class CoverageOracle {
public:
    CoverageOracle(const std::vector<double>& pixel_probs,
                   const std::vector<std::vector<int>>& member_pixels);

    int nodeCount() const { return static_cast<int>(member_.size()); }

    //incremental interfaces for greedy node/tile selection
    void reset();
    void commit(int node);
    void commitAll(const std::vector<int>& nodes);
    double marginal(int node) const;      // gain of node over the committed set
    double committed() const { return committed_; }

    //scoring a finished path
    double rewardOf(const std::vector<int>& nodes) const;

    //gain of `node` over `nodes`, ignoring the committed state.
    double marginalOver(int node, const std::vector<int>& nodes) const;

private:
    std::vector<double> pixel_probs_;
    std::vector<std::vector<int>> member_;
    std::vector<int> cover_count_;//how many times a pixel is covered
    double committed_ = 0.0;

    // Epoch-stamped scratch so rewardOf() need not clear a pixel-sized vector.
    mutable std::vector<unsigned> stamp_;
    mutable unsigned epoch_ = 0;
};



