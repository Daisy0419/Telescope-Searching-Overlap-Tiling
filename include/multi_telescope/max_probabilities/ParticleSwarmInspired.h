#pragma once
#include <functional>

// Same alias you already use
using SplitterFn = std::function<
    std::pair<std::vector<std::vector<int>>, double>(
        const Instance&, double, const std::vector<int>&, int)>;

std::pair<std::vector<std::vector<int>>, double>
run_psoia_team_op(const Instance& I, double B, int m_routes,
                  int swarm_size, int k_stop, unsigned seed,
                  const SplitterFn& splitter,
                  double ph=0.10, double d=1e-2,
                  double w0=0.9, double c1=0.5, double c2=0.5);