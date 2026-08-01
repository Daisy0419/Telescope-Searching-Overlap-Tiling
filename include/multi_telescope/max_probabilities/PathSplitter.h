#pragma once

#include "ReadData.h"


std::tuple<bool, std::vector<int>, double, double>
best_subsequence_dp(const Instance& I, double B,
                    const std::vector<int>& perm,
                    int bins = 200);

std::pair<std::vector<std::vector<int>>, double>
split_segment_dp_with_subseq(const Instance& I, double B,
                             const std::vector<int>& perm,
                             int m_routes,
                             int bins_for_subseq = 2000);

std::pair<std::vector<std::vector<int>>, double>
split_segment_dp(const Instance& I, double B,
                 const std::vector<int>& perm,
                 int m_routes);


std::pair<std::vector<std::vector<int>>, double>
split_interval_graph(const Instance& I, double B,
                     const std::vector<int>& perm,
                     int m_routes);

std::pair<std::vector<std::vector<int>>, double>
naive_split(const Instance& I, double B,
            const std::vector<int>& perm,
            int m_routes);

std::pair<std::vector<std::vector<int>>, double>
naive_split(const std::vector<std::vector<double>>& costs,
            const std::vector<double>& prizes,
            int s, int t,
            double B,
            const std::vector<int>& perm,
            int m_routes);

// std::tuple<bool, std::vector<int>, double, double>
// best_subsequence_rcsp(const Instance& I, double B, const std::vector<int>& sub);


// std::pair<std::vector<std::vector<int>>, double>
// split_segment_dp_with_subseq_exact(const Instance& I, double B,
//                                    const std::vector<int>& perm,
//                                    int m_routes);

std::pair<std::vector<std::vector<int>>, double>
greedy_split_prune_until_feasible(const std::vector<std::vector<double>>& costs,
                                  const std::vector<double>& prizes,
                                  const std::vector<int>& starts, int t,
                                  double B,
                                  const std::vector<int>& perm_in,
                                  int m_routes);