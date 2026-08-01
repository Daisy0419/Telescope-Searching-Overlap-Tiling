#pragma once

#include <vector>
#include <string>
#include <utility>
 
/**
 * Problem data for a Team Orienteering / Min-Routes-Cover problem
 * with a single common start node and a single common end node.
 */
struct ProblemDataMinKCommonR {
    int N;                                      // number of nodes
    int K;                                      // upper bound on number of teams/routes
    int startIndex;                             // common start node s
    int endIndex;                               // common end node t
    double BudgetPerTeam;                       // per-team travel budget
    std::vector<std::vector<double>> Cost;      // NxN cost/distance matrix
    std::vector<double> Prize;                  // N-length prize vector (unused in min-routes,
                                                // but kept for interface compatibility)
 
    ProblemDataMinKCommonR(int N_, int K_, int s_, int t_, double B_,
                      const std::vector<std::vector<double>>& C_,
                      const std::vector<double>& P_)
        : N(N_), K(K_), startIndex(s_), endIndex(t_),
          BudgetPerTeam(B_), Cost(C_), Prize(P_) {}
};
 

std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverCommonR(const ProblemDataMinKCommonR& data,
                    double mipGap       = -1.0,
                    double timeLimit    = -1.0,
                    const std::vector<std::vector<int>>& initialRoutes = {});


struct ProblemDataMinKDiffR {
    int N;                                      // number of nodes
    int K;                                      // upper bound on number of teams/routes
    std::vector<int> startIndices;              // per-route start node s_k, size K
    int endIndex;                               // common end node t
    double BudgetPerTeam;                       // per-team travel budget
    std::vector<std::vector<double>> Cost;      // NxN cost/distance matrix
    std::vector<double> Prize;                  // N-length prize vector

    ProblemDataMinKDiffR(int N_, int K_,
                      const std::vector<int>& starts_, int t_, double B_,
                      const std::vector<std::vector<double>>& C_,
                      const std::vector<double>& P_)
        : N(N_), K(K_), startIndices(starts_), endIndex(t_),
          BudgetPerTeam(B_), Cost(C_), Prize(P_) {}
};


std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverDiffR(const ProblemDataMinKDiffR& data,
                    double mipGap       = -1.0,
                    double timeLimit    = -1.0,
                    const std::vector<std::vector<int>>& initialRoutes = {});



struct ProblemDataMinKUnrooted {
    int N;                                       // number of nodes
    int K;                                       // upper bound on number of routes
    double BudgetPerTeam;                        // per-route budget
    std::vector<std::vector<double>> SlewCost;   // NxN edge (slew) cost matrix
    std::vector<double> DwellCost;               // N-length node (dwell) cost vector

    ProblemDataMinKUnrooted(int N_, int K_, double B_,
                            const std::vector<std::vector<double>>& SC_,
                            const std::vector<double>& DC_)
        : N(N_), K(K_), BudgetPerTeam(B_), SlewCost(SC_), DwellCost(DC_) {}
};


std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverUnrooted(const ProblemDataMinKUnrooted& data,
                            double mipGap       = -1.0,
                            double timeLimit    = -1.0,
                            const std::vector<std::vector<int>>& initialRoutes = {});



struct ProblemDataMinKRootedOpen {
    int N;                                       // number of customer nodes
    int K;                                       // upper bound on number of routes
    std::vector<int> startIndices;               // per-route start node, size K
    double BudgetPerTeam;                        // per-route budget
    std::vector<std::vector<double>> SlewCost;   // NxN edge (slew) cost matrix
    std::vector<double> DwellCost;               // N-length node (dwell) cost vector

    ProblemDataMinKRootedOpen(int N_, int K_,
                              const std::vector<int>& starts_, double B_,
                              const std::vector<std::vector<double>>& SC_,
                              const std::vector<double>& DC_)
        : N(N_), K(K_), startIndices(starts_), BudgetPerTeam(B_),
          SlewCost(SC_), DwellCost(DC_) {}
};


std::pair<std::vector<std::vector<int>>, double>
solveMinRoutesCoverRootedOpen(const ProblemDataMinKRootedOpen& data,
                              double mipGap       = -1.0,
                              double timeLimit    = -1.0,
                              const std::vector<std::vector<int>>& initialRoutes = {});

std::pair<std::vector<std::vector<int>>, double>
MinKCoverUnrooted_ILP(double budget,
                      const std::vector<std::vector<double>>& slew_costs,
                      const std::vector<double>& dwell_costs,
                      double mipGap       = -1.0,
                      double timeLimit    = -1.0,
                      const std::vector<std::vector<int>>& initialRoutes = {});

std::pair<std::vector<std::vector<int>>, double>
MinKCoverRootedOpen_ILP(const std::vector<int>& starts, double budget,
                        const std::vector<std::vector<double>>& slew_costs,
                        const std::vector<double>& dwell_costs,
                            double mipGap       = -1.0,
                            double timeLimit    = -1.0,
                            const std::vector<std::vector<int>>& initialRoutes = {});
