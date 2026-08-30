#pragma once
#include "SubmodularSelection.h"


//methods greedy pick. 
enum class CoverageScore {
    Marginal,   // f(v | S)
    Boosted,    // asin(f(v | S)^4) / (minDist + eps)
    Ratio       // f(v | S) / (minDist + eps)
};

// with a submodular objective the best next node depends on what is already chosen
// greedily fills positions [range.first, range.second) 
void selectNodesByCoverageInRange(CoverageOracle& oracle,
                                  std::vector<int>& indices,
                                  std::pair<int, int> range);

                                  
void selectNodesByCoverageScore(const std::vector<std::vector<double>>& costs,
                                CoverageOracle& oracle,
                                std::vector<int>& indices,
                                std::vector<bool>& selected,
                                std::vector<double>& minDist,
                                int num,
                                CoverageScore score);

void selectNodesByCoverageBoosting(const std::vector<std::vector<double>>& costs,
                                   CoverageOracle& oracle,
                                   std::vector<int>& indices,
                                   std::vector<bool>& selected,
                                   std::vector<double>& minDist,
                                   int num);

void selectNodesByCoverageRatio(const std::vector<std::vector<double>>& costs,
                                CoverageOracle& oracle,
                                std::vector<int>& indices,
                                std::vector<bool>& selected,
                                std::vector<double>& minDist,
                                int num);

// union coverage of the path, counting each pixel once.
double pathCoverage(const CoverageOracle& oracle, const std::vector<int>& path);

// Repeatedly takes the unvisited node with the largest marginal over 
// the path's current coverage and splices it in at its cheapest interior position, 
void insertHighestMarginalCheapest(const std::vector<std::vector<double>>& costs,
                                   const CoverageOracle& oracle,
                                   std::vector<int>& path,
                                   double& path_cost,
                                   double budget,
                                   std::vector<bool>& visited);

void binarySearchBestCoveragePathTwoFixed(const std::vector<std::vector<double>>& costs,
                                    const CoverageOracle& oracle,
                                    const std::vector<int>& indices,
                                    std::vector<int>& best_mst_path,
                                    int N_min, int N_max, double budget);

std::vector<int> GCP_BestCoverageTwoFix(const std::vector<std::vector<double>>& costs,
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int s_idx, int t_idx);

std::vector<int> gcpCoverageScoredTwoFix(const std::vector<std::vector<double>>& costs,
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int s_idx, int t_idx, CoverageScore score);

std::vector<int> GCP_CoverageBoostTwoFix(const std::vector<std::vector<double>>& costs,
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int s_idx, int t_idx);

std::vector<int> GCP_CoverageRatioTwoFix(const std::vector<std::vector<double>>& costs,
                                    const std::vector<double>& pixel_probs,
                                    const std::vector<std::vector<int>>& member_pixels,
                                    double budget, int s_idx, int t_idx);                                    


std::vector<int> GCP_Coverage(const std::vector<std::vector<double>>& costs,
                              const std::vector<double>& pixel_probs,
                              const std::vector<std::vector<int>>& member_pixels,
                              double budget, int s_idx, int t_idx);

                              