#pragma once
#include <vector>
#include <string>
#include <list>
#include <unordered_map>

#include <lemon/list_graph.h>
#include <lemon/kruskal.h>
#include <lemon/matching.h>

// std::vector<int> mstNaive(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);             

// std::vector<int> mstNaiveUpdate(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);                  

// std::vector<int> mstLemon3(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);

// std::vector<int> mstLemon2(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);

// std::vector<int> mstLemon(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);


// std::vector<int> mstLemonUpdate(const std::vector<std::vector<double>>& costs, 
//                           std::vector<double>& prizes, 
//                           double budget, 
//                           int init_pos);


void sortNodes(const std::vector<double>& prizes, std::vector<double>& prizes_sorted, 
               std::vector<int>& indices, int init_pos);
void sortNodesInRange(const std::vector<double>& prizes, std::vector<int>& indices, std::pair<int, int> range);
void selectNodesPrizeRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<double>& prizes_selected, 
                            std::vector<std::vector<double>>& costs_selected, 
                            std::vector<int>& indices, 
                            int init_pos, double cost_threshold);
// void chooseFrontNodes(const std::vector<int>& indices,
//                     const std::vector<std::vector<double>>& costs,
//                     std::vector<std::vector<double>>& costs_filtered, int topN);
// void chooseFrontNodes(const std::vector<std::vector<double>>& costs,
//                     std::vector<std::vector<double>>& costs_filtered, int topN);
// std::list<std::pair<int, int>> buildMSTKruskalUpdate(const std::vector<std::vector<double>>& costs, 
//                                                 const std::vector<int>& indices, int N_nodes);       
// double mstCosts(const std::vector<std::vector<double>>& costs,const std::list<std::pair<int, int>>& mst); 
// std::unordered_map<int, std::list<int>> buildAdjacencyList(const std::list<std::pair<int, int>>& edges);
// std::vector<int> findEulerianTour(std::unordered_map<int, std::list<int>>& adj, int start);
// std::vector<int> shortcutTour(const std::vector<int>& eulerTour, int startNode);        

typedef lemon::ListGraph Graph;
typedef Graph::Node Node;
typedef Graph::Edge Edge;

void buildGraphUpdate(const std::vector<std::vector<double>>& costs, 
                      const std::vector<int>& indices, int N_nodes,
                      Graph& g, Graph::EdgeMap<double>& weight);
void buildMST(const Graph& g, const Graph::EdgeMap<double>& weight, std::vector<Edge>& mst_edges);
double computeMSTWeight(const Graph& g, const Graph::EdgeMap<double>& weight, 
                        const std::vector<Edge>& mst_edges);


std::vector<int> GCP_BestPrizeOneFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);
// std::vector<int> mstHoogeveen(const std::vector<std::vector<double>>& costs, 
//                            const std::vector<double>& prizes, 
//                            double budget, 
//                            int init_pos);

std::vector<int> GCP_PrizeBoostOneFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, 
                          int init_pos);

void buildGraph(const std::vector<std::vector<double>>& costs, Graph& g, Graph::EdgeMap<double>& weight);

std::vector<int> GCP_BestPrizeTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx);

std::vector<int> GCP_PrizeBoostTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx);

std::vector<int> GCP_PrizeRatioTwoFix(const std::vector<std::vector<double>>& costs, 
                          const std::vector<double>& prizes, 
                          double budget, int s_idx, int t_idx);

std::vector<int> GCP(const std::vector<std::vector<double>>& costs, 
                                    const std::vector<double>& prizes, 
                                    double budget, int s_idx, int t_idx);



double sumPrizes(const std::vector<double>& prizes, const std::vector<int>& path);

std::vector<int> christofidesPathTwoFixed(const Graph& g,
                         const std::vector<std::vector<double>>& costs,
                         const std::vector<int>& indices,
                         int N_nodes,
                         const std::vector<Edge>& mst_edges,
                         int s_idx,
                         int t_idx);

void greedyExtendPath(const std::vector<std::vector<double>>& costs, const std::vector<double>& prizes,
                      std::vector<int>& path, double& total_cost, double budget);


// void selectNodesPrizeRatio_2(const std::vector<std::vector<double>>& costs,
//                             const std::vector<double>& prizes,
//                             std::vector<int>& indices, 
//                             std::vector<bool>& selected,
//                             std::vector<double>& minDist,
//                             std::priority_queue<std::pair<double, int>>& heap,
//                             int num);

// void selectNodesPrizeRatio_(const std::vector<std::vector<double>>& costs,
//                             const std::vector<double>& prizes,
//                             std::vector<int>& indices, 
//                             std::vector<bool>& selected,
//                             std::vector<double>& minDist,
//                             std::priority_queue<std::pair<double, int>>& heap,
//                             int num);

std::pair<std::vector<std::vector<int>>, double> GCP2(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes, 
                    double budget, int s, int t, int num_paths);

std::pair<std::vector<std::vector<int>>, double> GCP2(const std::vector<std::vector<double>>& costs, 
                    const std::vector<double>& prizes, 
                    double budget, const std::vector<int>& s, int t, int num_paths);