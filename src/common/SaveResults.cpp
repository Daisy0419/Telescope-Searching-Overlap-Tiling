#include "SaveResults.h"
#include "ReadData.h"
#include "Genetic.h"
#include "Greedy.h"
#include "ReadData.h"
#include "ILP_gurobi_top.h"
#include "GCP.h"


#include <string>
#include <iostream>
#include <sstream>
#include <filesystem> 
#include <numeric>
#include <stdexcept> 
#include <vector>
#include <iomanip>


// //save result to stdout
// double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
//                     std::vector<int> rank, std::vector<int> best_path, double padding) {

//     if (best_path.empty()) {
//         throw std::runtime_error(
//             "No valid path found. Please check your Gurobi license is mounted correctly."
//         );
//     }

//     if (rank.empty()) {
//         rank.resize(best_path.size());
//         std::iota(rank.begin(), rank.end(), 1); 
//     }
    
//     double sum_probability = 0.0;
//     std::cout << "Path(tile rank): ";
//     for (int tile : best_path) {
//         std::cout << tile << " "; 
//         // std::cout << rank[tile] << " "; //ranks in tiling file
//         sum_probability += probability[tile]; 
//     }
//     std::cout << std::endl;

//     std::cout << "Path(in Original idx): ";
//     for (int tile : best_path) {
//         // std::cout << tile << " "; 
//         std::cout << rank[tile] << " "; //ranks in tiling file
//     }
//     std::cout << std::endl;

//     double total_cost = -padding;
//     for (size_t i = 0; i < best_path.size() - 1; i++) {
//         total_cost += costs[best_path[i]][best_path[i + 1]]; 
//     }
//     std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
//     std::cout << "Sum Probability: " << sum_probability << std::endl; 
//     std::cout << "Total Cost: " << total_cost << std::endl;

//     return sum_probability;
// }


// //save result to a csv file
// void save_result(const std::string& filename, 
//                   const std::string& method, 
//                   const std::string& data, 
//                   double budget, double slew_rate, 
//                   const std::vector<std::vector<double>>& costs, 
//                   const std::vector<double>& probability, 
//                   std::vector<int> rank,
//                   std::vector<int> best_path, 
//                   double elapsed_time, double padding) {

//     bool file_exists = std::filesystem::exists(filename); 

//     std::ofstream outfile(filename, std::ios::app);  

//     if (!outfile.is_open()) {
//         std::cerr << "Error: Could not open results file: " << filename << std::endl;
//         return;
//     }

//     //Write header
//     if (!file_exists) {
//         outfile << "Method,Dataset,Budget,SlewRate,NumTiles,SumProb,TotalCost,TimeSec,Path\n";
//     }

//     double sum_probability = 0.0;
//     double total_cost = -padding;

//     if (!best_path.empty()) {
//         for (size_t i = 0; i < best_path.size(); i++) {
//             sum_probability += probability[best_path[i]];
//             if (i < best_path.size() - 1) {
//                 total_cost += costs[best_path[i]][best_path[i + 1]];
//             }
//         }
//     }

//     if (rank.empty()) {
//         rank.resize(best_path.size());
//         std::iota(rank.begin(), rank.end(), 1); 
//     }

//     // Write result to CSV
//     outfile << method << ","
//             << std::filesystem::path(data).filename().string() << ","
//             << budget << ","
//             << slew_rate << ","
//             << best_path.size() << ","
//             << sum_probability << ","
//             << total_cost << ","
//             << elapsed_time << ",";

//     // Write path
//     for (size_t i = 0; i < best_path.size(); i++) {
//         outfile << rank[best_path[i]];
//         if (i < best_path.size() - 1) {
//             outfile << " ";  
//         }
//     }
//     outfile << "\n";

//     outfile.close();

// }

// // save the result of multi deadline result in csv
// void save_result2(const std::string& filename, 
//                   const std::string& method, 
//                   const std::vector<double>& budgets,
//                   const std::vector<double>& dwell_time, 
//                   const std::vector<std::vector<double>>& costs, 
//                   const std::vector<double>& probability, 
//                   std::vector<int> rank,
//                   std::vector<int> best_path, 
//                   double elapsed_time) {

//     bool file_exists = std::filesystem::exists(filename); 

//     std::ofstream outfile(filename, std::ios::app);  

//     if (!outfile.is_open()) {
//         std::cerr << "Error: Could not open results file: " << filename << std::endl;
//         return;
//     }

//     //Write header
//     if (!file_exists) {
//         outfile << "Method,NumTiles,TimeSec,Budgets,Costs,FOMs,Path\n";
//     }

//     double sum_probability = 0.0;
//     double total_cost = 0.0;
//     std::vector<double> foms;
//     std::vector<double> accu_costs;
//     int b = 0;

//     if (!best_path.empty()) {
//         for (size_t i = 1; i < best_path.size(); i++) {
//             total_cost += costs[best_path[i]][best_path[i - 1]];

//             if(b < budgets.size() && total_cost + 0.5*dwell_time[best_path[i]] > budgets[b]){
//                 foms.push_back(sum_probability);
//                 sum_probability = 0.0;
//                 double accumulated_cost=(total_cost - costs[best_path[i]][best_path[i - 1]] + 0.5*dwell_time[best_path[i-1]]);
//                 accu_costs.push_back(accumulated_cost);
//                 b++;
//             }

//             sum_probability += probability[best_path[i]];
//         }
//     }

//     if (rank.empty()) {
//         rank.resize(best_path.size());
//         std::iota(rank.begin(), rank.end(), 1); 
//     }

//     // Write result to CSV
//     outfile << method << ","
//             << best_path.size() << ","
//             << elapsed_time << ",";


//     // Write budgets
//     for (size_t i = 0; i < budgets.size(); i++) {
//         outfile << budgets[i];
//         if (i < budgets.size() - 1) {
//             outfile << " ";  
//         }
//     }
//     outfile << ",";

//     // Write Costs
//     for (size_t i = 0; i < accu_costs.size(); i++) {
//         outfile << accu_costs[i];
//         if (i < accu_costs.size() - 1) {
//             outfile << " ";  
//         }
//     }
//     outfile << ",";

//     // Write foms
//     for (size_t i = 0; i < foms.size(); i++) {
//         outfile << foms[i];
//         if (i < foms.size() - 1) {
//             outfile << " ";  
//         }
//     }
//     outfile << ",";

//     // Write path
//     for (size_t i = 0; i < best_path.size(); i++) {
//         outfile << rank[best_path[i]];
//         if (i < best_path.size() - 1) {
//             outfile << " ";  
//         }
//     }
//     outfile << "\n";

//     outfile.close();

// }


