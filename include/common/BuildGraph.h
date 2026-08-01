#pragma once

#include <vector>
#include <string>
#include <limits>
#include <optional>
#include <tuple>

std::tuple<int, int, double> buildGraphOrienteering(const std::string& filename, 
                std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, std::vector<int>& ranks, 
                std::vector<double>& dwell_times, double dwell_zenith,
                double w_max, double w_acc, double settle_time, 
                std::optional<std::pair<double, double>> zenith_pos, 
                bool is_SlowDeep = false, 
                int init_pos=-1);

std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& filename, 
                std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, std::vector<int>& ranks, 
                std::vector<double>& dwell_times, double dwell_zenith,
                double w_max, double w_acc, double settle_time, 
                bool is_SlowDeep, int m, 
                std::optional<std::pair<double, double>> zenith_pos, 
                std::vector<std::pair<double, double>> init_poses);

                
std::tuple<int, int, double> buildGraphOrienteering(const std::string& tilefile, 
                const std::string& pro_map_file,
                std::vector<std::vector<double>>& costs, 
                std::vector<double>& probability, std::vector<int>& ranks, 
                std::vector<double>& dwell_times, double dwell_zenith,
                double w_max, double w_acc, double settle_time, 
                bool is_SlowDeep, int& num_tiles, int init_pos=-1);

std::tuple<int, int, double> buildGraphOrienteeringDummyRoot(const std::string& tilefile, 
            const std::string& pro_map_file,
            std::vector<std::vector<double>>& costs, 
            std::vector<double>& probability, std::vector<int>& ranks, 
            std::vector<double>& dwell_times, double dwell_zenith,
            double w_max, double w_acc, double settle_time, 
            bool is_SlowDeep, int& num_tiles);

                
std::tuple<std::vector<int>, int, double>
buildGraphTeamOrienteering(const std::string& tilefile,
                           const std::string& pro_map_file,
                           std::vector<std::vector<double>>& costs,
                           std::vector<double>& probability,
                           std::vector<int>& ranks,
                           std::vector<double>& dwell_times,
                           std::vector<double>& ra,
                           std::vector<double>& dec,
                           double dwell_zenith,
                           double w_max, double w_acc, double settle_time,
                           bool is_SlowDeep, int m, int& num_tiles,
                           std::optional<std::pair<double, double>> zenith_pos = std::nullopt,
                           std::vector<std::pair<double, double>> init_poses = {});

// std::tuple<std::vector<int>, int, double> buildGraphTeamOrienteering(const std::string& tilefile,
//                                 const std::string& pro_map_file,
//                                 std::vector<std::vector<double>>& costs,
//                                 std::vector<double>& probability,
//                                 std::vector<int>& ranks,
//                                 std::vector<double>& dwell_times,
//                                 double dwell_zenith,
//                                 double w_max, double w_acc, double settle_time,
//                                 bool is_SlowDeep, int m, int& num_tiles,
//                                 std::vector<int> init_poses,
//                                 std::optional<std::pair<double, double>> zenith_pos= std::nullopt);


std::pair<int, int> removeNodesOrienteering(const std::vector<double>& prizes,
                const std::vector<std::vector<double>>& costs,
                const std::vector<double>& dwell_times,
                const std::vector<int>& nodes_to_remove,
                std::vector<std::vector<double>>& costs_removed,
                std::vector<double>& prizes_removed,
                std::vector<int>& ranks_mapping,
                int start_idx, int end_idx);


std::pair<std::vector<int>, int> removeVisitedNodes(std::vector<double>& prizes,
                   std::vector<std::vector<double>>& costs,
                   std::vector<int>& ranks,
                   const std::vector<int>& nodes_to_remove,
                   std::vector<int>& start_indices,
                   int end_idx);


std::vector<int> buildGraphMultiStarts(const std::string& tilefile,
                                const std::string& pro_map_file,
                                std::vector<std::vector<double>>& costs,
                                std::vector<double>& probability,
                                std::vector<int>& ranks,
                                std::vector<double>& dwell_times,
                                std::vector<double>& ra,
                                std::vector<double>& dec,
                                double dwell_zenith,
                                double w_max, double w_acc, double settle_time,
                                bool is_SlowDeep, int m, int& num_tiles,
                                std::optional<std::pair<double, double>> zenith_pos = std::nullopt,
                                std::vector<std::pair<double, double>> init_poses = {});


std::tuple<int, int, double> buildGraphwithPixelsOP(const std::string& tilefile,
                                        const std::string& pro_map_file,
                                        std::vector<std::vector<double>>& costs,
                                        std::vector<double>& dwell_times,
                                        std::vector<double>& tile_probs,
                                        std::vector<int>& ranks,
                                        std::vector<std::vector<int>>& member_pixels_trimmed,
                                        std::vector<double>& pixel_probs_trimmed,
                                        double dwell_zenith,
                                        double w_max,
                                        double w_acc,
                                        double settle_time,
                                        bool is_SlowDeep,
                                        int& num_tiles,
                                        int init_pos); 
                                        
                                        
int buildGraphwithPixelsOriginal(const std::string& tilefile, 
                                const std::string& pro_map_file,
                                std::vector<std::vector<double>>& slew_times, 
                                std::vector<double>& dwell_times,
                                std::vector<double>& tile_probs, std::vector<int>& ranks, 
                                std::vector<std::vector<int>>& member_pixels_trimmed,
                                std::vector<double>& pixel_probs_trimmed,
                                double w_max, double w_acc, double settle_time, 
                                bool is_SlowDeep, int init_pos);                                        