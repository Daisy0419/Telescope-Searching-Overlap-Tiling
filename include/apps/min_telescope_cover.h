#pragma once
#include <string>
#include <vector>

void run_mink_unrooted(const std::string& tilefile,
                       const std::string& mapfile,
                       const std::string& out_file,
                       double budget, double dwell_zenith,
                       double w_max, double w_acc, double settle_time,
                       bool is_deepslow, int m,
                       double mipGap, double timeLimit);

                       

std::vector<std::pair<double, double>> run_mink_rooted(const std::string& tilefile,
                     const std::string& mapfile,
                     const std::string& out_file,
                     double budget, double dwell_zenith,
                     double w_max, double w_acc, double settle_time,
                     bool is_deepslow, int m,
                     double mipGap, double timeLimit);
                     
void mink_annealing_rooted(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m, 
                            std::vector<std::pair<double, double>> init_poses);


void mink_annealing_rooted_best(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m,
                            std::vector<std::pair<double, double>> init_poses);

void mink_annealing_unrooted(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow);


void mink_annealing_rooted_random(const std::string& tilefile,
                           const std::string& mapfile,
                           const std::string& out_file,
                           double budget, double dwell_zenith,
                           double w_max, double w_acc, double settle_time,
                           bool is_deepslow, int m,
                           std::vector<std::pair<double, double>> init_poses);