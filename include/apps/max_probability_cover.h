#pragma once
#include <string>

void run_top(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m,
                            double accu_thr = 0.0001,
                            double time_limit = 1800);

void run_greedy_multistarts_gcp(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m);

void run_greedy_multistarts_annealing(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m);


void run_greedy_multistarts_ilp(const std::string& tilefile,
                            const std::string& mapfile,
                            const std::string& out_file,
                            double budget, double dwell_zenith,
                            double w_max, double w_acc, double settle_time,
                            bool is_deepslow, int m,
                            double mipGap = 0.0001, 
                            double timeLimit = 1800);