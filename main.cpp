
#include "max_probability_cover.h"
#include "min_telescope_cover.h"
#include "BuildGraph.h"

#include <chrono>
#include <fstream>
#include <filesystem> 
#include <numeric> 
#include <iostream>



// // void test_mmktc (std::string file, std::string out_file, double budget, 
// //     double w_max, double w_acc, double settle_time, double dwell_time,
// //     bool is_deepslow, int m) {

// //     std::vector<std::vector<double>> costs;
// //     std::vector<double> probability;
// //     std::vector<int> ranks;
// //     std::vector<double> dwell_times;
   
// //     double time_limit = 18000;
// //     double accu_thr = 0.001;
// //     int init_pos_idx = 0;
//     bool is_valid = false;

//     // auto [start_indices, end_idx, padding] = buildGraphOrienteeringMultiS(file, costs, probability, ranks, dwell_times, slew_rate, m);
//     auto [start_indices, end_idx, padding] = buildGraphOrienteeringMultiS(file, costs, probability, ranks, dwell_times, dwell_time, w_max, w_acc, settle_time, is_deepslow, init_pos_idx);
    
//     std::cout << "padding: " << padding << "\n";
//     double B = (budget + padding) * m + m;
//     std::chrono::high_resolution_clock::time_point start;
//     std::chrono::high_resolution_clock::time_point end;
//     std::chrono::duration<double> elapsed_seconds;
                            
//     std::cout << "*********running mmktc*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     auto [paths_mmktc, prize_mmktc] = mmktc_wraper(costs, probability, budget + padding, start_indices, end_idx, m);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_mmktc, start_indices, end_idx, costs.size(), budget + padding);
//     elapsed_seconds = end - start;
//     double p = multiPathsSumPrize(probability, paths_mmktc);
//     if(std::abs(p-prize_mmktc)>10e9) {
//         std::cout << "wrong prize get\n";
//         prize_mmktc = p;
//     }
//     std::cout << "total prize = " << prize_mmktc << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "MMkTC", file, paths_mmktc, ranks, budget, slew_rate, m, prize_mmktc, elapsed_seconds.count(), is_valid);


//     std::cout << "*********running Path Cover*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     auto [paths_pc, prize_pc] = PathCoverWraper(costs, probability, budget + padding, start_indices, end_idx, m);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths_any_order(costs, paths_pc, start_indices, end_idx, costs.size(), budget + padding);
//     elapsed_seconds = end - start;
//     double pc = multiPathsSumPrize(probability, paths_pc);
//     if(std::abs(pc-prize_pc)>10e9) {
//         std::cout << "wrong prize get\n";
//         prize_pc = pc;
//     }
//     std::cout << "total prize = " << prize_pc << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "PathCover", file, paths_pc, ranks, budget, slew_rate, m, prize_pc, elapsed_seconds.count(), is_valid);


//     // GCP
//     std::cout << "*********running GCP2*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     auto [paths_GCP2, prize_GCP2] = GCP2(costs, probability, budget + padding, start_indices, end_idx, m);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_GCP2, start_indices, end_idx, costs.size(), budget + padding);
//     // print_path(costs, probability, ranks, mst_pathH, padding);
//     elapsed_seconds = end - start;
//     double p2 = multiPathsSumPrize(probability, paths_GCP2);
//     if(std::abs(p2 - prize_GCP2) > 10e9) {
//         std::cout << "wrong prize get\n";
//         prize_GCP2 = p2;
//     }
//     std::cout << "total prize = " << prize_GCP2 << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "MultiGCP", file, paths_GCP2, ranks, budget, slew_rate, m, prize_GCP2, elapsed_seconds.count(), is_valid);

//     // path relinking
//     std::cout << "*********running path relinking*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     auto [paths_pr, prize_pr] = PROptimization_parallel(costs, probability, budget + padding, m, start_indices, end_idx);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_pr, start_indices, end_idx, costs.size(), budget + padding);
//     elapsed_seconds = end - start;

//     double p3 = multiPathsSumPrize(probability, paths_pr);
//     if(std::abs(p3 - prize_pr) > 10e9) {
//         std::cout << "wrong prize get\n";
//         prize_pr = p3;
//     }
//     std::cout << "total prize = " << prize_pr << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "PathRelink", file, paths_pr, ranks, budget, slew_rate, m, prize_pr, elapsed_seconds.count(), is_valid);


//     // genetic
//     std::cout << "*********running genetic*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int>  genetic_path1 = genetic_optimization_st(costs, probability, B, start_indices[0], end_idx);
//     auto [paths_genetic1, prize_genetic1] = greedy_split_prune_until_feasible(costs, probability, start_indices, end_idx, budget + padding, genetic_path1, m);
    
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_genetic1, start_indices, end_idx, costs.size(), budget + padding);
//     elapsed_seconds = end - start;

//     std::cout << "total prize = " << prize_genetic1 << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "genetic", file, paths_genetic1, ranks, budget, slew_rate, m, prize_genetic1, elapsed_seconds.count(), is_valid);


//     // annealing
//     std::cout << "*********running annealing*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     // genetic_optimization(costs, budget, 0);
//     std::vector<int>  annealing_path = simulated_annealing_optimization(costs, probability, B, start_indices[0], end_idx);
//     auto [paths_annealing, prize_annealing] = greedy_split_prune_until_feasible(costs, probability, start_indices, end_idx, budget + padding, annealing_path, m);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_annealing, start_indices, end_idx, costs.size(), budget + padding);
//     // print_path(costs, probability, ranks, annealing_path, padding);
    
//     elapsed_seconds = end - start;
//     std::cout << "total prize = " << prize_annealing << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "annealing", file, paths_annealing, ranks, budget, slew_rate, m, prize_annealing, elapsed_seconds.count(), is_valid);


//     // pso
//     std::cout << "*********running top annealing*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     // genetic_optimization(costs, budget, 0);
//     auto [paths_top_annealing, prize_top_annealing] = TopAnnealing(costs, probability, budget + padding, m, start_indices, end_idx);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_top_annealing, start_indices, end_idx, costs.size(), budget + padding);
//     // print_path(costs, probability, ranks, pso_path, padding);
    
//     elapsed_seconds = end - start;
//     std::cout << "total prize = " << prize_top_annealing << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "pso", file, paths_top_annealing, ranks, budget, slew_rate, m, prize_top_annealing, elapsed_seconds.count(), is_valid);



//     // // pso
//     // std::cout << "*********running pso*********" << std::endl;
//     // start = std::chrono::high_resolution_clock::now();
//     // // genetic_optimization(costs, budget, 0);
//     // auto [paths_pso, prize_pso] = run_pso_team_op_multi(costs, probability, budget + padding, m, start_indices, end_idx);
//     // end = std::chrono::high_resolution_clock::now();
//     // is_valid = verify_paths(costs, paths_pso, start_indices, end_idx, costs.size(), budget + padding);
//     // // print_path(costs, probability, ranks, pso_path, padding);
    
//     // elapsed_seconds = end - start;
//     // std::cout << "total prize = " << prize_pso << "\n";
//     // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     // save_result_with_paths(out_file, "pso", file, paths_pso, ranks, budget, slew_rate, m, prize_pso, elapsed_seconds.count(), is_valid);


//     // GCP
//     std::cout << "*********running GCP*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> GCP_path = GCP(costs, probability, B, start_indices[0], end_idx);
//     remove_st(GCP_path);
//     auto [paths_GCP, prize_GCP] = greedy_split_prune_until_feasible(costs, probability, start_indices, end_idx, budget + padding, GCP_path, m);
//     end = std::chrono::high_resolution_clock::now();
//     is_valid = verify_paths(costs, paths_GCP, start_indices, end_idx, costs.size(), budget + padding);
//     // print_path(costs, probability, ranks, mst_pathH, padding);
//     elapsed_seconds = end - start;
//     std::cout << "total prize = " << prize_GCP << "\n";
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     save_result_with_paths(out_file, "GCP", file, paths_GCP, ranks, budget, slew_rate, m, prize_GCP, elapsed_seconds.count(), is_valid);

//     // std::cout << "*********running gurobi*********" << std::endl;
//     // start = std::chrono::high_resolution_clock::now();
//     // auto [routes_gurobi, total_prize_gurobi] = gurobiSolveTeamSTBound(costs, probability, start_indices, end_idx, 
//     //                                                             budget + padding, m, accu_thr, time_limit, {});
//     // end = std::chrono::high_resolution_clock::now();
//     // elapsed_seconds = end - start;
//     // is_valid = verify_paths(costs, routes_gurobi, start_indices, end_idx, costs.size(), budget + padding);
    
//     // // double p4 = multiPathsSumPrize(prizes,routes_gurobi);
//     // // if(p4 != prize_pr + 10e9) {
//     // //     std::cout << "wrong prize get\n";
//     // //     prize_pr = p4;
//     // // }
    
//     // std::cout << "gurobi best total prize = " << total_prize_gurobi << "\n";
//     // std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     // // save_result(out_file, "Gurobi", file, budget, slew_rate, m, total_prize_gurobi, elapsed_seconds.count(), is_valid);
//     // save_result_with_paths(out_file, "Gurobi", file, routes_gurobi, ranks, budget,slew_rate, m, total_prize_gurobi, elapsed_seconds.count(), is_valid);
                 
// }





// int main(int argc, char** argv) {
//     auto total_start = std::chrono::high_resolution_clock::now();

//     std::string out_file = "../results2/out.csv";
//     std::string tilefile = "../data/tilings/4.0x2.0_tiling.csv";
//     std::string mapfile = "../data/LIGO_maps/txt_files_4.0x2.0/medium_maps/GW191103_012549_638.txt";

//     double w_max = 10;
//     double w_acc = 10;
//     double dwell_zenith = 2;
//     double settle_time = 0.0;
//     bool is_deepslow = false;
//     double budget = 90;
//     double mipGap = 0.001, timeLimit = 7200;

//     if (argc > 2) {
//         mapfile = std::string(argv[1]);
//         tilefile = std::string(argv[2]);
//     }
//     if (argc > 7) {
//         budget = std::stod(argv[3]);
//         w_max = std::stod(argv[4]);
//         w_acc = std::stod(argv[5]);
//         dwell_zenith = std::stod(argv[6]);
//         is_deepslow = std::stoi(argv[7]);
//     }

//     std::cout << "Input Parameters:\n";
//     std::cout << "  file        = " << mapfile << "\n";
//     std::cout << "  slew_rate   = " << w_max << "\n";
//     std::cout << "  dwell_zenith= " << dwell_zenith << "\n\n";

//     int m = 0;
//     // run_mink_unrooted(tilefile, mapfile, out_file, budget, dwell_zenith,
//     //                   w_max, w_acc, settle_time, is_deepslow, m,
//     //                   mipGap, timeLimit);

//     // mink_annealing_unrooted(tilefile, mapfile, out_file, budget, dwell_zenith,
//     //                   w_max, w_acc, settle_time, is_deepslow);

//     if(budget < 19)
//         m = 1500;
//     else if(budget < 49)
//         m = 400;
//     else if(budget < 99)
//         m = 200;
//     else if(budget < 199)
//         m = 100;
//     else if(budget < 399)
//         m = 50;
//     else if(budget < 799)
//         m = 20;
//     else m = 10;
//     std::vector<std::pair<double, double>> starts = std::move(run_mink_rooted(tilefile, mapfile, out_file, budget, dwell_zenith,
//                     w_max, w_acc, settle_time, is_deepslow, m,
//                     mipGap, timeLimit));

//     mink_annealing_rooted_random(tilefile, mapfile, out_file, budget, dwell_zenith,
//                     w_max, w_acc, settle_time, is_deepslow, m, starts);

//     // mink_annealing_rooted_best(tilefile, mapfile, out_file, budget, dwell_zenith,
//     //                 w_max, w_acc, settle_time, is_deepslow, m, starts);


//     auto total_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed_seconds = total_end - total_start;
//     std::cout << "total time (wallclock): " << elapsed_seconds.count()
//               << " seconds\n";

//     return 0;
// }


int main(int argc, char** argv) {
    auto total_start = std::chrono::high_resolution_clock::now();
    // default parameters
    std::string out_file = "../results5/out.csv";
    std::string tilefile = "../data/tilings_small/6.9x6.9_tiling.csv";
    // std::string mapfile = "../data/ADAPT_maps_new/test/adapt_p0_a0_51_map.h5";
    // std::string mapfile = "../data/LIGO_maps/csv_files/GW191103_012549_pixels.csv";
    std::string mapfile = "../data/LIGO_maps/txt_files_6.9x6.9/small_maps/GW200316_215756_50.txt";
    // std::string runtime_file = "../results/time.csv";

    int source_tile = -2;
    double w_max = 10;
    double w_acc = 10;
    double dwell_zenith = 2;
    double settle_time = 0.0;
    bool is_deepslow = false;
    double budget = 60;
    int m = 4;

    double mipGap = 0.001, timeLimit = 36000;

    if (argc > 2) {
        mapfile = std::string(argv[1]);
        tilefile = std::string(argv[2]);
    }
    if (argc > 8) {
        budget = std::stod(argv[3]);
        w_max = std::stod(argv[4]);
        w_acc = std::stod(argv[5]);
        dwell_zenith = std::stod(argv[6]);
        is_deepslow = std::stoi(argv[7]);
        m = std::stoi(argv[8]);
    }

    std::cout << "Input Parameters:\n";
    std::cout << "  file        = " << mapfile << "\n";
    std::cout << "  source_tile = " << source_tile << "\n";
    std::cout << "  slew_rate   = " << w_max << "\n";
    std::cout << "  dwell_time  = " << dwell_zenith << "\n\n";

    // run_greedy_multistarts_gcp(tilefile, mapfile, out_file,
    //                         budget, dwell_zenith, w_max, w_acc, settle_time,
    //                         is_deepslow, m);

    run_greedy_multistarts_annealing(tilefile, mapfile, out_file,
                        budget, dwell_zenith, w_max, w_acc, settle_time,
                        is_deepslow, m);

    // run_greedy_multistarts_ilp(tilefile, mapfile, out_file,
    //                     budget, dwell_zenith, w_max, w_acc, settle_time,
    //                     is_deepslow, m, mipGap, timeLimit);


    // // double mipGapOP = 0.01, timeLimitOP = 1800;

    run_top(tilefile, mapfile, out_file,
                            budget, dwell_zenith, w_max, w_acc, settle_time,
                            is_deepslow, m, mipGap, timeLimit);


    // // test
    // std::vector<std::vector<double>> slew_costs;
    // std::vector<double> probability;
    // std::vector<int> ranks;
    // std::vector<double> dwell_times;

    // int num_tiles;
    // int k = 10;
    // auto start_indices = buildGraphMultiStarts(
    //     tilefile, mapfile, slew_costs, probability, ranks, dwell_times,
    //     dwell_zenith, w_max, w_acc, settle_time, is_deepslow, k, num_tiles);

    // std::cout << "num_tiles:" << num_tiles << "\n";
    // std::cout << "dwell_times.size():" << dwell_times.size() << "\n";



    // int num_tiles2;
    // std::vector<std::vector<double>> costs;
    // std::vector<double> probability2;
    // std::vector<int> ranks2;
    // std::vector<double> dwell_times2;
    // auto [start_indices_, end_idx, padding] = buildGraphTeamOrienteering(
    //     tilefile, mapfile, costs, probability2, ranks2, dwell_times2,
    //     dwell_zenith, w_max, w_acc, settle_time, is_deepslow, m, num_tiles2);
    // std::cout << "num_tiles2:" << num_tiles2 << "\n";
    // std::cout << "dwell_times2.size():" << dwell_times2.size() << "\n";


    auto total_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_seconds = total_end - total_start;
    std::cout << "total time (wallclock): " << elapsed_seconds.count() << " seconds" << std::endl;

    
    return 0;
}



// //wrapper function on small instances test
// //test algorihms: greedy, genetic, GCP, gurobi
// void test_paths_cover (std::string tilefile, std::string mapfile, std::string out_file, 
//                     std::string overlap, 
//                     double budget, double dwell_time, 
//                     double w_max, double w_acc, double settle_time,
//                     bool is_deepslow) {

//     std::vector<std::vector<double>> slew_times;
//     std::vector<double> tile_probs;
//     std::vector<int> ranks;
//     std::vector<double> dwell_times;
//     std::vector<std::vector<int>> member_pixels;
//     std::vector<double> pixel_probs;

//     std::string file = mapfile;
   
//     double time_limit = 3600;
//     double accu_thr = 0.000001;
//     int init_pos_idx = 0;


//     int start_idx = buildGraph(tilefile, mapfile, slew_times, dwell_times, tile_probs, ranks, member_pixels, pixel_probs,
//                         w_max, w_acc, settle_time, is_deepslow, 0);

//     std::chrono::high_resolution_clock::time_point start;
//     std::chrono::high_resolution_clock::time_point end;
//     std::chrono::duration<double> elapsed_seconds;

//     // std::vector<int> path = gurobiPathCover(slew_times, dwell_times, pixel_probs,member_pixels, start_idx, budget, accu_thr, time_limit)
    
//     // gurobi with time limit
//     std::cout << "*********ILP Path Cover*********" << std::endl;
//     start = std::chrono::high_resolution_clock::now();
//     std::vector<int> path = gurobiPathCover(slew_times, dwell_times, pixel_probs,member_pixels, start_idx, budget, accu_thr, time_limit);
    
//     end = std::chrono::high_resolution_clock::now();
//     print_path_overlap(slew_times, dwell_times, member_pixels, pixel_probs, ranks, path);
//     elapsed_seconds = end - start;
//     std::cout << "running time (wallclock): " << elapsed_seconds.count() << "seconds" << std::endl;
//     // save_result(out_file, "Gurobi_Flow", file, budget, w_max,costs, probability, ranks, ilp_path2, elapsed_seconds.count(), padding);
//     save_result_overlap(out_file, overlap+"path+cover", tilefile, mapfile, budget, w_max, slew_times, dwell_times, member_pixels, pixel_probs, ranks, path, elapsed_seconds.count());
    
// }
