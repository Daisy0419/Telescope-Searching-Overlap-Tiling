
#include "helplers.h"
#include "PathOperations.h"

#include <iostream>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <fstream>
#include <filesystem>
#include <vector>
#include <numeric>
#include <iomanip>
#include <unordered_map>
#include <limits>
#include <cmath>
#include <stdexcept>


double random_double(double min, double max) {
    thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<> dis(min, max);
    return dis(gen);
}

int random_int(int min, int max) {
    thread_local std::mt19937 gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(min, max);
    return dis(gen);
}

std::vector<int> unique_random_ints(int min, int max, int n) {
    if (n > (max - min + 1)) {
        throw std::invalid_argument("Cannot generate more unique numbers than the range size.");
    }
    std::vector<int> numbers(max - min + 1);
    std::iota(numbers.begin(), numbers.end(), min);

    thread_local std::mt19937 gen(std::random_device{}());
    std::shuffle(numbers.begin(), numbers.end(), gen);

    return std::vector<int>(numbers.begin(), numbers.begin() + n);
}

// double random_double(double min, double max) {
//     static std::random_device rd;
//     static std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dis(min, max);
//     return dis(gen);
// }

// int random_int(int min, int max) {
//     static std::random_device rd;
//     static std::mt19937 gen(rd());
//     std::uniform_real_distribution<> dis(min, max);
//     return dis(gen);
// }

// std::vector<int> unique_random_ints(int min, int max, int n) {
//     if (n > (max - min + 1)) {
//         throw std::invalid_argument("Cannot generate more unique numbers than the range size.");
//     }
//     std::vector<int> numbers(max - min + 1);
//     std::iota(numbers.begin(), numbers.end(), min);

//     static std::random_device rd;
//     static std::mt19937 gen(rd());
//     std::shuffle(numbers.begin(), numbers.end(), gen);

//     return std::vector<int>(numbers.begin(), numbers.begin() + n);
// }


// 2-opt, checks all valid pairs of edges in the path
void fix_cross(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    if (costs.empty() || costs.size() != costs[0].size()) {
        throw std::invalid_argument("Cost matrix must be non-empty and square.");
    }
    if (path.size() > costs.size()) {
        throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
    }
    if (path.size() < 4) { // at least 4 nodes
        return;  
    }

    for (size_t i = 1; i < path.size() - 2; ++i) {
        for (size_t j = i + 1; j < path.size() - 2; ++j) {
            int a = path[i - 1]; 
            int b = path[i];   
            int c = path[j];    
            int d = path[j + 1]; 

            double original_cost = costs[a][b] + costs[c][d];
            double uncrossed_cost = costs[a][c] + costs[b][d];

            if (uncrossed_cost < original_cost - 1e-9) {
                std::reverse(path.begin() + i, path.begin() + j + 1);
            }
        }
    }
}

// 2-opt, ensures that the path remains a valid open path from fixed start to fixed end
void fix_cross_st_path(std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    if (costs.empty() || costs.size() != costs[0].size()) {
        throw std::invalid_argument("Cost matrix must be non-empty and square.");
    }
    if (path.size() > costs.size()) {
        throw std::invalid_argument("Path contains nodes outside cost matrix bounds.");
    }
    if (path.size() < 4) {
        return;  // nothing to optimize
    }

    // Avoid modifying edges adjacent to start (s) or end (t)
    for (size_t i = 1; i < path.size() - 2; ++i) {
        for (size_t j = i + 1; j < path.size() - 1; ++j) {
            int a = path[i - 1];
            int b = path[i];
            int c = path[j];
            int d = path[j + 1];

            if (i == 1 && j + 1 == path.size() - 1) continue; // affect both s and t
            if (i == 1 && a == path[0]) continue;             // affects s
            if (j + 1 == path.size() - 1 && d == path.back()) continue; // affects t

            double original_cost = costs[a][b] + costs[c][d];
            double uncrossed_cost = costs[a][c] + costs[b][d];

            if (uncrossed_cost < original_cost - 1e-9) {
                std::reverse(path.begin() + i, path.begin() + j + 1);
            }
        }
    }
}


// compute total cost of a given path
double calculate_path_cost(const std::vector<int>& path, const std::vector<std::vector<double>>& costs) {
    double total_cost = 0.0;
    for (size_t i = 1; i < path.size(); ++i) {
        total_cost += costs[path[i - 1]][path[i]];
    }
    return total_cost;
}

// compute total cost of a given path
double calculate_paths_prize(const std::vector<std::vector<int>>& paths, 
                             const std::vector<double>& prizes) {

    double total_prize = 0.0;
    for (auto path : paths) 
        for (size_t i = 1; i < path.size()-1; ++i) {
            total_prize += prizes[path[i]];
        }
    return total_prize;
}


double calculate_paths_pixel_prize(
    const std::vector<std::vector<int>>& paths,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels
) {
    return pathsPixelPrize(pixel_probs, member_pixels, paths);
}



void print_paths (const std::vector<std::vector<int>>& paths) {
    std::cout << "Printing paths: " << std::endl;
    for(auto path : paths) {
        for (auto p : path) {
            std::cout << p << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

/****************************** */

void save_result(const std::string& filename, 
                  const std::string& method, 
                  const std::string& data, 
                  double budget, 
                  double slew_rate, 
                  int npaths,
                  const double prize,
                  double elapsed_time,
                  bool is_valid) {

    bool file_exists = std::filesystem::exists(filename); 
    std::ofstream outfile(filename, std::ios::app);  

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    //Write header
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,SlewRate,nPath, SumProb,TimeSec,IsValid\n";
    }

    // Write result to CSV
    outfile << method << ","
            << data << ","
            << budget << ","
            << slew_rate << ","
            << npaths << ","
            << prize << ","
            << elapsed_time << ",";

    if(is_valid)
        outfile << "true";
    else    
        outfile << "false";
    outfile << "\n";
    outfile.close();

}

static std::string csv_escape(const std::string& s) {
    bool need_quotes = s.find_first_of(",\"\n\r") != std::string::npos;
    if (!need_quotes) return s;
    std::string out; out.reserve(s.size() + 2);
    out.push_back('"');
    for (char c : s) out.append(c == '"' ? "\"\"" : std::string(1, c));
    out.push_back('"');
    return out;
}

void save_result_with_paths(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            const std::vector<std::vector<int>>& paths,
                            const std::vector<int>& rank,
                            double budget, 
                            double w_max, double w_acc,
                            int npaths,
                            int nTiles,
                            const double prize,
                            const double best_prize,
                            double tiling_time,
                            double planning_time,
                            bool IsValid) {
    bool file_exists = std::filesystem::exists(filename);
    std::ofstream outfile(filename, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    // Header (no trailing comma)
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,w_max,w_acc,nPath,mapTiles,SumProb,SumProbBound,TilingTime,PlanningTime,IsValid";
        for (size_t i = 0; i < paths.size(); ++i) {
            outfile << ",path " << i;
        }
        outfile << "\n";
    }

    // Fixed columns
    outfile << csv_escape(method) << ","
            << std::filesystem::path(data).filename().string() << ","
            << budget             << ","
            << w_max          << ","
            << w_acc          << ","
            << npaths             << ","
            << nTiles             << ","
            << prize              << ","
            << best_prize              << ","
            << tiling_time    << ","   
            << planning_time    << ","    
            << (IsValid ? "true" : "false");

    // Each path as a single CSV field containing space-separated node ids (rank-mapped)
    for (const auto& path : paths) {
        std::ostringstream ps;
        for (size_t i = 0; i < path.size(); ++i) {
            int idx = path[i];
            int out_id = (idx >= 0 && idx < (int)rank.size()) ? rank[idx] : idx; // guard
            if (i) ps << ' ';
            ps << out_id;
        }
        outfile << "," << csv_escape(ps.str());
    }

    outfile << "\n";
    // ofstream dtor closes & flushes
}


void save_result_with_paths2(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            const std::vector<std::vector<int>>& paths,
                            double budget, 
                            double w_max, double w_acc,
                            int npaths,
                            int nTiles,
                            const double prize,
                            const double best_prize,
                            double tiling_time,
                            double planning_time,
                            bool IsValid) {

    bool file_exists = std::filesystem::exists(filename);
    std::ofstream outfile(filename, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    // Header (no trailing comma)
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,w_max,w_acc,nPath,mapTiles,SumProb,SumProbBound,TilingTime,PlanningTime,IsValid";
        for (size_t i = 0; i < paths.size(); ++i) {
            outfile << ",path " << i;
        }
        outfile << "\n";
    }

    // Fixed columns
    outfile << csv_escape(method) << ","
            << std::filesystem::path(data).filename().string() << ","
            << budget             << ","
            << w_max          << ","
            << w_acc          << ","
            << npaths             << ","
            << nTiles             << ","
            << prize              << ","
            << best_prize              << ","
            << tiling_time    << ","   
            << planning_time    << ","    
            << (IsValid ? "true" : "false");

    // Each path as a single CSV field containing space-separated node ids (rank-mapped)
    for (const auto& path : paths) {
        std::ostringstream ps;
        for (size_t i = 0; i < path.size(); ++i) {
            int out_id = path[i];
            if (i) ps << ' ';
            ps << out_id;
        }
        outfile << "," << csv_escape(ps.str());
    }

    outfile << "\n";
    // ofstream dtor closes & flushes
}

void save_result_without_paths(const std::string& filename, 
                            const std::string& method, 
                            const std::string& data, 
                            double budget, 
                            double w_max, double w_acc,
                            double num_tiles,
                            double covered_tiles,
                            double npaths,
                            double npaths_bound,
                            double tiling_time,
                            double planning_time,
                            const std::vector<std::vector<int>>& paths,
                            bool IsValid) {

    bool file_exists = std::filesystem::exists(filename);
    std::ofstream outfile(filename, std::ios::app);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    // Header (no trailing comma)
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,w_max,w_acc,nTiles,CoveredTiles,nPath,nPathBound,TilingTime,PlanningTime,IsValid\n";
    }

    // Fixed columns
    outfile << csv_escape(method) << ","
            << std::filesystem::path(data).filename().string() << ","
            << budget             << ","
            << w_max          << ","
            << w_acc          << ","
            << num_tiles             << ","
            << covered_tiles             << ","
            << npaths             << ","
            << npaths_bound             << ","
            << tiling_time    << ","   
            << planning_time    << ","    
            << (IsValid ? "true" : "false");

    // Each path as a single CSV field containing space-separated node ids (rank-mapped)
    for (const auto& path : paths) {
        std::ostringstream ps;
        for (size_t i = 0; i < path.size(); ++i) {
            int out_id = path[i];
            if (i) ps << ' ';
            ps << out_id;
        }
        outfile << "," << csv_escape(ps.str());
    }

    outfile << "\n";
    // ofstream dtor closes & flushes
}

//print path to stdout
double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
                    std::vector<int> ranks, std::vector<int> best_path, double padding) {

    if (best_path.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }

    if (ranks.empty()) {
        ranks.resize(best_path.size());
        std::iota(ranks.begin(), ranks.end(), 1); 
    }
    
    double sum_probability = 0.0;
    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << ranks[tile] << " "; 
        // std::cout << rank[tile] << " "; //ranks in tiling file
        sum_probability += probability[tile]; 
    }
    std::cout << std::endl;

    double total_cost = -padding;
    for (size_t i = 0; i < best_path.size() - 1; i++) {
        total_cost += costs[best_path[i]][best_path[i + 1]]; 
    }
    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; 
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}


double print_path_pixels(
    const std::vector<std::vector<double>>& costs,
    const std::vector<double>& pixel_probs,
    const std::vector<std::vector<int>>& member_pixels,
    std::vector<int> ranks,
    const std::vector<int>& best_path,
    double padding ) {
    if (best_path.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }
    validatePixelCoverageData(
        costs.size(), pixel_probs, member_pixels);
    if (ranks.empty()) {
        ranks.resize(costs.size());
        std::iota(ranks.begin(), ranks.end(), 0);
    }
    if (ranks.size() != costs.size()) {
        throw std::invalid_argument(
            "print_path_pixels: ranks must contain one value per node.");
    }

    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        if (tile < 0 || static_cast<std::size_t>(tile) >= costs.size()) {
            throw std::out_of_range(
                "print_path_pixels: invalid node in path.");
        }
        std::cout << ranks[tile] << " ";
    }
    std::cout << std::endl;

    const double union_probability = pathPixelPrize(
        pixel_probs, member_pixels, best_path);
    const double total_cost = pathCost(costs, best_path) - padding;
    std::cout << "Num Tiles in Path: "
              << (best_path.size() >= 2 ? best_path.size() - 2 : 0)
              << std::endl;
    std::cout << "Union Pixel Probability: "
              << union_probability << std::endl;
    std::cout << "Total Cost: " << total_cost << std::endl;
    return union_probability;
}


void remove_st (std::vector<int>& path) {
    if (path.size() >= 2) {
        path.erase(path.begin()); // O(n)
        path.pop_back();     
    }
}


bool verify_routes(const Instance& I,
                   const std::vector<std::vector<int>>& routes,
                   double Budget,
                   int m_expected,                 // −1 don’t enforce
                   double tol) {
    const int N   = I.size;
    const int s   = I.s;
    const int t   = I.t;

    bool ok = true;

    if (m_expected >= 0 && (int)routes.size() != m_expected) {
        std::cerr << "[verify] expected " << m_expected
                  << " routes, got " << routes.size() << "\n";
        ok = false;
    }

    std::unordered_set<int> seen_global;   // all intermediate nodes used

    auto route_cost = [&](const std::vector<int>& r)->double{
        double c=0.0;
        for (size_t i=0;i+1<r.size();++i) c += I.costs[r[i]][r[i+1]];
        return c;
    };

    for (size_t r_idx = 0; r_idx < routes.size(); ++r_idx) {
        const auto& path = routes[r_idx];

        // ---- structural checks ----
        if (path.empty()) {
            std::cerr << "[verify] route " << r_idx << " is empty\n";
            ok = false; continue;
        }
        if (path.front() != s || path.back() != t) {
            std::cerr << "[verify] route " << r_idx
                      << " does not start at s or end at t\n";
            ok = false;
        }

        // ---- budget check ----
        double cost = route_cost(path);
        if (cost > Budget + tol) {
            std::cerr << "[verify] route " << r_idx << " cost "
                      << cost << " exceeds Budget " << Budget << "\n";
            ok = false;
        }

        // ---- duplicate checks ----
        std::unordered_set<int> seen_local;
        for (size_t pos = 0; pos < path.size(); ++pos) {
            int v = path[pos];
            if (v < 0 || v >= N) {
                std::cerr << "[verify] node id " << v
                          << " out of range in route " << r_idx << "\n";
                ok = false;
            }
            if (v == s || v == t) continue;

            if (!seen_local.insert(v).second) {               // appears twice in same route
                std::cerr << "[verify] node " << v
                          << " appears twice in route " << r_idx << "\n";
                ok = false;
            }
        }
        for (int v : seen_local) {
            if (!seen_global.insert(v).second) {              // appears in another route
                std::cerr << "[verify] node " << v
                          << " appears in multiple routes\n";
                ok = false;
            }
        }
    }

    if (ok)
        std::cout << "[verify] all checks passed for "
                  << routes.size() << " route(s)\n";
    return ok;
}


bool verify_paths(const std::vector<std::vector<double>>& costs,
                  const std::vector<std::vector<int>>& routes,
                  const std::vector<int>& s,  
                  const int t, const int N,
                  double Budget,
                  int m_expected,
                  double tol) {
    bool ok = true;

    if (m_expected >= 0 && static_cast<int>(routes.size()) != m_expected) {
        std::cerr << "[verify] expected " << m_expected
                  << " routes, got " << routes.size() << "\n";
        ok = false;
    }

    if (s.size() != routes.size()) {
        std::cerr << "[verify] s.size() (" << s.size()
                  << ") does not match routes.size() (" << routes.size() << ")\n";
        ok = false;
    }

    // Validate starts (allow duplicates), and build a set to exempt them from dup checks.
    std::unordered_set<int> start_set;
    for (size_t i = 0; i < s.size(); ++i) {
        int si = s[i];
        if (si < 0 || si >= N) {
            std::cerr << "[verify] s[" << i << "]=" << si
                      << " out of range [0," << (N-1) << "]\n";
            ok = false;
        }
        if (si == t) {
            std::cerr << "[verify] s[" << i << "] equals terminal t=" << t << "\n";
            ok = false;
        }
        start_set.insert(si); // duplicates are fine; set is only for membership tests
    }

    auto route_cost = [&](const std::vector<int>& r)->double {
        double c = 0.0;
        for (size_t i = 0; i + 1 < r.size(); ++i) 
            c += costs[r[i]][r[i+1]];
        return c;
    };

    std::unordered_set<int> seen_global; // track intermediate nodes across routes

    for (size_t r_idx = 0; r_idx < routes.size(); ++r_idx) {
        const auto& path = routes[r_idx];
        const int s_i = (r_idx < s.size()) ? s[r_idx] : -1;

        // ---- structural checks ----
        if (path.empty()) {
            std::cerr << "[verify] route " << r_idx << " is empty\n";
            ok = false; 
            continue;
        }
        if (s_i < 0 || s_i >= N) {
            ok = false; // already reported above
        } else if (path.front() != s_i) {
            std::cerr << "[verify] route " << r_idx
                      << " does not start at s[" << r_idx << "]=" << s_i
                      << " (got " << path.front() << ")\n";
            ok = false;
        }
        if (path.back() != t) {
            std::cerr << "[verify] route " << r_idx
                      << " does not end at t=" << t << " (got " << path.back() << ")\n";
            ok = false;
        }

        // ---- budget check ----
        double cost = route_cost(path);
        if (cost > Budget + tol) {
            std::cerr << "[verify] route " << r_idx << " cost "
                      << cost << " exceeds Budget " << Budget << "\n";
            ok = false;
        }

        // ---- duplicate & range checks (skip any start node and t) ----
        std::unordered_set<int> seen_local;
        for (size_t pos = 0; pos < path.size(); ++pos) {
            int v = path[pos];

            if (v < 0 || v >= N) {
                std::cerr << "[verify] node id " << v
                          << " out of range in route " << r_idx << "\n";
                ok = false;
            }

            if (start_set.count(v) || v == t) continue; // exempt starts and terminal

            if (!seen_local.insert(v).second) {
                std::cerr << "[verify] node " << v
                          << " appears twice in route " << r_idx << "\n";
                ok = false;
            }
        }

        for (int v : seen_local) {
            if (!seen_global.insert(v).second) {
                std::cerr << "[verify] node " << v
                          << " appears in multiple routes\n";
                ok = false;
            }
        }
    }

    if (ok) {
        std::cout << "[verify] all checks passed for "
                  << routes.size() << " route(s)\n";
    }
    return ok;
}

bool verify_paths_any_order(const std::vector<std::vector<double>>& costs,
                            const std::vector<std::vector<int>>& routes,
                            const std::vector<int>& s,   // allowed starts (multiset semantics)
                            int t, int N,
                            double Budget,
                            int m_expected,
                            double tol) {
    bool ok = true;

    if (m_expected >= 0 && static_cast<int>(routes.size()) != m_expected) {
        std::cerr << "[verify] expected " << m_expected
                  << " routes, got " << routes.size() << "\n";
        ok = false;
    }

    // Order no longer matters; size mismatch is a warning (not a failure)
    if (s.size() != routes.size()) {
        std::cerr << "[verify] warning: s.size() (" << s.size()
                  << ") != routes.size() (" << routes.size() << ")\n";
    }

    // Validate starts and build multiset (counts) of allowable starts
    std::unordered_map<int,int> start_count;
    start_count.reserve(s.size() * 2);

    for (size_t i = 0; i < s.size(); ++i) {
        int si = s[i];
        if (si < 0 || si >= N) {
            std::cerr << "[verify] s[" << i << "]=" << si
                      << " out of range [0," << (N-1) << "]\n";
            ok = false;
        }
        if (si == t) {
            std::cerr << "[verify] s[" << i << "] equals terminal t=" << t << "\n";
            ok = false;
        }
        ++start_count[si];
    }

    // For duplicate checks we’ll exempt any start node and t.
    std::unordered_set<int> start_keys;
    start_keys.reserve(start_count.size());
    for (const auto& kv : start_count) start_keys.insert(kv.first);

    auto route_cost = [&](const std::vector<int>& r)->double {
        double c = 0.0;
        for (size_t i = 0; i + 1 < r.size(); ++i) c += costs[r[i]][r[i+1]];
        return c;
    };

    std::unordered_set<int> seen_global; // track non-start, non-t nodes across routes

    for (size_t r_idx = 0; r_idx < routes.size(); ++r_idx) {
        const auto& path = routes[r_idx];

        // ---- structural checks ----
        if (path.empty()) {
            std::cerr << "[verify] route " << r_idx << " is empty\n";
            ok = false;
            continue;
        }
        if (path.front() < 0 || path.front() >= N) {
            std::cerr << "[verify] route " << r_idx
                      << " invalid start id " << path.front() << "\n";
            ok = false;
        }
        if (path.back() != t) {
            std::cerr << "[verify] route " << r_idx
                      << " does not end at t=" << t << " (got " << path.back() << ")\n";
            ok = false;
        }

        // ---- start check (order-independent, multiset semantics) ----
        int s_node = path.front();
        auto it = start_count.find(s_node);
        if (it == start_count.end() || it->second == 0) {
            std::cerr << "[verify] route " << r_idx
                      << " starts at unexpected or overused start " << s_node << "\n";
            ok = false;
        } else {
            --(it->second); // consume one copy of this start
        }

        // ---- budget check ----
        double cst = route_cost(path);
        if (cst > Budget + tol) {
            std::cerr << "[verify] route " << r_idx << " cost "
                      << cst << " exceeds Budget " << Budget << "\n";
            ok = false;
        }

        // ---- duplicate & range checks (skip any start node and t) ----
        std::unordered_set<int> seen_local;
        for (size_t pos = 0; pos < path.size(); ++pos) {
            int v = path[pos];

            if (v < 0 || v >= N) {
                std::cerr << "[verify] node id " << v
                          << " out of range in route " << r_idx << "\n";
                ok = false;
            }

            if (start_keys.count(v) || v == t) continue; // exempt starts and terminal

            if (!seen_local.insert(v).second) {
                std::cerr << "[verify] node " << v
                          << " appears twice in route " << r_idx << "\n";
                ok = false;
            }
        }

        for (int v : seen_local) {
            if (!seen_global.insert(v).second) {
                std::cerr << "[verify] node " << v
                          << " appears in multiple routes\n";
                ok = false;
            }
        }
    }

    // Any starts left unused? Warn (not fail). Flip to an error if you require full coverage.
    for (const auto& kv : start_count) {
        if (kv.second != 0) {
            std::cerr << "[verify] warning: start " << kv.first
                      << " unused (remaining count=" << kv.second << ")\n";
            // If strict coverage is desired, uncomment:
            // ok = false;
        }
    }

    if (ok) {
        std::cout << "[verify] all checks passed for "
                  << routes.size() << " route(s) [order-agnostic]\n";
    }
    return ok;
}


// print multiple paths to stdout
// Returns: total sum of probabilities over all tiles in all paths
double print_paths(const std::vector<std::vector<double>>& costs,
                   const std::vector<double>& probability,
                   std::vector<int> rank,                      // optional; can be empty
                   const std::vector<std::vector<int>>& paths, // multiple paths
                   double padding) {
    if (paths.empty()) {
        std::cout << "No valid path found." << std::endl;
        return 0.0;
    }

    // Find max tile id to validate/size rank if needed
    int max_tile = -1;
    for (const auto& p : paths)
        for (int v : p) max_tile = std::max(max_tile, v);

    if (max_tile < 0) {
        std::cout << "No tiles in provided paths." << std::endl;
        return 0.0;
    }

    // If rank not provided, synthesize a simple 1..(max_tile+1) mapping (safe for indexing by tile id)
    if (rank.empty()) {
        rank.resize(static_cast<size_t>(max_tile + 1));
        std::iota(rank.begin(), rank.end(), 1);
    }

    // Helper to compute path cost with your padding rule
    auto path_cost = [&](const std::vector<int>& p)->double {
        double c = -padding;
        for (size_t i = 0; i + 1 < p.size(); ++i) c += costs[p[i]][p[i+1]];
        return c;
    };

    // Totals across all paths
    double total_sum_prob = 0.0;
    double total_cost     = 0.0;
    size_t total_tiles    = 0;

    // Print each path
    for (size_t i = 0; i < paths.size(); ++i) {
        const auto& path = paths[i];
        if (path.empty()) {
            std::cout << "Path " << i << ": <empty>" << std::endl;
            continue;
        }

        // Print tile ids
        std::cout << "Path " << i << " (tile id): ";
        for (int tile : path) std::cout << tile << " ";
        std::cout << "\n";

        // Optionally print ranks (only if we can index safely)
        bool rank_ok = (static_cast<int>(rank.size()) > max_tile);
        if (rank_ok) {
            std::cout << "Path " << i << " (tile rank): ";
            for (int tile : path) std::cout << rank[tile] << " ";
            std::cout << "\n";
        }

        // Sum probability for this path
        double sum_prob = 0.0;
        for (int tile : path) sum_prob += probability[tile];

        // Cost for this path
        double c = path_cost(path);

        // Per-path stats
        std::cout << "  Num Tiles: " << path.size() << "\n";
        std::cout << "  Sum Probability: " << sum_prob << "\n";
        std::cout << "  Total Cost: " << c << "\n";

        // Accumulate totals
        total_sum_prob += sum_prob;
        total_cost     += c;
        total_tiles    += path.size();
    }

    // Print totals
    std::cout << "---- Totals ----\n";
    std::cout << "Num Paths: " << paths.size() << "\n";
    std::cout << "Total Tiles: " << total_tiles << "\n";
    std::cout << "Total Sum Probability: " << total_sum_prob << "\n";
    std::cout << "Total Cost: " << total_cost << "\n";

    return total_sum_prob;
}


// only keep one start from a team orienteering graph
std::pair<int,int> extractSingleStartSubgraph(
        const std::vector<std::vector<double>>& costs_full,
        const std::vector<int>& start_indices,
        int start_to_keep,
        int end_index,
        std::vector<std::vector<double>>& costs_out) {

    const int s_full = start_to_keep;
    const int t_full = end_index;
    int num_tiles = (int)costs_full.size() - (int)start_indices.size() - 1;
    const int M = num_tiles + 2;        // tiles + 1 start + t
    const int s_new = num_tiles;
    const int t_new = num_tiles + 1;

    costs_out.assign(M, std::vector<double>(M, 0.0));

    // Tile-to-tile block
    for (int i = 0; i < num_tiles; ++i)
        for (int j = 0; j < num_tiles; ++j)
            costs_out[i][j] = costs_full[i][j];

    // Start row/column <-> tiles
    for (int i = 0; i < num_tiles; ++i) {
        costs_out[s_new][i] = costs_full[s_full][i];
        costs_out[i][s_new] = costs_full[i][s_full];
    }
    costs_out[s_new][s_new] = costs_full[s_full][s_full];

    // t row/column <-> tiles
    for (int i = 0; i < num_tiles; ++i) {
        costs_out[t_new][i] = costs_full[t_full][i];
        costs_out[i][t_new] = costs_full[i][t_full];
    }
    costs_out[t_new][t_new] = costs_full[t_full][t_full];

    // Start <-> t
    costs_out[s_new][t_new] = costs_full[s_full][t_full];
    costs_out[t_new][s_new] = costs_full[t_full][s_full];

    return {s_new, t_new};
}



std::vector<double> extractSingleStartProbability(
        const std::vector<double>& probability_full,
        const std::vector<int>& start_indices) {

    int num_tiles = (int)probability_full.size() - (int)start_indices.size() - 1;
    std::vector<double> prob_out(num_tiles + 2, 0.0);

    for (int i = 0; i < num_tiles; ++i)
        prob_out[i] = probability_full[i];

    // prob_out[num_tiles]   = 0.0  (start, already zero)
    // prob_out[num_tiles+1] = 0.0  (end,   already zero)

    return prob_out;
}



double coveredProbabilityUnion(
    const std::vector<int>& path_tiles, 
    const std::vector<std::vector<int>>& member_pixels, 
    const std::vector<double>& pixel_probs) {

    const int nTiles  = static_cast<int>(member_pixels.size());
    const int nPixels = static_cast<int>(pixel_probs.size());
    if (nPixels == 0) return 0.0;

    std::vector<char> seen(nPixels, 0);

    double total = 0.0;

    for (int tile : path_tiles) {
        if (tile < 0 || tile >= nTiles) {
            throw std::runtime_error("coveredProbabilityUnion: tile index out of range.");
        }
        for (int p : member_pixels[tile]) {
            if (p < 0 || p >= nPixels) {
                continue;
            }
            if (!seen[p]) {
                seen[p] = 1;
                total += pixel_probs[p];
            }
        }
    }

    return total;
}

double print_path_overlap(const std::vector<std::vector<double>>& slew_times, 
                        const std::vector<double>& dwell_times, 
                        const std::vector<std::vector<int>>& member_pixels, 
                        const std::vector<double>& pixel_probs, 
                        std::vector<int> rank, 
                        std::vector<int> best_path) {

    if (best_path.empty()) {
        throw std::runtime_error(
            "No valid path found. Please check your Gurobi license is mounted correctly."
        );
    }

    if (rank.empty()) {
        rank.resize(best_path.size());
        std::iota(rank.begin(), rank.end(), 1); 
    }
    
    double sum_probability = coveredProbabilityUnion(best_path, member_pixels, pixel_probs);

    double total_cost = dwell_times[best_path[0]];
    for (size_t i = 1; i < best_path.size(); i++) {
        total_cost = total_cost + slew_times[best_path[i-1]][best_path[i]]+dwell_times[best_path[i]]; 
    }

    std::cout << "Path(tile rank): ";
    for (int tile : best_path) {
        std::cout << tile << " "; 
    }
    std::cout << std::endl;

    std::cout << "Num Tiles in Path: " << best_path.size() << std::endl;
    std::cout << "Sum Probability: " << sum_probability << std::endl; 
    std::cout << "Total Cost: " << total_cost << std::endl;

    return sum_probability;
}


void save_result_overlap(const std::string& filename,
                        const std::string& method,
                        const std::string& tiling_data,
                        const std::string& likelihood_data,
                        double budget, double slew_rate,
                        const std::vector<std::vector<double>>& slew_times, 
                        const std::vector<double>& dwell_times, 
                        const std::vector<std::vector<int>>& member_pixels, 
                        const std::vector<double>& pixel_probs, 
                        std::vector<int> rank, 
                        std::vector<int> best_path,
                        double elapsed_time) {

    bool file_exists = std::filesystem::exists(filename); 

    std::ofstream outfile(filename, std::ios::app);  

    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open results file: " << filename << std::endl;
        return;
    }

    //Write header
    if (!file_exists) {
        outfile << "Method,Dataset,Budget,SlewRate,NumTiles,SumProb,TotalCost,TimeSec,Path\n";
    }

    if (best_path.empty()) {
        throw std::runtime_error(
            "No valid path found. Please check your Gurobi license is mounted correctly."
        );
    }

    if (rank.empty()) {
        rank.resize(best_path.size());
        std::iota(rank.begin(), rank.end(), 1); 
    }
    
    double sum_probability = coveredProbabilityUnion(best_path, member_pixels, pixel_probs);

    double total_cost = dwell_times[best_path[0]];
    for (size_t i = 1; i < best_path.size(); i++) {
        total_cost = total_cost + slew_times[best_path[i-1]][best_path[i]]+dwell_times[best_path[i]]; 
    }


    outfile << method << ","
            << std::filesystem::path(tiling_data).filename().string() << ","
            << std::filesystem::path(likelihood_data).filename().string() << ","
            << budget << ","
            << slew_rate << ","
            << best_path.size() << ","
            << sum_probability << ","
            << total_cost << ","
            << elapsed_time << ",";

    // Write path
    for (size_t i = 0; i < best_path.size(); i++) {
        outfile << rank[best_path[i]];
        if (i < best_path.size() - 1) {
            outfile << " ";  
        }
    }
    outfile << "\n";

    outfile.close();
}