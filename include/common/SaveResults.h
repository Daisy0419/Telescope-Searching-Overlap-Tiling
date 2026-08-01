#pragma once
#include <string> 
#include <vector>
#include <filesystem> 
#include <numeric>
#include <stdexcept> 
#include <iomanip>
#include <chrono>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <thread>
#include <atomic>
#include <unordered_map>
#include <memory>


class TimingLog {
public:
    using Clock = std::chrono::high_resolution_clock;
    
    struct Entry {
        std::string label;
        double elapsed_ms;
        double cumulative_ms;
    };

    void start() { 
        start_time_ = Clock::now(); 
        last_time_ = start_time_;
    }

    void mark(const std::string& label) {
        auto now = Clock::now();
        double elapsed = std::chrono::duration<double, std::milli>(now - last_time_).count();
        double cumulative = std::chrono::duration<double, std::milli>(now - start_time_).count();
        entries_.push_back({label, elapsed, cumulative});
        last_time_ = now;
    }

    // Get elapsed time for a specific label (returns -1 if not found)
    double get_elapsed(const std::string& label) const {
        for (const auto& e : entries_) {
            if (e.label == label) return e.elapsed_ms;
        }
        return -1.0;
    }

    // Get cumulative time for a specific label
    double get_cumulative(const std::string& label) const {
        for (const auto& e : entries_) {
            if (e.label == label) return e.cumulative_ms;
        }
        return -1.0;
    }

    // Check if a label exists
    bool has(const std::string& label) const {
        for (const auto& e : entries_) {
            if (e.label == label) return true;
        }
        return false;
    }

    void print(std::ostream& os = std::cerr) const {
        os << "\n=== Timing Summary ===\n";
        os << std::left << std::setw(30) << "Label" 
           << std::right << std::setw(12) << "Elapsed(ms)" 
           << std::setw(14) << "Cumulative(ms)\n";
        os << std::string(56, '-') << "\n";
        for (const auto& e : entries_) {
            os << std::left << std::setw(30) << e.label
               << std::right << std::setw(12) << std::fixed << std::setprecision(3) << e.elapsed_ms
               << std::setw(14) << e.cumulative_ms << "\n";
        }
    }

    void save_csv(const std::string& filename) const {
        std::ofstream f(filename);
        f << "label,elapsed_ms,cumulative_ms\n";
        for (const auto& e : entries_) {
            f << e.label << "," << e.elapsed_ms << "," << e.cumulative_ms << "\n";
        }
    }

    void save_csv(const std::string& filename, const std::string& dataset, double budget, int num_tiles) const {
        bool file_exists = std::ifstream(filename).good();
        std::ofstream f(filename, std::ios::app);
        
        if (!file_exists) {
            f << "dataset,budget,num_tiles,label,elapsed_ms,cumulative_ms\n";
        }
        
        for (const auto& e : entries_) {
            f << dataset << "," << budget << "," << num_tiles << "," << e.label << "," << e.elapsed_ms << "," << e.cumulative_ms << "\n";
        }
    }

private:
    Clock::time_point start_time_;
    Clock::time_point last_time_;
    std::vector<Entry> entries_;
};

inline TimingLog& get_timer() {
    static TimingLog instance;
    return instance;
}



class MemoryLog {
public:
    struct Entry {
        std::string label;
        long rss_peak_kb;       // Process-lifetime peak RSS (VmHWM)
        long vm_size_kb;        // Current virtual memory
        long rss_kb;            // Current physical memory
        long scope_peak_kb;     // Scoped peak physical memory (0 if using mark())
    };

    static MemoryLog& instance() {
        static MemoryLog inst;
        return inst;
    }

    void mark(const std::string& label) {
        Entry e;
        e.label = label;
        read_memory(e);
        e.scope_peak_kb = 0;
        entries_.push_back(e);
    }

    void sample_start(const std::string& label) {
        // Stop any existing tracker with same label
        if (active_trackers_.count(label)) {
            sample_stop(label);  // BUG FIX: was stop(), should be sample_stop()
        }

        auto tracker = std::make_unique<Tracker>();
        tracker->label = label;
        tracker->rss_before = read_rss();
        tracker->rss_peak = tracker->rss_before;
        tracker->done = false;

        // Start sampling thread
        tracker->sampler = std::thread([t = tracker.get(), this]() {  // BUG FIX: capture 'this' for sample_interval_ms_
            while (!t->done.load()) {
                long current = read_rss();
                long expected = t->rss_peak.load();
                while (current > expected && 
                       !t->rss_peak.compare_exchange_weak(expected, current));
                std::this_thread::sleep_for(std::chrono::milliseconds(sample_interval_ms_));
            }
        });

        active_trackers_[label] = std::move(tracker);
    }

    void sample_stop(const std::string& label) {
        auto it = active_trackers_.find(label);
        if (it == active_trackers_.end()) {
            std::cerr << "Warning: No active tracker for '" << label << "'\n";
            return;
        }

        auto& tracker = it->second;
        tracker->done.store(true);
        if (tracker->sampler.joinable()) {
            tracker->sampler.join();
        }

        long rss_peak = tracker->rss_peak.load();
        long rss_after = read_rss();  // BUG FIX: rss_after was undefined
        
        // Update peak with final reading
        if (rss_after > rss_peak) rss_peak = rss_after;
        
        Entry e;
        e.label = label;
        read_memory(e);
        e.scope_peak_kb = rss_peak;
        entries_.push_back(e);

        active_trackers_.erase(it);
    }

    void set_sample_interval(int ms) { sample_interval_ms_ = ms; }

    void print(std::ostream& os = std::cerr) const {
        os << "\n=== Memory Summary ===\n";
        os << std::left << std::setw(25) << "Label"
           << std::right << std::setw(14) << "RssPeak(MB)"
           << std::setw(12) << "Vm(MB)"
           << std::setw(12) << "RSS(MB)"
           << std::setw(18) << "ScopedPeak(MB)" << "\n";
        os << std::string(81, '-') << "\n";
        
        for (const auto& e : entries_) {
            os << std::left << std::setw(25) << e.label
               << std::right << std::setw(14) << std::fixed << std::setprecision(2) 
               << e.rss_peak_kb / 1024.0
               << std::setw(12) << e.vm_size_kb / 1024.0
               << std::setw(12) << e.rss_kb / 1024.0
               << std::setw(18) << (e.scope_peak_kb > 0 ? std::to_string(e.scope_peak_kb / 1024.0) : "-") 
               << "\n";  // BUG FIX: was scope_rss_kb, should be scope_peak_kb
        }
    }

    void save_csv(const std::string& filename) const {
        bool file_exists = std::ifstream(filename).good();
        std::ofstream f(filename, std::ios::app);
        
        if (!file_exists) {
            f << "label,rss_peak_kb,vm_kb,rss_kb,scope_peak_kb\n";
        }
        
        for (const auto& e : entries_) {
            f << e.label << ","
              << e.rss_peak_kb << ","
              << e.vm_size_kb << ","
              << e.rss_kb << ","
              << e.scope_peak_kb << "\n";  // BUG FIX: was scope_rss_kb
        }
    }

    void save_csv(const std::string& filename, const std::string& dataset, double budget, int num_tiles) const {
        bool file_exists = std::ifstream(filename).good();
        std::ofstream f(filename, std::ios::app);
        
        if (!file_exists) {
            f << "dataset,budget,num_tiles,label,rss_peak_kb,vm_kb,rss_kb,scope_peak_kb\n";
        }
        
        for (const auto& e : entries_) {
            f << dataset << ","
              << budget << ","
              << num_tiles << ","
              << e.label << ","
              << e.rss_peak_kb << ","
              << e.vm_size_kb << ","
              << e.rss_kb << ","
              << e.scope_peak_kb << "\n";  
        }
    }

    const Entry* get_entry(const std::string& label) const {
        for (const auto& e : entries_) {
            if (e.label == label) return &e;
        }
        return nullptr;
    }

    // Get memory delta between two labels
    long get_rss_delta(const std::string& from, const std::string& to) const {
        auto* e1 = get_entry(from);
        auto* e2 = get_entry(to);
        if (e1 && e2) return e2->rss_kb - e1->rss_kb;  // BUG FIX: was vm_rss_kb
        return 0;
    }

    void clear() { 
        // Stop all active trackers first
        for (auto& [label, tracker] : active_trackers_) {
            tracker->done.store(true);
            if (tracker->sampler.joinable()) {
                tracker->sampler.join();
            }
        }
        active_trackers_.clear();
        entries_.clear(); 
    }

private:
    struct Tracker {
        std::string label;
        long rss_before;
        std::atomic<long> rss_peak;
        std::atomic<bool> done;
        std::thread sampler;
    };

    std::vector<Entry> entries_;
    std::unordered_map<std::string, std::unique_ptr<Tracker>> active_trackers_;
    int sample_interval_ms_ = 10;  

    void read_memory(Entry& e) {
        e.rss_peak_kb = e.vm_size_kb = e.rss_kb = 0; 
        
        std::ifstream status("/proc/self/status");
        if (!status.is_open()) {
            std::cerr << "Warning: Cannot open /proc/self/status\n";
            return;
        }

        std::string line;
        while (std::getline(status, line)) {
            std::istringstream iss(line);
            std::string key;
            long value;
            
            if (!(iss >> key >> value)) continue;
            
            if (key == "VmHWM:")       e.rss_peak_kb = value;
            else if (key == "VmSize:") e.vm_size_kb = value;
            else if (key == "VmRSS:")  e.rss_kb = value;
        }
    }

    static long read_rss() {  
        std::ifstream status("/proc/self/status");
        std::string line;
        while (std::getline(status, line)) {
            if (line.compare(0, 6, "VmRSS:") == 0) {
                std::istringstream iss(line.substr(6));
                long value;
                iss >> value;
                return value;
            }
        }
        return 0;
    }
};

// Convenience function
inline MemoryLog& get_memory_log() {
    return MemoryLog::instance();
}

// //print path to stdout
// double print_path(const std::vector<std::vector<double>>& costs, const std::vector<double>& probability, 
//                     std::vector<int> rank, std::vector<int> best_path, double dwell_time);

// //save result to a csv file
// void save_result(const std::string& filename, 
//                   const std::string& method, 
//                   const std::string& data, 
//                   double budget, double slew_rate, 
//                   const std::vector<std::vector<double>>& costs, 
//                   const std::vector<double>& probability, 
//                   std::vector<int> rank,
//                   std::vector<int> best_path, 
//                   double elapsed_time, double padding);

// // save the result of multi deadline to a csv file
// void save_result2(const std::string& filename, 
//                   const std::string& method, 
//                   const std::vector<double>& budgets,
//                   const std::vector<double>& dwell_time, 
//                   const std::vector<std::vector<double>>& costs, 
//                   const std::vector<double>& probability, 
//                   std::vector<int> rank,
//                   std::vector<int> best_path, 
//                   double elapsed_time);
                  
