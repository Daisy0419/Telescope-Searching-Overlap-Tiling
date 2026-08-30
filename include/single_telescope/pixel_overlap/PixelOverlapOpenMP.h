#pragma once

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// Return a practical number of workers. A requested value of zero means use
// the OpenMP runtime default. Without OpenMP, this always returns one.
inline int pixelOverlapWorkerCount(int requested_threads, int work_count) {
    int available_threads = 1;
#ifdef _OPENMP
    available_threads = std::max(1, omp_get_max_threads());
#endif

    int worker_count = requested_threads > 0
        ? std::min(requested_threads, available_threads)
        : available_threads;
    if (work_count > 0) {
        worker_count = std::min(worker_count, work_count);
    }
    return std::max(1, worker_count);
}

// Identify the calling worker. This lets the solver select the private
// CoverageOracle owned by that worker.
inline int pixelOverlapThreadNumber() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}