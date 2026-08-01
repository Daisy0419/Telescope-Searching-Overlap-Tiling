#include <queue>
#include <algorithm> 

std::pair<std::vector<std::vector<int>>, double> MultiGreedy(
        const std::vector<std::vector<double>>& costs,
        const std::vector<double>& prizes,
        double budget, int nPath,
        const std::vector<int>& starts, int end) {

    const int n = static_cast<int>(prizes.size());

    std::vector<std::pair<double, int>> prize(n);
    for (int i = 0; i < n; ++i) prize[i] = {prizes[i], i};
    std::sort(prize.begin(), prize.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    std::vector<bool> visited(n, false);
    std::vector<std::vector<int>> paths(nPath);
    double totalPrize = 0.0;

    for (int k = 0; k < nPath; ++k) {
        int s = starts[k];
        paths[k].push_back(s);
        if (!visited[s]) {
            visited[s] = true;
            totalPrize += prizes[s];
        }
    }
    if (end >= 0 && end < n) visited[end] = true;

    // Heap key: accumulated time INCLUDING the closing leg to `end`.
    std::priority_queue<std::pair<double, int>, std::vector<std::pair<double, int>>, std::greater<std::pair<double, int>>> heap;
    for (int k = 0; k < nPath; ++k) {
        double init = costs[starts[k]][end];  // start -> end closing leg
        if (init <= budget) {
            heap.push({init, k});
        } else {
            paths[k].push_back(end);
        }
    }
    std::vector<std::size_t> scanFrom(nPath, 0);

    while (!heap.empty()) {
        auto [acc, k] = heap.top();   
        heap.pop();

        int cur = paths[k].back();

        // Try to swap the closing leg for (cur -> next -> end).
        int best = -1;
        double bestAcc = 0.0;
        // for (const auto& [p, i] : prize) {
        //     if (visited[i]) continue;
        //     double newAcc = acc - costs[cur][end] + costs[cur][i] + costs[i][end];
        //     if (newAcc <= budget) {
        //         best = i;
        //         bestAcc = newAcc;
        //         break;
        //     }
        // }
        for (std::size_t q = scanFrom[k]; q < prize.size(); ++q) {
            const int i = prize[q].second;
            if (visited[i]) continue;

            const double newAcc =
                acc - costs[cur][end] + costs[cur][i] + costs[i][end];

            if (newAcc <= budget) {
                best = i;
                bestAcc = newAcc;
                scanFrom[k] = q + 1;   // never revisit earlier prizes for this path
                break;
            }
        }

        if (best == -1) {
            // Nothing affordable: commit the closing leg.
            paths[k].push_back(end);
            continue;
        }

        visited[best] = true;
        totalPrize += prizes[best];
        paths[k].push_back(best);
        heap.push({bestAcc, k});
    }

    return {paths, totalPrize};
}


