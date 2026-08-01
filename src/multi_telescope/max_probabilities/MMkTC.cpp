#include "GCP.h"
#include "LocalSearch.h"
#include "PathOperations.h"
#include "LocalSearch.h"

#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <numeric>
#include <iostream>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <cmath>
#include <future>


const double INF = std::numeric_limits<double>::infinity();

struct EDGE {
    int u, v;
    double weight;
    EDGE() : u(-1), v(-1), weight(INF) {}
    EDGE(int u, int v, double w) : u(u), v(v), weight(w) {}
    bool operator<(const EDGE& other) const { return weight < other.weight; }
};

struct Tree {
    std::vector<EDGE> edges;
    double weight;
    int root;
    std::unordered_set<int> nodes;

    Tree() : weight(0), root(-1) {}
};

// Hopcroft–Karp for maximum cardinality bipartite matching
struct HopcroftKarp {
    int L, R;
    std::vector<std::vector<int>> adj; // left u -> list of right v
    std::vector<int> dist, matchL, matchR;
    HopcroftKarp(int L_=0, int R_=0): L(L_), R(R_), adj(L_), dist(L_,0), matchL(L_,-1), matchR(R_,-1) {}
    
    void reset(int L_, int R_) { 
        L=L_; 
        R=R_; 
        adj.assign(L,{}); 
        dist.assign(L,0); 
        matchL.assign(L,-1); 
        matchR.assign(R,-1); 
    }

    void addEDGE(int u,int v) { 
        adj[u].push_back(v); 
    }

    bool bfs() {
        std::queue<int> q; bool reachableFree=false;
        for (int u=0; u<L; ++u){
            if (matchL[u]==-1){ dist[u]=0; q.push(u); }
            else dist[u]=-1;
        }
        while(!q.empty()) {
            int u=q.front(); q.pop();
            for (int v: adj[u]){
                int mu=matchR[v];
                if (mu!=-1 && dist[mu]==-1){ dist[mu]=dist[u]+1; q.push(mu); }
                if (mu==-1) reachableFree=true;
            }
        }
        return reachableFree;
    }

    bool dfs(int u) {
        for (int v: adj[u]){
            int mu=matchR[v];
            if (mu==-1 || (dist[mu]==dist[u]+1 && dfs(mu))){
                matchL[u]=v; matchR[v]=u; return true;
            }
        }
        dist[u]=-1; return false;
    }

    int maxMatching() {
        int m=0;
        while(bfs()){
            for (int u=0; u<L; ++u) if (matchL[u]==-1) if (dfs(u)) ++m;
        }
        return m;
    }
};


template <class SetA, class SetB>
static EDGE cheapest_bridge(const std::vector<std::vector<double>>& costs,
                            const SetA& A, const SetB& B) {
    EDGE best;
    for (int u : A) {
        for (int v : B) {
            double w = costs[u][v];
            if (w < best.weight) best = EDGE(u, v, w);
        }
    }
    return best;
}


class RootedTreeCover {
private:
    int n;
    std::vector<std::vector<double>> costs;
    std::vector<int> roots;

    // Union-Find for Kruskal
    std::vector<int> parentUF, rankUF;
    int findUF(int x){ return parentUF[x]==x? x: parentUF[x]=findUF(parentUF[x]); }
    void uniteUF(int a,int b){
        a=findUF(a); b=findUF(b);
        if (a==b) return;
        if (rankUF[a]<rankUF[b]) parentUF[a]=b;
        else if (rankUF[a]>rankUF[b]) parentUF[b]=a;
        else { parentUF[b]=a; rankUF[a]++; }
    }

    // Build adjacency list of edges <= B
    std::vector<std::vector<std::pair<int,double>>> buildPrunedAdj(double B) const {
        std::vector<std::vector<std::pair<int,double>>> adj(n);
        for (int i=0;i<n;++i){
            for (int j=i+1;j<n;++j){
                double w = costs[i][j];
                if (w<INF && w<=B){
                    adj[i].push_back({j,w});
                    adj[j].push_back({i,w});
                }
            }
        }
        return adj;
    }

    // Dijkstra on the pruned graph
    std::pair<std::vector<double>, std::vector<int>> dijkstra(const std::vector<std::vector<std::pair<int,double>>>& adj, int s) {
        std::vector<double> d(n, INF); std::vector<int> par(n,-1);
        using P = std::pair<double,int>;
        std::priority_queue<P, std::vector<P>, std::greater<P>> pq;
        d[s]=0.0; pq.push({0.0,s});
        while(!pq.empty()){
            auto [du,u]=pq.top(); pq.pop();
            if (du!=d[u]) continue;
            for (auto [v,w] : adj[u]){
                if (d[v] > du + w){
                    d[v] = du + w; par[v]=u; pq.push({d[v], v});
                }
            }
        }
        return {d, par};
    }

    // Compute MST on: (all edges <=B) U (0-edges between every root pair)
    std::vector<EDGE> mst_with_root_contraction(double B){
        std::vector<EDGE> edges;
        edges.reserve(n*n/2 + (int)roots.size()*((int)roots.size()-1)/2);

        // thresholded edges
        for (int i=0;i<n;++i){
            for (int j=i+1;j<n;++j){
                double w = costs[i][j];
                if (w<INF && w<=B) edges.emplace_back(i,j,w);
            }
        }
        // 0-weight edges among roots
        for (int i=0;i<(int)roots.size();++i)
            for (int j=i+1;j<(int)roots.size();++j)
                edges.emplace_back(roots[i], roots[j], 0.0);

        // Kruskal
        parentUF.resize(n); rankUF.assign(n,0);
        std::iota(parentUF.begin(), parentUF.end(), 0);
        std::sort(edges.begin(), edges.end());

        std::vector<EDGE> mst;
        mst.reserve(n-1);
        for (const auto& e : edges){
            int a=findUF(e.u), b=findUF(e.v);
            if (a!=b){ uniteUF(a,b); mst.push_back(e); if ((int)mst.size()==n-1) break; }
        }
        return mst;
    }

    // Build adjacency from edges
    std::vector<std::vector<std::pair<int,double>>>
    buildAdjFromEDGEs(int n, const std::vector<EDGE>& edges){
        std::vector<std::vector<std::pair<int,double>>> adj(n);
        for (auto &e : edges){
            adj[e.u].push_back({e.v, e.weight});
            adj[e.v].push_back({e.u, e.weight});
        }
        return adj;
    }

    std::vector<Tree> split_root_components(const std::vector<EDGE>& mst) {
        std::vector<char> isRoot(n,0); for (int r: roots) isRoot[r]=1;

        std::vector<std::vector<std::pair<int,double>>> adj(n);
        for (auto &e : mst){
            if (e.weight==0.0 && isRoot[e.u] && isRoot[e.v]) continue; // remove root–root 0-edges
            adj[e.u].push_back({e.v, e.weight});
            adj[e.v].push_back({e.u, e.weight});
        }

        std::vector<char> vis(n,0);
        std::vector<Tree> root_comps; 
        root_comps.reserve(roots.size());

        for (int r : roots){
            if (vis[r]) continue;
            Tree rc; rc.root = r;
            std::vector<int> nodes;
            std::queue<int> q; q.push(r); vis[r]=1;
            while(!q.empty()){
                int u=q.front(); q.pop();
                nodes.push_back(u);
                for (auto [v,w] : adj[u]){
                    if (!vis[v]){ vis[v]=1; q.push(v); }
                }
            }
            // collect edges inside this component
            std::unordered_set<int> S(nodes.begin(), nodes.end());
            for (auto &e : mst){
                if (e.weight==0.0 && isRoot[e.u] && isRoot[e.v]) continue;
                if (S.count(e.u) && S.count(e.v)) rc.edges.push_back(e);
            }
            rc.nodes = std::unordered_set<int>(nodes.begin(), nodes.end()); // <-- FIX
            root_comps.push_back(std::move(rc));
        }
        return root_comps;
    }


    struct Bundle { double w = 0.0; std::vector<EDGE> es; };
    std::pair<std::vector<Tree>, Tree>
    rooted_decompose(const std::vector<EDGE>& treeEDGEs, int root, double B) {
        // Build adjacency for this tree
        auto adj = buildAdjFromEDGEs(n, treeEDGEs);
        std::vector<Tree> pieces;

        auto emit_piece = [&](const std::vector<EDGE>& es) {
            Tree P; P.edges = es; P.weight = 0.0;
            for (const auto& e : P.edges) {
                P.nodes.insert(e.u); P.nodes.insert(e.v);
                P.weight += e.weight;
            }
            pieces.push_back(std::move(P));
        };

        // dfs(u,p) returns the leftover (<B) bundle from the subtree of u (toward parent p).
        std::function<Bundle(int,int)> dfs = [&](int u, int p) -> Bundle {
            // Gather child bundles (child leftover + edge (u,v))
            std::vector<Bundle> bundles;
            for (auto [v, w] : adj[u]) {
                if (v == p) continue;
                Bundle ch = dfs(v, u);
                ch.w += w;
                ch.es.push_back(EDGE(u, v, w));
                bundles.push_back(std::move(ch));
            }

            // Emit any large child bundles (≥ B) as pieces; keep "small" ones (< B)
            std::vector<Bundle> small;
            small.reserve(bundles.size());
            for (auto& b : bundles) {
                if (b.w >= B) {
                    emit_piece(b.es);
                } else {
                    small.push_back(std::move(b)); // b.w < B
                }
            }

            // Greedy packing of smalls so that each emitted pack is in [B,2B)
            std::sort(small.begin(), small.end(),
                    [](const Bundle& a, const Bundle& b){ return a.w > b.w; });

            Bundle carry; // < B
            for (auto& b : small) {
                if (carry.w + b.w < B) {
                    carry.w += b.w;
                    carry.es.insert(carry.es.end(), b.es.begin(), b.es.end());
                } else {
                    std::vector<EDGE> pack = carry.es;
                    pack.insert(pack.end(), b.es.begin(), b.es.end());
                    emit_piece(pack);
                    carry = Bundle{}; // reset leftover
                }
            }
            return carry; 
        };

        Bundle rootLeft = dfs(root, -1);

        Tree L;
        L.root = root;
        L.edges = std::move(rootLeft.es);
        L.weight = 0.0;
        for (const auto& e : L.edges) {
            L.nodes.insert(e.u); L.nodes.insert(e.v);
            L.weight += e.weight;
        }
        L.nodes.insert(root); //leftover always contains the root
        return {pieces, L};
    }

    std::pair<std::vector<Tree>, Tree>
    rooted_decompose2(const std::vector<EDGE>& treeEDGEs, int root, double B) {
        auto adj = buildAdjFromEDGEs(n, treeEDGEs);

        std::vector<Tree> pieces;
        auto emit_piece = [&](const std::vector<EDGE>& es) {
            Tree P; P.root = root; P.edges = es; P.weight = 0.0;
            for (const auto& e : P.edges) { P.nodes.insert(e.u); P.nodes.insert(e.v); P.weight += e.weight; }
            pieces.push_back(std::move(P));
        };

        std::function<Bundle(int,int)> dfs = [&](int u, int p) -> Bundle {
            // Collect child bundles; add (u,v) edge so each child becomes (0, 2B)
            std::vector<Bundle> items;
            for (auto [v, w] : adj[u]) {
                if (v == p) continue;
                Bundle b = dfs(v, u);   // b.w < B guaranteed from recursion
                b.w += w;               // now 0 < b.w < 2B
                b.es.push_back(EDGE(u, v, w));
                items.push_back(std::move(b));
            }

            // Sort heavy-to-light to fill bigs first
            std::sort(items.begin(), items.end(),
                    [](const Bundle& a, const Bundle& b){ return a.w > b.w; });

            std::vector<char> used(items.size(), 0);

            // Pass 1: emit each big (>=B), topping up with smalls while staying < 2B
            for (size_t i = 0; i < items.size(); ++i) {
                if (used[i] || items[i].w < B) continue;
                used[i] = 1;
                double accW = items[i].w;
                std::vector<EDGE> accE = items[i].es;
                for (size_t j = i + 1; j < items.size(); ++j) {
                    if (used[j] || items[j].w >= B) continue;     // only use smalls to top-up
                    if (accW + items[j].w < 2.0 * B) {
                        accW += items[j].w;
                        accE.insert(accE.end(), items[j].es.begin(), items[j].es.end());
                        used[j] = 1;
                    }
                }
                // accW ∈ [B, 2B)
                emit_piece(accE);
            }

            // Pass 2: pack remaining smalls into [B,2B), keep exactly one leftover < B
            Bundle carry; // the only leftover returned upward (kept < B)
            for (size_t i = 0; i < items.size(); ++i) {
                if (used[i] || items[i].w >= B) continue; // skip used or bigs
                // Try to absorb into carry
                if (carry.w + items[i].w < B) {
                    carry.w += items[i].w;
                    carry.es.insert(carry.es.end(), items[i].es.begin(), items[i].es.end());
                    used[i] = 1;
                    continue;
                }

                // Form a pack: carry + this small crosses B ⇒ emit piece in [B,2B)
                std::vector<EDGE> accE = carry.es; double accW = carry.w;
                accE.insert(accE.end(), items[i].es.begin(), items[i].es.end());
                accW += items[i].w;
                used[i] = 1;
                for (size_t j = i + 1; j < items.size(); ++j) {
                    if (used[j] || items[j].w >= B) continue;
                    if (accW + items[j].w < 2.0 * B) {
                        accW += items[j].w;
                        accE.insert(accE.end(), items[j].es.begin(), items[j].es.end());
                        used[j] = 1;
                    }
                }
                emit_piece(accE);
                carry = Bundle{}; // reset leftover
            }

            // Guaranteed carry.w < B here
            return carry;
        };

        Bundle rootLeft = dfs(root, -1);

        // Leftover: ensure it includes the root even if no edges
        Tree L; L.root = root; L.edges = std::move(rootLeft.es); L.weight = 0.0;
        for (const auto& e : L.edges) { L.nodes.insert(e.u); L.nodes.insert(e.v); L.weight += e.weight; }
        L.nodes.insert(root);
        return {pieces, L};
    }


public:
    RootedTreeCover(const std::vector<std::vector<double>>& costs_,
                    const std::vector<int>& rootIndices)
        : costs(costs_), roots(rootIndices), n((int)costs_.size()) {}

    std::pair<bool, std::vector<Tree>> solve(double B) {
        // ---- 1) Reachability in the pruned graph (edges <= B) ----
        auto prunedAdj = buildPrunedAdj(B);
        {
            std::vector<char> seen(n,0);
            std::queue<int> q;
            for (int r : roots) if (!seen[r]) { seen[r]=1; q.push(r); }
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (auto [v,w] : prunedAdj[u]) if (!seen[v]) { seen[v]=1; q.push(v); }
            }
            for (int v=0; v<n; ++v) if (!seen[v]) return {false, {}};
        }

        // ---- 2) MST on (edges <= B) with roots contracted via 0-edges ----
        auto mst = mst_with_root_contraction(B);
        if ((int)mst.size() < n-1) return {false, {}}; // disconnected under <=B

        // ---- 3) Remove root-root 0-edges -> forest of rooted components ----
        auto comps = split_root_components(mst);
        if ((int)comps.size() > (int)roots.size()) return {false, {}};

        // Build root -> component index map by membership
        std::vector<int> root_to_comp(roots.size(), -1);
        for (int ci=0; ci<(int)comps.size(); ++ci) {
            for (int rj=0; rj<(int)roots.size(); ++rj) {
                if (comps[ci].nodes.count(roots[rj])) root_to_comp[rj] = ci;
            }
        }
        for (int rj=0; rj<(int)roots.size(); ++rj)
            if (root_to_comp[rj] == -1) return {false, {}};

        // ---- 4) Decompose each component into [B,2B) pieces + leftover < B ----
        std::vector<Tree> pieces;      pieces.reserve(n);
        std::vector<Tree> leftovers;   leftovers.reserve(comps.size());

        for (int ci=0; ci<(int)comps.size(); ++ci) {
            auto [pieceTrees, leftoverTree] = rooted_decompose2(comps[ci].edges, comps[ci].root, B);
            for (auto &PT : pieceTrees) pieces.push_back(std::move(PT)); // pieces are Trees
            leftovers.push_back(std::move(leftoverTree));                // leftover is a Tree (< B)
        }

        // Necessary condition
        if ((int)pieces.size() > (int)roots.size()) return {false, {}};

        // ---- 6) Build bipartite graph using CHEAPEST EDGE to leftover and a 4B budget test ----
        // For reconstruction, stash the best bridge for every (piece i, root j).
        std::vector<std::vector<EDGE>> bestBridge(pieces.size(),
                                                std::vector<EDGE>(roots.size(), EDGE()));

        HopcroftKarp hk((int)pieces.size(), (int)roots.size());
        for (int i = 0; i < (int)pieces.size(); ++i) {
            const Tree& P = pieces[i];
            for (int j = 0; j < (int)roots.size(); ++j) {
                int cj = root_to_comp[j];                   // component index for root j
                const Tree& L = leftovers[cj];              // leftover tree for that component
                EDGE br = cheapest_bridge(costs, P.nodes, L.nodes);
                // Accept match if we can connect via a single edge and total per-root cost ≤ 4B
                if (br.weight < INF && (P.weight + L.weight + br.weight) <= 4.0 * B + 1e-12) {
                    hk.addEDGE(i, j);
                    bestBridge[i][j] = br;                 // remember how to connect
                }
            }
        }
        if (hk.maxMatching() != (int)pieces.size()) return {false, {}};

        // Map root j -> assigned piece index (or -1)
        std::vector<int> rootToPiece(roots.size(), -1);
        for (int i = 0; i < (int)pieces.size(); ++i) {
            int rj = hk.matchL[i];
            if (rj >= 0) rootToPiece[rj] = i;
        }

        // ---- 7) Assemble final trees: leftover + (optional) piece + single bridge edge ----
        std::vector<Tree> result; result.reserve(roots.size());

        for (int j = 0; j < (int)roots.size(); ++j) {
            int r = roots[j];
            int cj = root_to_comp[j];
            const Tree& L = leftovers[cj];

            Tree T; T.root = r;

            // Add leftover (always)
            T.edges.insert(T.edges.end(), L.edges.begin(), L.edges.end());
            for (const auto& e : L.edges) {
                T.weight += e.weight;
                T.nodes.insert(e.u); T.nodes.insert(e.v);
            }

            // If this root received a piece, add the piece and the stored bridge edge
            if (rootToPiece[j] != -1) {
                int i = rootToPiece[j];
                const Tree& P = pieces[i];
                const EDGE& br = bestBridge[i][j];

                // (Defensive) Check again our 4B budget
                if (br.weight >= INF || (T.weight + P.weight + br.weight) > 4.0 * B + 1e-12)
                    return {false, {}}; // should not happen if we built the graph correctly

                // add piece edges
                T.edges.insert(T.edges.end(), P.edges.begin(), P.edges.end());
                for (const auto& e : P.edges) {
                    T.weight += e.weight;
                    T.nodes.insert(e.u); T.nodes.insert(e.v);
                }

                // add the single bridge edge
                T.edges.push_back(br);
                T.weight += br.weight;
                T.nodes.insert(br.u); T.nodes.insert(br.v);
            }

            if (T.edges.empty()) T.nodes.insert(r); // singleton root (degenerate)
            // (Optional) assert(T.weight <= 4.0 * B + 1e-9);

            result.push_back(std::move(T));
        }

        return {true, result};
    }


    // Binary search for minimal feasible B
    std::vector<Tree> findOptimalCover(double budget) {
        // Get an initial high bound: max finite edge
        // double hi = 0.0;
        // for (int i=0;i<n;++i) for (int j=0;j<n;++j)
        //     if (costs[i][j]<INF) hi = std::max(hi, costs[i][j]);
        // if (hi <= 0.0) hi = 1.0;
        double hi = budget;

        // Exponential grow if needed
        auto probe = solve(hi);
        int guard=0;
        while(!probe.first && guard<25){ hi *= 2.0; probe = solve(hi); ++guard; }
        if (!probe.first) return {}; // give up

        double lo = 0.0;
        std::vector<Tree> best = probe.second; // feasible at hi
        for (int it=0; it<100; ++it){
            double mid = (lo+hi)/2.0;
            auto fr = solve(mid);
            if (fr.first){ best = std::move(fr.second); hi = mid; }
            else lo = mid;
            if (hi - lo <= 1e-6 * std::max(1.0, hi)) break;
        }
        return best;
    }
};

// Build an s–t path from one local tree.
std::vector<int> buildPathFromLocalTree(const Tree& tree,
                                        const std::vector<std::vector<double>>& costs,
                                        const std::vector<int>& indices, 
                                        const std::vector<int>& local_roots, 
                                        int t_global) {
    // using Graph = lemon::ListGraph;
    // using Node  = Graph::Node;
    // using Edge  = Graph::Edge;

    std::vector<int> cf_indices;
    cf_indices.reserve(tree.nodes.size() + 2);

    if (std::find(local_roots.begin(), local_roots.end(), tree.root) == local_roots.end())
        throw std::runtime_error("Wrong root in a tree!");

    const int root_global = indices[tree.root];
    cf_indices.push_back(root_global);
    cf_indices.push_back(t_global);
    // if (t_global != root_global) cf_indices.push_back(t_global);

    for (int u_local : tree.nodes) {
        int g = indices[u_local];
        if (g == root_global) continue;
        if (cf_indices.size() >= 2 && g == cf_indices[1]) continue;
        cf_indices.push_back(g);
    }

    // Trivial
    if (cf_indices.size() == 1) return cf_indices;
    if (cf_indices.size() == 2) return cf_indices;  // [root, t]

    // Build complete graph on positions 0..M-1
    const int M = (int)cf_indices.size();
    Graph g;
    std::vector<Node> node(M);
    for (int i = 0; i < M; ++i) node[i] = g.addNode();

    Graph::EdgeMap<double> w(g);
    for (int i = 0; i < M; ++i)
        for (int j = i + 1; j < M; ++j) {
            Edge e = g.addEdge(node[i], node[j]);
            w[e] = costs[cf_indices[i]][cf_indices[j]];
            // std::cout << w[e] << " ";
        }
    // std::cout << std::endl;   


    // MST over the compact set
    std::vector<Edge> mst_edges;
    buildMST(g, w, mst_edges);

    int s_idx = 0;
    int t_idx = 1;

    // Christofides on compact set
    std::vector<int> path_pos = christofidesPathTwoFixed(g, costs, cf_indices, M, mst_edges, s_idx, t_idx);

    std::vector<int> path_global; 
    path_global.reserve(path_pos.size());
    for (int p : path_pos) path_global.push_back(cf_indices[p]);
    return path_global;
}


std::vector<std::vector<double>> UpdateGraph(const std::vector<std::vector<double>>& costs, 
                                             const std::vector<int>& indices,
                                             const std::vector<int>& local_roots,
                                             int N) {
    std::vector<std::vector<double>> sub(N, std::vector<double>(N, INF));
    std::unordered_set<int> rootSet(local_roots.begin(), local_roots.end());

    for (int i = 0; i < N; ++i) {
        sub[i][i] = 0.0;
        for (int j = i + 1; j < N; ++j) {
            double w = costs[indices[i]][indices[j]];
            sub[i][j] = sub[j][i] = w;
            // Optional: 0-weight among any pair of roots to simulate contraction
            if (rootSet.count(i) && rootSet.count(j)) {
                sub[i][j] = sub[j][i] = 0.0;
            }
        }
    }
    return sub;
}


double binarySearchBestCover(const std::vector<std::vector<double>>& costs, 
                             const std::vector<double>& prizes,
                             const std::vector<int>& indices,      // global ids in fixed order
                             const std::vector<int>& local_roots,   // local ids (0..R-1) of roots
                             int t,                                 // global target id
                             std::vector<std::vector<int>>& best_paths, 
                             double budget) {

    int N_min = (int)local_roots.size();  
    int N_max = (int)indices.size() + 1; 
    double best_prize = -1.0;
    best_paths.clear();

    // std::vector<int> global_roots;
    // for(int i : local_roots) {
    //     global_roots.push_back(indices[i]);
    // }

    if (N_min > (int)indices.size()) return best_prize;

    while (N_max - N_min > 1) {
        int N = N_min + (N_max - N_min) / 2;
        // std::cout << "N: " << N << std::endl;

        // Build N-node subgraph (0-weight among roots)
        std::vector<std::vector<double>> sub_costs = UpdateGraph(costs, indices, local_roots, N);

        // Solve rooted cover on the N-node subgraph (LOCAL ids)
        RootedTreeCover rtc(sub_costs, local_roots);
        auto trees = rtc.findOptimalCover(budget);

        bool feasible = true;
        std::vector<std::vector<int>> paths;
        if (trees.empty()) {
            feasible = false;
        } else {
            // Build s–t paths (GLOBAL ids). Use the same N-node mapping via `indices`.
            paths.reserve(trees.size());
            try {
                for (const auto& tree : trees) {
                    // buildPathFromLocalTree may throw (e.g., disconnected compact graph).
                    auto path = buildPathFromLocalTree(tree, costs, indices, local_roots, t);
                    if (path.empty()) { feasible = false; break; }
                    two_opt_st(path, costs);
                    paths.push_back(std::move(path));
                }
                // repairMultiPaths(costs, prizes, paths, budget, global_roots,t);
                // insertHighValueCheapestMulti(costs, prizes, paths, budget, global_roots, t);
            } catch (const std::exception&) {
                feasible = false;
            }

            if (feasible) {
                double max_weight = multiPathsMaxCost(costs, paths);
                feasible = (max_weight <= budget);
            }
        }

        if (!feasible) {
            // Shrink upper bound when infeasible at N
            N_max = N;
            continue;
        }

        // repairMultiPaths(costs,prizes,paths,budget)
        // Feasible at N: record best prize and grow lower bound
        double prize = multiPathsSumPrize(prizes, paths);
        if (prize > best_prize) {
            best_prize = prize;
            best_paths = std::move(paths);
        }
        N_min = N;
    }

    return best_prize; 
}


std::pair<std::vector<std::vector<int>>, double> bestPrizeMMkTC(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m){

    std::vector<int> indices =  start_indices;
    std::vector<int> local_roots(start_indices.size());
    std::iota (local_roots.begin(), local_roots.end(), 0);

    std::vector<int> indices_tail(prizes.size());
    std::iota (indices_tail.begin(), indices_tail.end(), 0);
    std::sort(indices_tail.begin(), indices_tail.end(),
            [&](const int& a, const int& b){ return prizes[a] > prizes[b]; });


    for(int i : indices_tail) {
        if(std::count(start_indices.begin(), start_indices.end(), i) > 0)
            continue;
        indices.push_back(i);
    }

    std::vector<std::vector<int>> best_paths;
    double best_prize = binarySearchBestCover(costs, prizes, indices, local_roots, t, best_paths, budget*8);
    for(auto& path : best_paths) {
        bool is_improve = true;
        while(is_improve) {
            is_improve = false;
            is_improve |= two_opt_st(path, costs);
        }
    }
    repairMultiPaths(costs, prizes, best_paths, budget, start_indices,t);
    best_prize = multiPathsSumPrize( prizes, best_paths);
    return {best_paths, best_prize};

}


void selectNodesPrizeRatio(const std::vector<std::vector<double>>& costs,
                            const std::vector<double>& prizes,
                            std::vector<int>& indices,
                            const std::vector<int>& start_indices,
                            int t) {
    
    const int n = costs.size();
    std::vector<bool> selected(n, false);
    std::vector<double> minDist(n, std::numeric_limits<double>::max());
    std::priority_queue<std::pair<double, int>> heap; // Max-heap for prize/distance ratio
    const double eps = 1e-12;

    for(int s_idx : start_indices) {
        selected[s_idx] = true;
        indices.push_back(s_idx);
    }

    selected[t] = true;

    // init min dist and heap
    for (int v = 0; v < n; ++v) {
        if (std::count(start_indices.begin(), start_indices.end(), v) > 0 || v == t) 
            continue;

        for(int s_idx : start_indices)
            minDist[v] = std::min(costs[s_idx][v], minDist[v]);

        double score = prizes[v] / (minDist[v] + eps);
        heap.push({score, v});
    }
    

    // Helper to compute score
    auto safe_score = [&](int v) -> double {

        double d = minDist[v];
        if (!std::isfinite(d)) return -std::numeric_limits<double>::infinity();
        double denom = d + eps;
        if (!(denom > 0.0) || !std::isfinite(denom)) denom = eps;
        double s = prizes[v] / denom;
        if (!std::isfinite(s)) {
            s = -std::numeric_limits<double>::infinity();
        }
        return s;
    };

    while (!heap.empty()) {
        auto [best_score, u] = heap.top();
        heap.pop();

        if (selected[u]) continue; //stale entry

        selected[u] = true;
        indices.push_back(u);

        //update minDist for unselected nodes
        for (int v = 0; v < n; ++v) {
            if (!selected[v]) {
                double newDist = costs[u][v];
                if (newDist < minDist[v]) {
                    minDist[v] = newDist;
                }
                double new_score = prizes[v]/ (minDist[v] + eps);
                heap.push({new_score, v});
            }
        } 
    }
}


std::pair<std::vector<std::vector<int>>, double> bestPrizeRatioMMkTC(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m){

    std::vector<int> indices;
    std::vector<int> local_roots(start_indices.size());
    std::iota (local_roots.begin(), local_roots.end(), 0);

    selectNodesPrizeRatio(costs, prizes,indices,start_indices,t);

    std::vector<std::vector<int>> best_paths;
    double best_prize = binarySearchBestCover(costs, prizes, indices, local_roots, t, best_paths, budget*8);
    
    for(auto& path : best_paths) {
        bool is_improve = true;
        while(is_improve) {
            // is_improve = false;
            is_improve = two_opt_st(path, costs);
        }
    }
    // while(is_improve) {
    //     is_improve = false;
    //     for(auto path : best_paths) 
    //         is_improve |= two_opt_st(path, costs);
    // }
    repairMultiPaths(costs, prizes, best_paths, budget, start_indices,t);
    best_prize = multiPathsSumPrize( prizes, best_paths);
    return {best_paths, best_prize};

}


std::pair<std::vector<std::vector<int>>, double> mmktc_wraper(const std::vector<std::vector<double>>& costs, 
                                                    const std::vector<double>& prizes,
                                                    double budget,
                                                    const std::vector<int>& start_indices,
                                                    int t, int m) {

    auto f1 = std::async(std::launch::async, bestPrizeMMkTC,  costs, prizes, budget, start_indices, t, m);
    auto f2 = std::async(std::launch::async, bestPrizeRatioMMkTC, costs, prizes, budget, start_indices, t, m);
    

    std::pair<std::vector<std::vector<int>>, double> p1 = f1.get();
    std::pair<std::vector<std::vector<int>>, double> p2 = f2.get();

    if (p1.second < p2.second) {
        return p2;
    } else {
        return p1;
    }
}
