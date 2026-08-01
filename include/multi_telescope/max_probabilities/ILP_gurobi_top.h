#pragma once

#include <gurobi_c++.h>
#include <vector>
#include <fstream>
#include <iomanip>  
#include <cmath> 
#include <queue>

struct ProblemDataTeamST {
    int N;                                      // number of nodes
    int K;                                      // number of teams (exactly m)
    int startIndex;                             // common start s
    int endIndex;                               // common end t
    double BudgetPerTeam;                       // per-team budget (scalar)
    std::vector<std::vector<double>> Cost;      // NxN
    std::vector<double> Prize;                  // N, set Prize[s]=Prize[t]=0 if desired

    ProblemDataTeamST(int N_, int K_, int s_, int t_, double B_,
                      const std::vector<std::vector<double>>& C_,
                      const std::vector<double>& P_)
        : N(N_), K(K_), startIndex(s_), endIndex(t_), BudgetPerTeam(B_), Cost(C_), Prize(P_) {}
};


struct ProblemDataTeamMultiS {
    int N;                                      // number of nodes
    int K;                                      // number of teams (exactly m)
    std::vector<int> startIndices;                             // common start s
    int endIndex;                               // common end t
    double BudgetPerTeam;                       // per-team budget (scalar)
    std::vector<std::vector<double>> Cost;      // NxN
    std::vector<double> Prize;                  // N, set Prize[s]=Prize[t]=0 if desired

    ProblemDataTeamMultiS(int N_, int K_, const std::vector<int>& s_, int t_, double B_,
                      const std::vector<std::vector<double>>& C_,
                      const std::vector<double>& P_)
        : N(N_), K(K_), startIndices(s_), endIndex(t_), BudgetPerTeam(B_), Cost(C_), Prize(P_) {}
};

std::pair<std::vector<std::vector<int>>, double>
gurobiSolveTeamST(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  int start, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes={});


std::pair<std::vector<std::vector<int>>, double>
gurobiSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  int start, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes={});


// std::pair<std::vector<std::vector<int>>, double>
std::tuple<std::vector<std::vector<int>>, double, double>
gurobiSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                  const std::vector<double>& Prize,
                  std::vector<int> starts, int end, double BudgetPerTeam,
                  int mTeams,
                  double mipGap, double timeLimit,
                  const std::vector<std::vector<int>>& initialRoutes={});

struct TopConnCallback : public GRBCallback {
    const std::vector<std::vector<GRBVar>>& x;
    const std::vector<GRBVar>* y;   // nullptr for OP
    int N, s, t;
    double eps;
    TopConnCallback(const std::vector<std::vector<GRBVar>>& x_,
                    const std::vector<GRBVar>* y_,
                    int N_, int s_, int t_, double eps_=1e-6)
        : x(x_), y(y_), N(N_), s(s_), t(t_), eps(eps_) {}

    void bfsForward(int src,
                    const std::vector<std::vector<double>>& xv,
                    std::vector<char>& vis) {
        std::queue<int> q; q.push(src); vis[src]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v=0; v<N; ++v) if (u!=v && !vis[v] && xv[u][v] > eps){
                vis[v]=1; q.push(v);
            }
        }
    }
    void bfsReverse(int sink,
                    const std::vector<std::vector<double>>& xv,
                    std::vector<char>& vis) {
        std::queue<int> q; q.push(sink); vis[sink]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v=0; v<N; ++v) if (u!=v && !vis[v] && xv[v][u] > eps){
                vis[v]=1; q.push(v);
            }
        }
    }

    void addCutFromSet(const std::vector<char>& S, double rhs) {
        GRBLinExpr lhs = 0.0;
        for(int i=0;i<N;++i) if (S[i])
            for(int j=0;j<N;++j) if (!S[j] && i!=j)
                lhs += x[i][j];
        addLazy(lhs >= rhs);
    }

    void callback() override {
        if (where != GRB_CB_MIPSOL) return;

        // read solution
        std::vector<std::vector<double>> xv(N, std::vector<double>(N, 0.0));
        for(int i=0;i<N;++i) for(int j=0;j<N;++j) if (i!=j)
            xv[i][j] = getSolution(x[i][j]);

        // (1) s->t cut
        std::vector<char> S(N,0);
        bfsForward(s, xv, S);
        if (!S[t]) {
            addCutFromSet(S, 1.0);
            return; // one cut per callback is fine
        }

        // (2) visited-inside cuts (TOP only)
        if (y) {
            for (int h=0; h<N; ++h) {
                if (h==s || h==t) continue;
                double yh = getSolution((*y)[h]);
                if (yh <= eps) continue;

                std::vector<char> Sh(N,0);
                bfsReverse(h, xv, Sh); // nodes that can reach h
                if (Sh[s] || Sh[t]) continue;

                // compute capacity across delta+(S_h)
                double cap = 0.0;
                for(int i=0;i<N;++i) if (Sh[i])
                    for(int j=0;j<N;++j) if (!Sh[j] && i!=j)
                        cap += xv[i][j];

                if (cap + 1e-8 < yh) { // violated
                    addCutFromSet(Sh, yh);
                    return;
                }
            }
        }
    }
};


// //how to use it
// model.set(GRB_IntParam_LazyConstraints, 1);
// TopConnCallback cb(x, /* TOP: */ &y, /* OP: */ nullptr, N, s, t);
// model.setCallback(&cb);


struct BCCallback : public GRBCallback {
    // Problem data / variables
    const int N, s, t;
    const std::vector<std::vector<GRBVar>>& x;
    const std::vector<GRBVar>& y;

    // Tuning
    const double EPS = 1e-7;
    const int MAX_CUTS_PER_CALL = 16;

    BCCallback(int N_, int s_, int t_,
               const std::vector<std::vector<GRBVar>>& x_,
               const std::vector<GRBVar>& y_)
      : N(N_), s(s_), t(t_), x(x_), y(y_) {}

    // ----- Tiny directed max-flow (Edmonds–Karp) -----
    struct Edge { int to, rev; double cap; };
    struct MF {
        int n;
        std::vector<std::vector<Edge>> g;
        MF(int n_) : n(n_), g(n_) {}
        void addEdge(int u, int v, double c) {
            Edge a{v, (int)g[v].size(), c};
            Edge b{u, (int)g[u].size(), 0.0};
            g[u].push_back(a); g[v].push_back(b);
        }
        double maxflow(int s, int t) {
            double flow = 0.0;
            std::vector<int> parV(n), parE(n);
            while (true) {
                std::fill(parV.begin(), parV.end(), -1);
                std::queue<int> q; q.push(s); parV[s] = s; parE[s] = -1;
                while (!q.empty() && parV[t] == -1) {
                    int u = q.front(); q.pop();
                    for (int i = 0; i < (int)g[u].size(); ++i) {
                        auto &e = g[u][i];
                        if (parV[e.to] == -1 && e.cap > 1e-12) {
                            parV[e.to] = u; parE[e.to] = i; q.push(e.to);
                            if (e.to == t) break;
                        }
                    }
                }
                if (parV[t] == -1) break; // no augmenting path
                // augment
                double aug = std::numeric_limits<double>::infinity();
                for (int v = t; v != s; v = parV[v]) {
                    auto &e = g[parV[v]][parE[v]];
                    aug = std::min(aug, e.cap);
                }
                for (int v = t; v != s; v = parV[v]) {
                    auto &e = g[parV[v]][parE[v]];
                    auto &rev = g[e.to][e.rev];
                    e.cap -= aug; rev.cap += aug;
                }
                flow += aug;
            }
            return flow;
        }
        // Nodes reachable from s in the residual graph after maxflow
        std::vector<char> sourceSideReachable(int s) {
            std::vector<char> vis(n, 0);
            std::queue<int> q; q.push(s); vis[s] = 1;
            while (!q.empty()) {
                int u = q.front(); q.pop();
                for (auto &e : g[u]) if (!vis[e.to] && e.cap > 1e-12) {
                    vis[e.to] = 1; q.push(e.to);
                }
            }
            return vis;
        }
    };

    // Build cut for a subset S (customers only) and an index h in S
    void add_connectivity_cut(const std::vector<char>& inS, int h, bool asLazy) {
        GRBLinExpr lhs = 0.0;
        // sum_{(p,q) in δ+(S)} x_{pq}  −  y_h  >= 0
        for (int p = 0; p < N; ++p) {
            if (!inS[p]) continue;
            if (p == s || p == t) continue;  // S only over customers
            for (int q = 0; q < N; ++q) if (q != p) {
                if (!inS[q]) lhs += x[p][q];
            }
        }
        lhs -= y[h];

        if (asLazy) addLazy(lhs >= 0.0);
        else        addCut (lhs >= 0.0);
    }

    // Separation at a fractional LP node (user cuts)
    int separate_fractional() {
        // Pull current LP relaxation
        std::vector<std::vector<double>> X(N, std::vector<double>(N, 0.0));
        std::vector<double> Y(N, 0.0);
        for (int i = 0; i < N; ++i) {
            Y[i] = getNodeRel(y[i]);
            for (int j = 0; j < N; ++j) if (i != j)
                X[i][j] = getNodeRel(x[i][j]);
        }

        int cuts_added = 0;

        // Build support graph once (capacities = X_ij)
        MF mf(N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (i != j && X[i][j] > EPS)
                    mf.addEdge(i, j, X[i][j]);

        // For each customer i, run min s=i, t=t cut
        for (int src = 0; src < N && cuts_added < MAX_CUTS_PER_CALL; ++src) {
            if (src == s || src == t) continue;

            // We need a fresh residual for each maxflow; rebuild quickly
            MF local = mf;
            (void)local.maxflow(src, t);
            auto inS = local.sourceSideReachable(src);

            // Build S_customers and compute cut capacity and max y in S
            double cutCap = 0.0, maxY = 0.0;
            int h = -1;

            // collect maxY over customers in S
            for (int v = 0; v < N; ++v) {
                if (inS[v] && v != s && v != t) {
                    if (Y[v] > maxY + EPS) { maxY = Y[v]; h = v; }
                }
            }
            if (h == -1) continue; // no customer in S

            // capacity of δ+(S_customers)
            for (int p = 0; p < N; ++p) {
                if (!inS[p]) continue;
                if (p == s || p == t) continue;
                for (int q = 0; q < N; ++q) if (q != p) {
                    if (!inS[q]) cutCap += X[p][q];
                }
            }

            // Violation test: cutCap < max_{v in S} y_v
            if (cutCap + 1e-6 < maxY) {
                add_connectivity_cut(inS, h, /*asLazy=*/false);
                ++cuts_added;
            }
        }
        return cuts_added;
    }

    // Check integer incumbent and add lazy CCs if needed
    int separate_integer() {
        std::vector<std::vector<double>> X(N, std::vector<double>(N, 0.0));
        std::vector<double> Y(N, 0.0);
        for (int i = 0; i < N; ++i) {
            Y[i] = getSolution(y[i]);
            for (int j = 0; j < N; ++j) if (i != j)
                X[i][j] = getSolution(x[i][j]);
        }

        // Build support graph of selected arcs (cap=1 if x>0)
        MF mf(N);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) if (i != j && X[i][j] > 0.5)
                mf.addEdge(i, j, 1.0);

        int cuts_added = 0;

        // For each visited customer as a source, see if it can reach t
        for (int src = 0; src < N && cuts_added < MAX_CUTS_PER_CALL; ++src) {
            if (src == s || src == t) continue;
            if (Y[src] < 0.5) continue;

            MF local = mf;
            (void)local.maxflow(src, t);
            auto inS = local.sourceSideReachable(src);

            // If t is not reachable from src, then δ+(S) has capacity 0;
            // pick h in S with y_h=1 and add lazy cut: sum δ+(S) x >= y_h (=1)
            if (!inS[t]) {
                int h = -1;
                for (int v = 0; v < N; ++v) {
                    if (inS[v] && v != s && v != t && Y[v] > 0.5) { h = v; break; }
                }
                if (h != -1) {
                    add_connectivity_cut(inS, h, /*asLazy=*/true);
                    ++cuts_added;
                }
            }
        }
        return cuts_added;
    }

    void callback() override {
        // DON'T do: int where = getIntInfo(GRB_CB_WHERE);

        if (where == GRB_CB_MIPNODE) {
            // Fractional LP node: separate user cuts
            int stat = getIntInfo(GRB_CB_MIPNODE_STATUS);
            if (stat == GRB_OPTIMAL) {
                separate_fractional();
            }
        } else if (where == GRB_CB_MIPSOL) {
            // Integer incumbent: enforce CCs as lazy constraints
            separate_integer();
        }
    }

};

// model.set(GRB_IntParam_LazyConstraints, 1);
// BCCallback cb(N, s, t, x, y);
// model.setCallback(&cb);