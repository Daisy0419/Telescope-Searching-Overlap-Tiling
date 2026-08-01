#pragma once
#include <ilconcert/ilomodel.h>
#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <queue>

ILOSTLBEGIN

struct ProblemDataTeam {
    int N;                                      // #nodes
    int K;                                      // #teams
    int startIndex;                             // s
    int endIndex;                               // t
    double BudgetPerTeam;                       // per-team budget (scalar)
    std::vector<std::vector<double>> Cost;      // NxN
    std::vector<double> Prize;                  // N (typically Prize[s]=Prize[t]=0)

    ProblemDataTeam(int N_, int K_, int s_, int t_, double B_,
                      const std::vector<std::vector<double>>& C_,
                      const std::vector<double>& P_)
        : N(N_), K(K_), startIndex(s_), endIndex(t_), BudgetPerTeam(B_),
          Cost(C_), Prize(P_) {}
};

std::pair<std::vector<std::vector<int>>, double>
cplexSolveTeamST(const std::vector<std::vector<double>>& Cost,
                 const std::vector<double>& Prize,
                 int start, int end, double BudgetPerTeam,
                 int mTeams,
                 double mipGap, double timeLimitSec,
                 const std::vector<std::vector<int>>& initialRoutes={});


std::pair<std::vector<std::vector<int>>, double>
cplexSolveTeamSTBound(const std::vector<std::vector<double>>& Cost,
                 const std::vector<double>& Prize,
                 int start, int end, double BudgetPerTeam,
                 int mTeams,
                 double mipGap, double timeLimitSec,
                 const std::vector<std::vector<int>>& initialRoutes={});




class TopConnLazyI : public IloCplex::LazyConstraintCallbackI {
    const IloArray<IloNumVarArray>& x;
    const IloNumVarArray* y; // nullptr for OP
    int N, s, t;
    double eps;

public:
    ILOCOMMONCALLBACKSTUFF(TopConnLazy) // <-- macro base name 'TopConnLazy' => class must be TopConnLazyI

    TopConnLazyI(IloEnv env,
                 const IloArray<IloNumVarArray>& x_,
                 const IloNumVarArray* y_,
                 int N_, int s_, int t_, double eps_=1e-6)
        : IloCplex::LazyConstraintCallbackI(env),
          x(x_), y(y_), N(N_), s(s_), t(t_), eps(eps_) {}

    void bfsForward(int src, const std::vector<std::vector<double>>& xv,
                    std::vector<char>& vis){
        std::queue<int> q; q.push(src); vis[src]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v=0; v<N; ++v) if (u!=v && !vis[v] && xv[u][v]>eps){
                vis[v]=1; q.push(v);
            }
        }
    }

    void bfsReverse(int sink, const std::vector<std::vector<double>>& xv,
                    std::vector<char>& vis){
        std::queue<int> q; q.push(sink); vis[sink]=1;
        while(!q.empty()){
            int u=q.front(); q.pop();
            for(int v=0; v<N; ++v) if (u!=v && !vis[v] && xv[v][u]>eps){
                vis[v]=1; q.push(v);
            }
        }
    }

    void main() override {
        // read solution
        std::vector<std::vector<double>> xv(N, std::vector<double>(N,0.0));
        for (int i=0;i<N;++i) for (int j=0;j<N;++j) if (i!=j)
            xv[i][j] = getValue(x[i][j]);

        // (1) s->t cut
        std::vector<char> S(N,0);
        bfsForward(s, xv, S);
        if (!S[t]) {
            IloExpr lhs(getEnv());
            for (int i=0;i<N;++i) if (S[i])
                for (int j=0;j<N;++j) if (!S[j] && i!=j)
                    lhs += x[i][j];
            add(lhs >= 1.0);
            lhs.end();
            return;
        }

        // (2) visited-inside (TOP only)
        if (y) {
            for (int h=0; h<N; ++h) {
                if (h==s || h==t) continue;
                double yh = getValue((*y)[h]);
                if (yh <= eps) continue;

                std::vector<char> Sh(N,0);
                bfsReverse(h, xv, Sh);
                if (Sh[s] || Sh[t]) continue;

                double cap = 0.0;
                for (int i=0;i<N;++i) if (Sh[i])
                    for (int j=0;j<N;++j) if (!Sh[j] && i!=j)
                        cap += xv[i][j];

                if (cap + 1e-8 < yh) {
                    IloExpr lhs(getEnv());
                    for (int i=0;i<N;++i) if (Sh[i])
                        for (int j=0;j<N;++j) if (!Sh[j] && i!=j)
                            lhs += x[i][j];
                    add(lhs >= (*y)[h]);
                    lhs.end();
                    return;
                }
            }
        }
    }
};

// how to use
// IloCplex cplex(model);
// // ... set params ...
// cplex.use(new (env) TopConnLazyI(env, x, /* TOP: */ &y, /* OP: */ nullptr, N, s, t));

// Tiny max-flow for separation
struct CutFlow {
  struct Edge { int to, rev; double cap; };
  int n; std::vector<std::vector<Edge>> g;
  explicit CutFlow(int n_) : n(n_), g(n_) {}
  void addEdge(int u, int v, double c) {
    Edge a{v,(int)g[v].size(),c}, b{u,(int)g[u].size(),0.0};
    g[u].push_back(a); g[v].push_back(b);
  }
  double maxflow(int s, int t) {
    double flow = 0.0;
    std::vector<int> pV(n), pE(n);
    for (;;) {
      std::fill(pV.begin(), pV.end(), -1);
      std::queue<int> q; q.push(s); pV[s]=s; pE[s]=-1;
      while (!q.empty() && pV[t]==-1) {
        int u=q.front(); q.pop();
        for (int i=0;i<(int)g[u].size();++i) {
          auto &e=g[u][i];
          if (pV[e.to]==-1 && e.cap>1e-12) { pV[e.to]=u; pE[e.to]=i; q.push(e.to); if (e.to==t) break; }
        }
      }
      if (pV[t]==-1) break;
      double aug = std::numeric_limits<double>::infinity();
      for (int v=t; v!=s; v=pV[v]) aug = std::min(aug, g[pV[v]][pE[v]].cap);
      for (int v=t; v!=s; v=pV[v]) { auto &e=g[pV[v]][pE[v]]; auto &r=g[e.to][e.rev]; e.cap-=aug; r.cap+=aug; }
      flow += aug;
    }
    return flow;
  }
  std::vector<char> sourceSide(int s) {
    std::vector<char> vis(n,0);
    std::queue<int> q; q.push(s); vis[s]=1;
    while (!q.empty()) {
      int u=q.front(); q.pop();
      for (auto &e:g[u]) if (!vis[e.to] && e.cap>1e-12) { vis[e.to]=1; q.push(e.to); }
    }
    return vis;
  }
};

// ---------- USER CUT: fractional LP nodes ----------
class TopUserCutCB : public IloCplex::UserCutCallbackI {
  IloArray<IloNumVarArray> x;
  IloNumVarArray           y;
  int N, s, t;
  double EPS = 1e-7;
  int    MAX_CUTS = 16;

public:
  TopUserCutCB(IloEnv env,
               const IloArray<IloNumVarArray>& x_,
               const IloNumVarArray& y_,
               int N_, int s_, int t_)
  : IloCplex::UserCutCallbackI(env), x(x_), y(y_), N(N_), s(s_), t(t_) {}

  // Copy ctor required for duplicateCallback
  TopUserCutCB(const TopUserCutCB& cb)
  : IloCplex::UserCutCallbackI(cb),
    x(cb.x), y(cb.y), N(cb.N), s(cb.s), t(cb.t),
    EPS(cb.EPS), MAX_CUTS(cb.MAX_CUTS) {}

  IloCplex::CallbackI* duplicateCallback() const override {
    return (new (getEnv()) TopUserCutCB(*this));
  }

  void addConnectivityCut(const std::vector<char>& inS, int h) {
    IloEnv env = getEnv();
    IloExpr lhs(env);
    for (int p=0;p<N;++p) {
      if (!inS[p] || p==s || p==t) continue; // customers only
      for (int q=0;q<N;++q) if (q!=p && !inS[q]) lhs += x[p][q];
    }
    lhs -= y[h];
    add(lhs >= 0.0);
    lhs.end();
  }

  void main() override {
    // fractional values
    std::vector<std::vector<double>> X(N, std::vector<double>(N,0.0));
    std::vector<double> Y(N,0.0);
    for (int i=0;i<N;++i) {
      Y[i] = getValue(y[i]);
      for (int j=0;j<N;++j) if (i!=j) X[i][j] = getValue(x[i][j]);
    }

    // support graph
    CutFlow base(N);
    for (int i=0;i<N;++i)
      for (int j=0;j<N;++j) if (i!=j && X[i][j] > EPS) base.addEdge(i,j,X[i][j]);

    int added = 0;
    for (int src=0; src<N && added<MAX_CUTS; ++src) {
      if (src==s || src==t) continue;
      CutFlow G = base;
      (void)G.maxflow(src, t);
      auto inS = G.sourceSide(src);

      double maxY = 0.0; int h = -1;
      for (int v=0; v<N; ++v)
        if (inS[v] && v!=s && v!=t && Y[v] > maxY + EPS) { maxY = Y[v]; h = v; }
      if (h==-1) continue;

      double cutCap = 0.0;
      for (int p=0;p<N;++p) {
        if (!inS[p] || p==s || p==t) continue;
        for (int q=0;q<N;++q) if (q!=p && !inS[q]) cutCap += X[p][q];
      }
      if (cutCap + 1e-6 < maxY) { addConnectivityCut(inS, h); ++added; }
    }
  }
};

// ---------- LAZY: integer incumbents ----------
class TopLazyCB : public IloCplex::LazyConstraintCallbackI {
  IloArray<IloNumVarArray> x;
  IloNumVarArray           y;
  int N, s, t;

public:
  TopLazyCB(IloEnv env,
            const IloArray<IloNumVarArray>& x_,
            const IloNumVarArray& y_,
            int N_, int s_, int t_)
  : IloCplex::LazyConstraintCallbackI(env), x(x_), y(y_), N(N_), s(s_), t(t_) {}

  TopLazyCB(const TopLazyCB& cb)
  : IloCplex::LazyConstraintCallbackI(cb),
    x(cb.x), y(cb.y), N(cb.N), s(cb.s), t(cb.t) {}

  IloCplex::CallbackI* duplicateCallback() const override {
    return (new (getEnv()) TopLazyCB(*this));
  }

  void addConnectivityCut(const std::vector<char>& inS, int h) {
    IloEnv env = getEnv();
    IloExpr lhs(env);
    for (int p=0;p<N;++p) {
      if (!inS[p] || p==s || p==t) continue;
      for (int q=0;q<N;++q) if (q!=p && !inS[q]) lhs += x[p][q];
    }
    lhs -= y[h];
    add(lhs >= 0.0);
    lhs.end();
  }

  void main() override {
    // incumbent values
    std::vector<std::vector<double>> X(N, std::vector<double>(N,0.0));
    std::vector<double> Y(N,0.0);
    for (int i=0;i<N;++i) {
      Y[i] = getValue(y[i]);
      for (int j=0;j<N;++j) if (i!=j) X[i][j] = getValue(x[i][j]);
    }

    // selected-arc graph
    CutFlow base(N);
    for (int i=0;i<N;++i)
      for (int j=0;j<N;++j) if (i!=j && X[i][j] > 0.5) base.addEdge(i,j,1.0);

    for (int src=0; src<N; ++src) {
      if (src==s || src==t) continue;
      if (Y[src] < 0.5)     continue;

      CutFlow G = base;
      (void)G.maxflow(src, t);
      auto inS = G.sourceSide(src);

      if (!inS[t]) {
        int h=-1;
        for (int v=0; v<N; ++v)
          if (inS[v] && v!=s && v!=t && Y[v] > 0.5) { h=v; break; }
        if (h!=-1) addConnectivityCut(inS, h);
      }
    }
  }
};