#include <Rcpp.h>
#include <vector>
#include <thread>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

namespace {

// Log-sum-exp over the entries of L at the given indices.
inline double lse_idx(const std::vector<double>& L, const std::vector<int>& ids) {
  double mx = R_NegInf;
  for (int q : ids) if (R_finite(L[q]) && L[q] > mx) mx = L[q];
  if (!R_finite(mx)) return R_NegInf;
  double tot = 0.0;
  for (int q : ids) if (R_finite(L[q])) tot += std::exp(L[q] - mx);
  return mx + std::log(tot);
}

struct PanelJob {
  const double* LbT; const double* y; const int* perm; int n; int K; int C; int V;
  const double* mu; const double* sigma2;
  const int* c_obs; const int* v_obs;
  const int* comp; int ncomp; double var_floor;
  double* out;            // C x B
  int b0, b1;
};

void run_panel(const PanelJob& jb) {
  const int n = jb.n, K = jb.K, C = jb.C, V = jb.V;
  const double log2pi = std::log(2.0 * M_PI);
  std::vector<double> lognorm(K), inv2s(K);
  for (int q = 0; q < K; ++q) {
    lognorm[q] = -0.5 * (log2pi + std::log(jb.sigma2[q]));
    inv2s[q] = 0.5 / jb.sigma2[q];
  }
  std::vector<int> all_ids(K);
  for (int q = 0; q < K; ++q) all_ids[q] = q;
  std::vector<double> L(K), S0(K), S1(K), S2(K);
  std::vector<int> ids;

  for (int b = jb.b0; b < jb.b1; ++b) {
    std::fill(S0.begin(), S0.end(), 0.0);
    std::fill(S1.begin(), S1.end(), 0.0);
    std::fill(S2.begin(), S2.end(), 0.0);
    const int* pb = jb.perm + static_cast<size_t>(b) * n;   // 1-based source cell for each k
    for (int k = 0; k < n; ++k) {
      const double y = jb.y[pb[k] - 1];
      const double* Lk = jb.LbT + static_cast<size_t>(k) * K;   // cell k, contiguous over states
      const bool ck = jb.c_obs[k] != NA_INTEGER;
      const bool vk = jb.v_obs[k] != NA_INTEGER;
      if (ck && vk) {
        const int q = (jb.c_obs[k] - 1) + C * (jb.v_obs[k] - 1);
        S0[q] += 1.0; S1[q] += y; S2[q] += y * y;
        continue;
      }
      // only the states this cell can occupy are evaluated
      if (ck) {
        ids.assign(V, 0);
        for (int v = 0; v < V; ++v) ids[v] = (jb.c_obs[k] - 1) + C * v;
      } else if (vk) {
        ids.assign(C, 0);
        for (int c = 0; c < C; ++c) ids[c] = c + C * (jb.v_obs[k] - 1);
      } else {
        ids = all_ids;
      }
      for (int q : ids) {
        const double d = y - jb.mu[q];
        L[q] = Lk[q] + lognorm[q] - d * d * inv2s[q];
      }
      const double l = lse_idx(L, ids);
      if (!R_finite(l)) continue;
      for (int q : ids) {
        if (!R_finite(L[q])) continue;
        const double wq = std::exp(L[q] - l);
        S0[q] += wq; S1[q] += wq * y; S2[q] += wq * y * y;
      }
    }
    for (int c = 0; c < C; ++c) {
      double sum_n = 0.0, sum_s1 = 0.0;
      for (int i = 0; i < jb.ncomp; ++i) {
        const int q = c + C * (jb.comp[i] - 1);
        sum_n += S0[q]; sum_s1 += S1[q];
      }
      const double mu0 = sum_n > 0.0 ? sum_s1 / sum_n : 0.0;
      double lrt = 0.0;
      for (int i = 0; i < jb.ncomp; ++i) {
        const int q = c + C * (jb.comp[i] - 1);
        const double nn = S0[q];
        if (!(nn > 0.0)) continue;
        const double m1 = S1[q] / nn;
        const double ss1 = std::max(S2[q] - m1 * S1[q], 0.0);
        const double sig1 = std::max(ss1 / nn, jb.var_floor);
        const double ss0 = std::max(S2[q] - 2.0 * mu0 * S1[q] + mu0 * mu0 * nn, 0.0);
        const double sig0 = std::max(ss0 / nn, jb.var_floor);
        lrt += ss0 / sig0 + nn * std::log(sig0) - ss1 / sig1 - nn * std::log(sig1);
      }
      jb.out[c + static_cast<size_t>(C) * b] = lrt;
    }
  }
}

} // namespace

// Permutation null for a gene that entered the test-set weights.
//
// Lbase: n x (C*V) log-numerators of the E-step (including log pi) with the
//        contribution of the gene under test removed; column index c + C*(v-1) (0-based).
// y:     length-n expression of the gene; perm: n x B integer matrix of 1-based
//        source cells, so column b of the permuted gene is y[perm[, b]].
// mu, sigma2: length C*V fitted mean and variance of the gene, same ordering.
// c_obs, v_obs: 1-based observed labels (NA_INTEGER if unobserved), as in E_step_rcpp.
// Returns a C x B matrix of LRT statistics, one per cell type and permutation,
// computed exactly as de_from_suffstats() does for the compared statuses.
// Permutations are split across std::thread workers.
// [[Rcpp::export]]
NumericMatrix panel_perm_lrt_rcpp(const NumericMatrix& Lbase,
                                  const NumericVector& y,
                                  const IntegerMatrix& perm,
                                  const NumericVector& mu,
                                  const NumericVector& sigma2,
                                  const IntegerVector& c_obs,
                                  const IntegerVector& v_obs,
                                  int C, int V,
                                  const IntegerVector& comp_idx,
                                  double var_floor,
                                  int threads = 1) {
  const int n = Lbase.nrow();
  const int K = C * V;
  const int B = perm.ncol();
  if (Lbase.ncol() != K) stop("ncol(Lbase) must equal C*V.");
  if (perm.nrow() != n || y.size() != n) stop("y and perm must match nrow(Lbase).");
  if (c_obs.size() != n || v_obs.size() != n) stop("Label lengths must match nrow(Lbase).");
  if (threads < 1) threads = 1;

  // transpose the log-numerators once to K x n so each cell's states are contiguous
  std::vector<double> LbT(static_cast<size_t>(K) * n);
  for (int q = 0; q < K; ++q) {
    const double* col = Lbase.begin() + static_cast<size_t>(q) * n;
    for (int k = 0; k < n; ++k) LbT[static_cast<size_t>(k) * K + q] = col[k];
  }
  NumericMatrix out(C, B);
  const int T = std::min(threads, B);
  const int per = (B + T - 1) / T;
  std::vector<PanelJob> jobs(T);
  for (int t = 0; t < T; ++t) {
    jobs[t] = PanelJob{LbT.data(), y.begin(), perm.begin(), n, K, C, V, mu.begin(), sigma2.begin(),
                       c_obs.begin(), v_obs.begin(), comp_idx.begin(), static_cast<int>(comp_idx.size()),
                       var_floor, out.begin(), std::min(B, t * per), std::min(B, (t + 1) * per)};
  }
  if (T == 1) {
    run_panel(jobs[0]);
  } else {
    std::vector<std::thread> pool;
    for (int t = 0; t < T; ++t) pool.emplace_back(run_panel, std::cref(jobs[t]));
    for (auto& th : pool) th.join();
  }
  return out;
}
