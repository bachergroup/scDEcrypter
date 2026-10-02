#include <Rcpp.h>
#include <vector>
#include <thread>
#include <algorithm>
#include <cmath>
using namespace Rcpp;

// Weighted sufficient statistics for batched permutations, without BLAS.
//
// Y:   n x p expression (column-major, genes in columns).
// Bg:  n x Kb matrix of permuted weights, Kb = (#columns per permutation) x (#permutations).
// Returns S0 (length Kb), S1 and S2 (Kb x p) with
//   S1[q, j] = sum_k Bg[k, q] * Y[k, j],   S2[q, j] = sum_k Bg[k, q] * Y[k, j]^2.
//
// The weights are transposed once to Kb x n so that, for each cell, the update over
// the Kb columns is a contiguous axpy-type loop with no reduction, which compilers
// vectorize at -O2. Genes are processed in small blocks so each pass over the
// weights serves several genes, and cells are blocked so the weight block stays in
// cache. Work is split over genes across std::thread workers; no R API is touched
// inside the workers.

namespace {

struct Job {
  const double* Y; int n; int p;
  const double* Bt; int Kb;            // Bt is Kb x n (column-major => cell-major rows of length Kb)
  double* S1; double* S2;              // Kb x p
  int g0, g1;
};

void run_job(const Job& jb) {
  const int GB = 8;        // genes per block
  const int NB = 1024;     // cells per block
  const int Kb = jb.Kb, n = jb.n;
  std::vector<double> acc1(static_cast<size_t>(GB) * Kb), acc2(static_cast<size_t>(GB) * Kb);

  for (int g = jb.g0; g < jb.g1; g += GB) {
    const int ng = std::min(GB, jb.g1 - g);
    std::fill(acc1.begin(), acc1.end(), 0.0);
    std::fill(acc2.begin(), acc2.end(), 0.0);
    for (int k0 = 0; k0 < n; k0 += NB) {
      const int k1 = std::min(n, k0 + NB);
      for (int gi = 0; gi < ng; ++gi) {
        const double* y = jb.Y + static_cast<size_t>(g + gi) * n;
        double* a1 = acc1.data() + static_cast<size_t>(gi) * Kb;
        double* a2 = acc2.data() + static_cast<size_t>(gi) * Kb;
        for (int k = k0; k < k1; ++k) {
          const double yy = y[k];
          if (yy == 0.0) continue;                 // sparse-friendly: zeros contribute nothing
          const double y2 = yy * yy;
          const double* w = jb.Bt + static_cast<size_t>(k) * Kb;
          for (int q = 0; q < Kb; ++q) {
            a1[q] += w[q] * yy;
            a2[q] += w[q] * y2;
          }
        }
      }
    }
    for (int gi = 0; gi < ng; ++gi) {
      std::copy(acc1.begin() + static_cast<size_t>(gi) * Kb, acc1.begin() + static_cast<size_t>(gi + 1) * Kb,
                jb.S1 + static_cast<size_t>(g + gi) * Kb);
      std::copy(acc2.begin() + static_cast<size_t>(gi) * Kb, acc2.begin() + static_cast<size_t>(gi + 1) * Kb,
                jb.S2 + static_cast<size_t>(g + gi) * Kb);
    }
  }
}

} // namespace

// [[Rcpp::export]]
List perm_suffstats_rcpp(const NumericMatrix& Y, const NumericMatrix& Bg, int threads) {
  const int n = Y.nrow(), p = Y.ncol(), Kb = Bg.ncol();
  if (Bg.nrow() != n) stop("nrow(Bg) must equal nrow(Y).");
  if (threads < 1) threads = 1;

  // transpose weights to Kb x n
  std::vector<double> Bt(static_cast<size_t>(Kb) * n);
  for (int q = 0; q < Kb; ++q) {
    const double* col = Bg.begin() + static_cast<size_t>(q) * n;
    for (int k = 0; k < n; ++k) Bt[static_cast<size_t>(k) * Kb + q] = col[k];
  }
  NumericVector S0(Kb);
  for (int q = 0; q < Kb; ++q) {
    double s = 0.0; const double* col = Bg.begin() + static_cast<size_t>(q) * n;
    for (int k = 0; k < n; ++k) s += col[k];
    S0[q] = s;
  }
  NumericMatrix S1(Kb, p), S2(Kb, p);

  const int T = std::min(threads, std::max(1, p / 8));
  std::vector<std::thread> pool;
  std::vector<Job> jobs(T);
  const int per = (p + T - 1) / T;
  for (int t = 0; t < T; ++t) {
    jobs[t] = Job{Y.begin(), n, p, Bt.data(), Kb, S1.begin(), S2.begin(),
                  std::min(p, t * per), std::min(p, (t + 1) * per)};
  }
  if (T == 1) {
    run_job(jobs[0]);
  } else {
    for (int t = 0; t < T; ++t) pool.emplace_back(run_job, std::cref(jobs[t]));
    for (auto& th : pool) th.join();
  }
  return List::create(Named("S0") = S0, Named("S1") = S1, Named("S2") = S2);
}
