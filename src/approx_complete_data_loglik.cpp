#include <Rcpp.h>
using namespace Rcpp;

inline int idx3(int i, int j, int k, int d1, int d2) {
  return i + d1 * (j + d2 * k);   // R column-major indexing
}

inline NumericVector init_array3(int d1, int d2, int d3) {
  NumericVector out(d1 * d2 * d3);
  out.attr("dim") = IntegerVector::create(d1, d2, d3);
  return out;
}

inline NumericMatrix mean_weights_over_status(const NumericVector& W,
                                              const IntegerVector& dimW,
                                              const std::vector<int>& status_idx) {
  const int n_cells = dimW[0];
  const int c_dim = dimW[1];
  NumericMatrix weights(n_cells, c_dim);
  const double scale = 1.0 / static_cast<double>(status_idx.size());

  for (int vv : status_idx) {
    for (int cc = 0; cc < c_dim; ++cc) {
      for (int ii = 0; ii < n_cells; ++ii) {
        weights(ii, cc) += W[idx3(ii, cc, vv, dimW[0], dimW[1])] * scale;
      }
    }
  }

  return weights;
}

inline void fill_means_for_statuses(const NumericMatrix& Y,
                                    const NumericMatrix& weights,
                                    const std::vector<int>& status_idx,
                                    NumericVector& out,
                                    const IntegerVector& dimOut) {
  const int n_cells = Y.nrow();
  const int n_genes = Y.ncol();
  const int c_dim = weights.ncol();

  NumericVector denom(c_dim);
  for (int cc = 0; cc < c_dim; ++cc) {
    double total = 0.0;
    for (int ii = 0; ii < n_cells; ++ii) {
      total += weights(ii, cc);
    }
    denom[cc] = total;
  }

  for (int vv : status_idx) {
    for (int cc = 0; cc < c_dim; ++cc) {
      const double denom_cc = denom[cc];
      for (int gg = 0; gg < n_genes; ++gg) {
        double numer = 0.0;
        for (int ii = 0; ii < n_cells; ++ii) {
          numer += weights(ii, cc) * Y(ii, gg);
        }
        out[idx3(gg, cc, vv, dimOut[0], dimOut[1])] = numer / denom_cc;
      }
    }
  }
}

inline std::vector<int> complement_status_idx(int v_dim, const IntegerVector& comp_idx) {
  std::vector<bool> selected(v_dim, false);
  for (int idx : comp_idx) {
    selected[idx - 1] = true;
  }

  std::vector<int> rest_idx;
  rest_idx.reserve(v_dim - comp_idx.size());
  for (int vv = 0; vv < v_dim; ++vv) {
    if (!selected[vv]) {
      rest_idx.push_back(vv);
    }
  }

  return rest_idx;
}

// [[Rcpp::export]]
NumericMatrix approx_complete_data_loglik_rcpp(
    const NumericMatrix& Y,        // n_cells x n_genes
    const NumericVector& M,        // dim n_genes x c_dim x v_dim
    const NumericVector& W,        // dim n_cells x c_dim x v_dim
    const NumericVector& sigma2) { // dim n_genes x c_dim x v_dim

  IntegerVector dimM = M.attr("dim");
  IntegerVector dimW = W.attr("dim");
  IntegerVector dimS = sigma2.attr("dim");

  if (dimM.size() != 3 || dimW.size() != 3 || dimS.size() != 3) {
    stop("M, W, sigma2 must be 3D arrays.");
  }

  const int n_cells = Y.nrow();
  const int n_genes = Y.ncol();
  const int c_dim   = dimM[1];
  const int v_dim   = dimM[2];

  if (dimM[0] != n_genes) stop("dim(M)[1] must equal ncol(Y).");
  if (dimW[0] != n_cells) stop("dim(W)[1] must equal nrow(Y).");
  if (dimW[1] != c_dim || dimW[2] != v_dim) stop("W dims mismatch M.");
  if (dimS[0] != n_genes || dimS[1] != c_dim || dimS[2] != v_dim) stop("sigma2 dims mismatch M.");

  NumericMatrix out(n_genes, c_dim);
  const double log_2pi = std::log(2.0 * M_PI);
  const double* y_ptr = Y.begin();
  const double* m_ptr = M.begin();
  const double* w_ptr = W.begin();
  const double* s_ptr = sigma2.begin();

  for (int cc = 0; cc < c_dim; ++cc) {
    for (int kk = 0; kk < n_genes; ++kk) {
      double acc_k = 0.0;

      for (int vv = 0; vv < v_dim; ++vv) {
        const double mu = m_ptr[idx3(kk, cc, vv, dimM[0], dimM[1])];
        const double s2 = s_ptr[idx3(kk, cc, vv, dimS[0], dimS[1])];
        if (!(s2 > 0.0) || !R_finite(s2)) continue;
        const double inv_s2 = 1.0 / s2;
        const double log_norm_const = -0.5 * (log_2pi + std::log(s2));

        double sum_ll = 0.0;
        for (int ii = 0; ii < n_cells; ++ii) {
          const double x = y_ptr[ii + n_cells * kk];
          const double w = w_ptr[idx3(ii, cc, vv, dimW[0], dimW[1])];
          const double diff = x - mu;
          const double ll = (log_norm_const - 0.5 * diff * diff * inv_s2) * w;
          if (R_finite(ll)) sum_ll += ll;   // matches finite-filter logic
        }
        acc_k += sum_ll;
      }

      out(kk, cc) = acc_k;
    }
  }

  return out;
}

// [[Rcpp::export]]
List approx_complete_data_loglik_pair_rcpp(
    const NumericMatrix& Y,
    const NumericVector& M_alt,
    const NumericVector& sigma2_alt,
    const NumericVector& M_null,
    const NumericVector& sigma2_null,
    const NumericVector& W) {

  IntegerVector dimM_alt = M_alt.attr("dim");
  IntegerVector dimS_alt = sigma2_alt.attr("dim");
  IntegerVector dimM_null = M_null.attr("dim");
  IntegerVector dimS_null = sigma2_null.attr("dim");
  IntegerVector dimW = W.attr("dim");

  if (dimM_alt.size() != 3 || dimS_alt.size() != 3 || dimM_null.size() != 3 || dimS_null.size() != 3 || dimW.size() != 3) {
    stop("All inputs must be 3D arrays except Y.");
  }

  const int n_cells = Y.nrow();
  const int n_genes = Y.ncol();
  const int c_dim = dimM_alt[1];
  const int v_dim = dimM_alt[2];

  if (dimM_alt[0] != n_genes || dimS_alt[0] != n_genes || dimM_null[0] != n_genes || dimS_null[0] != n_genes) {
    stop("Gene dimensions of M/sigma2 must match ncol(Y).");
  }
  if (dimM_alt[1] != c_dim || dimM_alt[2] != v_dim || dimS_alt[1] != c_dim || dimS_alt[2] != v_dim ||
      dimM_null[1] != c_dim || dimM_null[2] != v_dim || dimS_null[1] != c_dim || dimS_null[2] != v_dim) {
    stop("Dimensions of M and sigma2 must match between alternative and null models.");
  }
  if (dimW[0] != n_cells || dimW[1] != c_dim || dimW[2] != v_dim) {
    stop("W dims must match Y and model dimensions.");
  }

  NumericMatrix ll_alt(n_genes, c_dim);
  NumericMatrix ll_null(n_genes, c_dim);
  const double log_2pi = std::log(2.0 * M_PI);
  const double* y_ptr = Y.begin();
  const double* m_alt_ptr = M_alt.begin();
  const double* s_alt_ptr = sigma2_alt.begin();
  const double* m_null_ptr = M_null.begin();
  const double* s_null_ptr = sigma2_null.begin();
  const double* w_ptr = W.begin();

  for (int cc = 0; cc < c_dim; ++cc) {
    for (int kk = 0; kk < n_genes; ++kk) {
      double acc_alt = 0.0;
      double acc_null = 0.0;

      for (int vv = 0; vv < v_dim; ++vv) {
        const double mu_alt = m_alt_ptr[idx3(kk, cc, vv, dimM_alt[0], dimM_alt[1])];
        const double s2_alt = s_alt_ptr[idx3(kk, cc, vv, dimS_alt[0], dimS_alt[1])];
        const double mu_null = m_null_ptr[idx3(kk, cc, vv, dimM_null[0], dimM_null[1])];
        const double s2_null = s_null_ptr[idx3(kk, cc, vv, dimS_null[0], dimS_null[1])];

        if (s2_alt > 0.0 && R_finite(s2_alt)) {
          const double inv_s2_alt = 1.0 / s2_alt;
          const double log_norm_const_alt = -0.5 * (log_2pi + std::log(s2_alt));
          double sum_alt = 0.0;
          for (int ii = 0; ii < n_cells; ++ii) {
            const double x = y_ptr[ii + n_cells * kk];
            const double w = w_ptr[idx3(ii, cc, vv, dimW[0], dimW[1])];
            const double diff = x - mu_alt;
            const double ll = (log_norm_const_alt - 0.5 * diff * diff * inv_s2_alt) * w;
            if (R_finite(ll)) sum_alt += ll;
          }
          acc_alt += sum_alt;
        }

        if (s2_null > 0.0 && R_finite(s2_null)) {
          const double inv_s2_null = 1.0 / s2_null;
          const double log_norm_const_null = -0.5 * (log_2pi + std::log(s2_null));
          double sum_null = 0.0;
          for (int ii = 0; ii < n_cells; ++ii) {
            const double x = y_ptr[ii + n_cells * kk];
            const double w = w_ptr[idx3(ii, cc, vv, dimW[0], dimW[1])];
            const double diff = x - mu_null;
            const double ll = (log_norm_const_null - 0.5 * diff * diff * inv_s2_null) * w;
            if (R_finite(ll)) sum_null += ll;
          }
          acc_null += sum_null;
        }
      }

      ll_alt(kk, cc) = acc_alt;
      ll_null(kk, cc) = acc_null;
    }
  }

  return List::create(
    Named("ll.alternative") = ll_alt,
    Named("ll.null") = ll_null
  );
}

// [[Rcpp::export]]
NumericVector de_mu_rcpp(const NumericMatrix& Y,
                         const NumericVector& W,
                         const IntegerVector& comp_idx) {
  IntegerVector dimW = W.attr("dim");
  if (dimW.size() != 3) {
    stop("W must be a 3D array.");
  }

  const int n_genes = Y.ncol();
  const int c_dim = dimW[1];
  const int v_dim = dimW[2];
  NumericVector out = init_array3(n_genes, c_dim, v_dim);
  IntegerVector dimOut = out.attr("dim");

  std::vector<int> rest_idx = complement_status_idx(v_dim, comp_idx);

  for (int raw_idx : comp_idx) {
    const int vv = raw_idx - 1;
    NumericMatrix weights(Y.nrow(), c_dim);
    for (int cc = 0; cc < c_dim; ++cc) {
      for (int ii = 0; ii < Y.nrow(); ++ii) {
        weights(ii, cc) = W[idx3(ii, cc, vv, dimW[0], dimW[1])];
      }
    }
    fill_means_for_statuses(Y, weights, std::vector<int>{vv}, out, dimOut);
  }

  if (!rest_idx.empty()) {
    NumericMatrix rest_weights = mean_weights_over_status(W, dimW, rest_idx);
    fill_means_for_statuses(Y, rest_weights, rest_idx, out, dimOut);
  }

  return out;
}

// [[Rcpp::export]]
NumericVector de_mu_null_rcpp(const NumericMatrix& Y,
                              const NumericVector& W,
                              const IntegerVector& comp_idx) {
  IntegerVector dimW = W.attr("dim");
  if (dimW.size() != 3) {
    stop("W must be a 3D array.");
  }

  const int n_genes = Y.ncol();
  const int c_dim = dimW[1];
  const int v_dim = dimW[2];
  NumericVector out = init_array3(n_genes, c_dim, v_dim);
  IntegerVector dimOut = out.attr("dim");

  std::vector<int> comp_zero_based;
  comp_zero_based.reserve(comp_idx.size());
  for (int idx : comp_idx) {
    comp_zero_based.push_back(idx - 1);
  }
  std::vector<int> rest_idx = complement_status_idx(v_dim, comp_idx);

  NumericMatrix shared_weights = mean_weights_over_status(W, dimW, comp_zero_based);
  fill_means_for_statuses(Y, shared_weights, comp_zero_based, out, dimOut);

  if (!rest_idx.empty()) {
    NumericMatrix rest_weights = mean_weights_over_status(W, dimW, rest_idx);
    fill_means_for_statuses(Y, rest_weights, rest_idx, out, dimOut);
  }

  return out;
}

// [[Rcpp::export]]
NumericVector de_sigma2_rcpp(const NumericMatrix& Y,
                             const NumericVector& W,
                             const NumericVector& M) {
  IntegerVector dimW = W.attr("dim");
  IntegerVector dimM = M.attr("dim");
  if (dimW.size() != 3 || dimM.size() != 3) {
    stop("W and M must be 3D arrays.");
  }
  if (dimM[0] != Y.ncol() || dimW[0] != Y.nrow() || dimM[1] != dimW[1] || dimM[2] != dimW[2]) {
    stop("Dimensions of Y, W, and M do not align.");
  }

  NumericVector out = init_array3(dimM[0], dimM[1], dimM[2]);

  for (int vv = 0; vv < dimM[2]; ++vv) {
    for (int cc = 0; cc < dimM[1]; ++cc) {
      double denom = 0.0;
      for (int ii = 0; ii < Y.nrow(); ++ii) {
        denom += W[idx3(ii, cc, vv, dimW[0], dimW[1])];
      }

      for (int gg = 0; gg < Y.ncol(); ++gg) {
        double numer = 0.0;
        const double mu = M[idx3(gg, cc, vv, dimM[0], dimM[1])];
        for (int ii = 0; ii < Y.nrow(); ++ii) {
          const double diff = Y(ii, gg) - mu;
          numer += W[idx3(ii, cc, vv, dimW[0], dimW[1])] * diff * diff;
        }

        double value = numer / denom;
        if (!R_finite(value) || value < 0.1) {
          value = 0.1;
        }
        out[idx3(gg, cc, vv, dimM[0], dimM[1])] = value;
      }
    }
  }

  return out;
}