#' Permutation-Calibrated Differential Expression Test
#'
#' Perform differential expression testing on test data.
#' The likelihood ratio statistic is computed in closed form from the weighted
#' sums of \eqn{y} and \eqn{y^2} per cell type and viral status, permutations
#' are batched into a small number of large matrix products, and (optionally)
#' the genes used to fit the model are handled by permuting the gene and
#' recomputing the test-set weights, so that their null distribution reflects
#' the dependence of the weights on the gene under test.
#'
#' @param mod_results A list containing model fitting results from \code{\link{fit_scDEcrypter}}.
#' @param testData A Seurat object or list containing test data.
#' @param testingGenes Character vector of gene names to test.
#' @param compGroups Vector (names or indices) of the viral statuses to compare.
#' @param method Testing method: \code{"permutation"} (default; \code{"perm"}
#'   is an alias) or \code{"chisq"} for the chi-square approximation.
#' @param nPerm Integer number of permutations. Default \code{1000}.
#' @param mc.cores Deprecated compatibility argument; sets the default for \code{threads}.
#' @param strata Optional vector of length equal to the number of test cells
#'   giving the permutation strata. Defaults to the observed partitioning
#'   labels, with unlabeled cells forming their own stratum. Supply, e.g.,
#'   \code{interaction(sample, celltype)} to also stratify on sample.
#' @param panelGenes How to treat tested genes that were also used to fit the
#'   model: \code{"recompute"} (default) permutes the gene within strata and
#'   recomputes the test-set weights for each permutation; \code{"exclude"}
#'   sets their p-values to \code{NA}; \code{"ignore"} treats them like any
#'   other gene (anti-conservative).
#' @param var_floor Minimum variance. Default \code{0.1}.
#' @param chunk Number of permutations per batch. Default \code{100}.
#' @param return_null Logical; if \code{TRUE}, return the permuted statistics
#'   (genes x cell types x nPerm), which allows pooled or tail-fitted p-values.
#' @param threads Integer number of C++ threads for the weighted sums and the
#'   panel-gene recomputation. Defaults to \code{mc.cores}.
#'
#' @return A list containing the likelihoods, test statistics, p-values, log
#'   fold changes, and testing method, plus
#'   \code{panel_genes} (the tested genes that were in the training panel) and,
#'   if requested, \code{null.stat}.
#' @export
deTest <- function(mod_results, testData, testingGenes, compGroups,
                        method = c("permutation", "perm", "chisq"),
                        nPerm = 1000, mc.cores = 1, strata = NULL,
                        panelGenes = c("recompute", "exclude", "ignore"),
                        var_floor = 0.1, chunk = 100, return_null = FALSE,
                        threads = mc.cores) {
  method <- match.arg(method)
  if (identical(method, "perm")) method <- "permutation"
  panelGenes <- match.arg(panelGenes)
  tc <- extract_test_components(testData)
  mc <- extract_model_components(mod_results)
  W <- mod_results$weights_test
  if (is.null(W)) {
    message("Calculating weights on test set...")
    W <- E_step(tc$Y_test[, rownames(mc$M), drop = FALSE], tc$c_obs, tc$v_obs,
                mc$M, mc$probs, mc$sigma2)
  }
  out <- deTest_core(Y_test = tc$Y_test, c_obs = tc$c_obs, v_obs = tc$v_obs, W = W,
                          M = mc$M, sigma2 = mc$sigma2, probs = mc$probs,
                          testingGenes = testingGenes, compGroups = compGroups,
                          nPerm = if (method == "chisq") 0L else nPerm,
                          strata = strata, panelGenes = if (method == "chisq") "ignore" else panelGenes,
                          var_floor = var_floor, chunk = chunk, return_null = return_null,
                          threads = threads)
  if (method == "chisq") {
    ## chi-square approximation with one degree of freedom
    out$pval <- apply(out$lrt.stat, 2, function(x) pchisq(x, df = 1, lower.tail = FALSE))
    out$adj.pval <- apply(out$pval, 2, function(x) p.adjust(x, method = "fdr"))
    dimnames(out$pval) <- dimnames(out$adj.pval) <- dimnames(out$lrt.stat)
    out$test.method <- "chisq"
  }
  out
}

#' @keywords internal
deTest_core <- function(Y_test, c_obs, v_obs, W, M, sigma2, probs,
                             testingGenes, compGroups, nPerm = 1000, strata = NULL,
                             panelGenes = "recompute", var_floor = 0.1, chunk = 100,
                             return_null = FALSE, threads = 1, .taus = NULL) {
  n <- dim(W)[1]; C <- dim(W)[2]; V <- dim(W)[3]; K <- C * V
  comp_idx <- resolve_comp_status_indices(W, compGroups)
  rest_idx <- setdiff(seq_len(V), comp_idx)
  ct_names <- dimnames(W)[[2]]; vs_names <- dimnames(W)[[3]]

  Y <- as.matrix(Y_test[, testingGenes, drop = FALSE]); storage.mode(Y) <- "double"
  p <- ncol(Y)

  if (is.null(strata)) {
    strata <- as.character(c_obs); strata[is.na(strata)] <- "unlabeled"
  }
  if (length(strata) != n) stop("strata must have one entry per test cell.")
  strata_idx <- split(seq_len(n), as.character(strata))

  ## observed statistic ------------------------------------------------------
  ss <- perm_suffstats_rcpp(Y, matrix(W, n, K), threads)
  obs <- de_from_suffstats(ss$S0, ss$S1, ss$S2, C, V, comp_idx, rest_idx, var_floor, full = TRUE)
  lrt <- obs$lrt
  dimnames(lrt) <- list(testingGenes, ct_names)
  dimnames(obs$M.alt) <- list(testingGenes, ct_names, vs_names)

  ## permutation of conditional viral-status weights --------------------------
  prep <- lapply(strata_idx, function(idx) perm_stratum_prep(W, idx, C, V, comp_idx))
  exceed <- matrix(0, p, C, dimnames = dimnames(lrt))
  null_stat <- if (return_null) array(NA_real_, c(p, C, nPerm), dimnames = list(testingGenes, ct_names, NULL)) else NULL

  ## genes in the training panel: precompute E-step log-numerators without each gene
  panel <- intersect(testingGenes, rownames(M))
  do_panel <- length(panel) > 0 && panelGenes == "recompute"
  if (do_panel) {
    Yp <- as.matrix(Y_test[, rownames(M), drop = FALSE]); storage.mode(Yp) <- "double"
    Mm <- matrix(M, nrow(M), K); Sm <- matrix(sigma2, nrow(M), K)
    logpi <- log(as.vector(probs)); logpi[!is.finite(logpi)] <- -Inf
    const <- -0.5 * nrow(M) * log(2 * pi) - 0.5 * colSums(log(Sm)) - 0.5 * colSums(Mm^2 / Sm) + logpi
    L_full <- -0.5 * ((Yp * Yp) %*% (1 / Sm) - 2 * Yp %*% (Mm / Sm))
    L_full <- sweep(L_full, 2, const, "+")
    Lbase_list <- lapply(panel, function(g) {
      gi <- match(g, rownames(M)); y <- Yp[, gi]
      Tj <- -0.5 * log(2 * pi) - 0.5 * rep(log(Sm[gi, ]), each = n) -
        (outer(y, Mm[gi, ], "-"))^2 / (2 * rep(Sm[gi, ], each = n))
      L_full - Tj
    })
    panel_ji <- match(panel, testingGenes); panel_gi <- match(panel, rownames(M))
  }

  done <- 0
  while (done < nPerm) {
    B <- min(chunk, nPerm - done)
    ## one set of within-stratum permutations per batch, shared by all genes
    taus <- lapply(seq_along(strata_idx), function(s) {
      if (is.null(.taus)) lapply(seq_len(B), function(b) sample.int(length(strata_idx[[s]])))
      else .taus[[s]][done + seq_len(B)]
    })
    A0 <- matrix(0, K, B); A1 <- matrix(0, K * B, p); A2 <- matrix(0, K * B, p)
    for (s in seq_along(prep)) {
      ps <- prep[[s]]
      if (length(ps$cols) == 0) next
      big <- perm_stratum_big(ps, B, taus[[s]])
      ss <- perm_suffstats_rcpp(Y[ps$idx, , drop = FALSE], big, threads)
      map <- as.vector(outer(ps$cols, (seq_len(B) - 1) * K, "+"))
      A0[map] <- A0[map] + ss$S0
      A1[map, ] <- A1[map, ] + ss$S1
      A2[map, ] <- A2[map, ] + ss$S2
      if (!is.null(ps$derived)) {
        ## pure stratum with all statuses compared: the dropped column is the
        ## stratum total minus the computed columns of the same cell type
        if (is.null(ps$tot)) { ps$tot <- stratum_totals(Y, ps$idx, threads); prep[[s]]$tot <- ps$tot }
        Knz <- length(ps$cols)
        for (b in seq_len(B)) {
          rows_b <- (b - 1) * Knz + seq_len(Knz)
          dmap <- (b - 1) * K + ps$derived
          A0[dmap] <- A0[dmap] + (ps$ns - sum(ss$S0[rows_b]))
          A1[dmap, ] <- A1[dmap, ] + (ps$tot$S1 - colSums(ss$S1[rows_b, , drop = FALSE]))
          A2[dmap, ] <- A2[dmap, ] + (ps$tot$S2 - colSums(ss$S2[rows_b, , drop = FALSE]))
        }
      }
    }
    for (b in seq_len(B)) {
      rows <- (b - 1) * K + seq_len(K)
      lrt_b <- de_from_suffstats(A0[, b], A1[rows, , drop = FALSE], A2[rows, , drop = FALSE],
                                 C, V, comp_idx, rest_idx, var_floor, full = FALSE)$lrt
      inc <- (lrt_b >= lrt)
      if (do_panel) inc[panel_ji, ] <- 0    # panel genes are counted in the recompute branch
      exceed <- exceed + inc
      if (return_null) null_stat[, , done + b] <- lrt_b
    }
    if (do_panel) {
      perm_idx <- matrix(seq_len(n), n, B)
      for (s in seq_along(strata_idx)) {
        idx <- strata_idx[[s]]
        for (b in seq_len(B)) perm_idx[idx, b] <- idx[taus[[s]][[b]]]
      }
      storage.mode(perm_idx) <- "integer"
      for (i in seq_along(panel)) {
        lrt_g <- panel_perm_lrt_rcpp(Lbase_list[[i]], Yp[, panel_gi[i]], perm_idx,
                                     Mm[panel_gi[i], ], Sm[panel_gi[i], ],
                                     as.integer(c_obs), as.integer(v_obs),
                                     C, V, as.integer(comp_idx), var_floor, threads)  # C x B
        exceed[panel_ji[i], ] <- exceed[panel_ji[i], ] + rowSums(lrt_g >= lrt[panel_ji[i], ])
        if (return_null) null_stat[panel_ji[i], , done + seq_len(B)] <- lrt_g
      }
    }
    done <- done + B
  }

  pval <- if (nPerm > 0) (exceed + 1) / (nPerm + 1) else exceed * NA_real_
  if (length(panel) > 0 && panelGenes == "exclude") pval[match(panel, testingGenes), ] <- NA
  pval.adjust <- apply(pval, 2, function(x) p.adjust(x, method = "fdr"))
  dimnames(pval.adjust) <- dimnames(pval)

  out <- list(
    ll.alternative = obs$ll.alt, ll.null = obs$ll.null, lrt.stat = lrt,
    pval = pval, adj.pval = pval.adjust,
    logfc = compute_logFC(obs$M.alt, compGroups),
    test.method = "permutation", panel_genes = panel
  )
  dimnames(out$ll.alternative) <- dimnames(out$ll.null) <- dimnames(lrt)
  if (return_null) out$null.stat <- null_stat
  out
}

## Closed-form log-likelihoods and LRT from weighted sufficient statistics.
## S0: length-K vector, S1/S2: K x p matrices, K = C*V indexed c + C*(v-1).
## Matches approx_complete_data_loglik_rcpp with the plug-in estimators of
## de_mu_rcpp / de_mu_null_rcpp / de_sigma2_rcpp (including the variance floor).
#' @keywords internal
de_from_suffstats <- function(S0, S1, S2, C, V, comp_idx, rest_idx, var_floor, full = FALSE) {
  p <- ncol(S1); nc <- length(comp_idx)
  log2pi <- log(2 * pi)
  ll_piece <- function(n_v, ss) {
    sig <- ss / n_v
    sig[!is.finite(sig) | sig < var_floor] <- var_floor
    ll <- -0.5 * n_v * log2pi - 0.5 * n_v * log(sig) - 0.5 * ss / sig
    ll[n_v == 0, ] <- 0
    ll
  }
  pooled_ss <- function(n_v, s1, s2, mu) {
    k <- length(n_v)
    pmax(s2 - 2 * s1 * rep(mu, each = k) + n_v * rep(mu^2, each = k), 0)
  }
  lrt <- matrix(0, p, C)
  if (full) {
    ll_alt <- ll_null <- matrix(0, p, C)
    M_alt <- array(0, c(p, C, V))
  }
  for (cc in seq_len(C)) {
    rows <- cc + C * (comp_idx - 1)
    n_v <- S0[rows]; s1 <- S1[rows, , drop = FALSE]; s2 <- S2[rows, , drop = FALSE]
    mu1 <- s1 / n_v; mu1[!is.finite(mu1)] <- 0
    ss1 <- pmax(s2 - mu1 * s1, 0)
    ll1 <- ll_piece(n_v, ss1)
    mu0 <- colSums(s1) / sum(n_v); mu0[!is.finite(mu0)] <- 0
    ll0 <- ll_piece(n_v, pooled_ss(n_v, s1, s2, mu0))
    lrt[, cc] <- -2 * colSums(ll0 - ll1)
    if (full) {
      ll_r <- 0; mu_r <- NULL
      if (length(rest_idx) > 0) {
        rrows <- cc + C * (rest_idx - 1)
        n_r <- S0[rrows]; s1r <- S1[rrows, , drop = FALSE]; s2r <- S2[rrows, , drop = FALSE]
        mu_r <- colSums(s1r) / sum(n_r); mu_r[!is.finite(mu_r)] <- 0
        ll_r <- colSums(ll_piece(n_r, pooled_ss(n_r, s1r, s2r, mu_r)))
        for (v in rest_idx) M_alt[, cc, v] <- mu_r
      }
      ll_alt[, cc] <- colSums(ll1) + ll_r
      ll_null[, cc] <- colSums(ll0) + ll_r
      for (i in seq_len(nc)) M_alt[, cc, comp_idx[i]] <- mu1[i, ]
    }
  }
  if (full) list(lrt = lrt, ll.alt = ll_alt, ll.null = ll_null, M.alt = M_alt) else list(lrt = lrt)
}

## Per-stratum pieces for the conditional-weight permutation: each cell keeps its
## cell-type mass R[k, c] = sum_v W[k, c, v]; the conditional weights omega are
## permuted among cells of the stratum. Only the compared statuses are needed.
#' @keywords internal
perm_stratum_prep <- function(W, idx, C, V, comp_idx) {
  ns <- length(idx)
  Ws <- W[idx, , , drop = FALSE]
  R <- matrix(rowSums(matrix(Ws, ns * C, V)), ns, C)          # ns x C
  cols_c <- rep(seq_len(C), times = length(comp_idx))
  cols_v <- rep(comp_idx, each = C)
  cols <- cols_c + C * (cols_v - 1)                             # K-indices of (c, v in comp)
  keep <- colSums(R)[cols_c] > 0
  cols <- cols[keep]; cols_c <- cols_c[keep]; cols_v <- cols_v[keep]
  Wc <- matrix(Ws, ns, C * V)[, cols, drop = FALSE]
  Rfac <- R[, cols_c, drop = FALSE]
  omega <- Wc / Rfac; omega[!is.finite(omega)] <- 0
  derived <- NULL
  ## If every cell in the stratum has all its mass on one cell type and all viral
  ## statuses are compared, the conditional weights sum to one across the compared
  ## columns, so the last column's sums follow from the stratum totals.
  if (length(comp_idx) == V && length(cols) == V && all(abs(Rfac - 1) < 1e-8)) {
    derived <- cols[V]
    cols <- cols[-V]; Rfac <- Rfac[, -V, drop = FALSE]; omega <- omega[, -V, drop = FALSE]
  }
  list(idx = idx, ns = ns, cols = cols, Rfac = Rfac, omega = omega, derived = derived, tot = NULL)
}

#' @keywords internal
stratum_totals <- function(Y, idx, threads) {
  ss <- perm_suffstats_rcpp(Y[idx, , drop = FALSE], matrix(1, length(idx), 1), threads)
  list(S1 = as.vector(ss$S1), S2 = as.vector(ss$S2))
}

#' @keywords internal
perm_stratum_big <- function(ps, B, taus = NULL) {
  Knz <- length(ps$cols)
  big <- matrix(0, ps$ns, Knz * B)
  for (b in seq_len(B)) {
    tau <- if (is.null(taus)) sample.int(ps$ns) else taus[[b]]
    big[, (b - 1) * Knz + seq_len(Knz)] <- ps$Rfac * ps$omega[tau, , drop = FALSE]
  }
  big
}
