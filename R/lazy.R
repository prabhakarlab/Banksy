
# Lazy BANKSY operator for implicit PCA
#
# Computes PCA on the BANKSY matrix without materializing the full
# [2*n_genes x n_cells] product. The BanksyLazy operator implements
# implicit forward (A %*% x) and adjoint (t(A) %*% x) multiplications
# that irlba needs, operating directly on the sparse expression matrix
# and weight matrix.
#
# Ported from SeuratWrappers::RunBanksy lazy path.

# Core lazy PCA: framework-agnostic (no Seurat, no SE dependencies).
# Both Banksy::runBanksyLazyPCA and SeuratWrappers::RunBanksy call this.
#
# @param data_own dgCMatrix of expression data (genes x cells)
# @param knn_list list of kNN data.tables (single-matrix path)
# @param group_knn per-group kNN list (grouped path, NULL if ungrouped)
# @param group_idx list of integer vectors mapping group -> cell indices
# @param lambda numeric spatial weight in [0,1]
# @param npcs integer number of PCs to compute
# @param split_scale logical whether to scale per group
# @param scale_max numeric max absolute z-score for clipping
# @param pca_backend character "cpp" or "r"
# @param verbose logical
#
# @return list(embeddings, loadings, stdev, total_var)
.banksy_lazy_pca_core <- function(data_own,
                                   knn_list = NULL,
                                   group_knn = NULL,
                                   group_idx = NULL,
                                   lambda,
                                   npcs = 50L,
                                   split_scale = FALSE,
                                   scale_max = 10,
                                   pca_backend = c("cpp", "r"),
                                   verbose = TRUE) {

    n_genes <- nrow(data_own)
    n_cells <- ncol(data_own)
    gene_names <- rownames(data_own)
    cell_names <- colnames(data_own)
    lambdas <- getLambdas(lambda, n_harmonics = 1)
    has_groups <- !is.null(group_knn)

    # Validate npcs
    max_npcs <- min(2L * n_genes, n_cells) - 1L
    if (npcs >= max_npcs) {
        stop('npcs (', npcs, ') must be < min(2*n_genes, n_cells) = ',
             max_npcs + 1L)
    }
    if (min(2L * n_genes, n_cells) < 6L) {
        stop('lazy PCA requires at least 3 genes and 6 cells')
    }

    if (has_groups) {
        # Per-group decomposition
        if (verbose) message('Building per-group weight matrices')
        W_list <- lapply(seq_along(group_idx), function(gr)
            .buildWeightMatrix(group_knn[[gr]][[1]],
                               length(group_idx[[gr]])))
        rm(group_knn); gc(verbose = FALSE)

        if (verbose) message('Building per-group expression matrices')

        # Materialize per-group dgCMatrix
        gcm_list <- lapply(group_idx, function(cid)
            as(data_own[, cid], 'dgCMatrix'))

        own_result <- .lazy_own_scaling(NULL, split_scale, group_idx,
                                        groups = NULL, ugroups = NULL,
                                        n_genes, n_cells, scale_max, verbose,
                                        gcm_list = gcm_list)

        h0_result <- .lazy_h0_scaling_grouped(
            NULL, W_list, split_scale, group_idx,
            n_genes, n_cells, scale_max, verbose,
            gcm_list = gcm_list)

        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) length(own_result$excess@x) else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) length(h0_result$excess@x) else 0
            message('Clipping corrections: own=', n_own_clip,
                    ' H0=', n_h0_clip, ' entries')
        }

        # Extract W CSC slots
        w_csc <- lapply(W_list, function(w)
            if (is(w, 'dgCMatrix')) w else as(w, 'dgCMatrix'))
        w_slots <- list(
            i = lapply(w_csc, slot, 'i'),
            p = lapply(w_csc, slot, 'p'),
            x = lapply(w_csc, slot, 'x'),
            ncol = vapply(w_csc, ncol, integer(1)))
        rm(w_csc)

        gcm_slots <- list(
            i = lapply(gcm_list, slot, 'i'),
            p = lapply(gcm_list, slot, 'p'),
            x = lapply(gcm_list, slot, 'x'))
        rm(gcm_list, data_own); gc(verbose = FALSE)

        banksy_op <- structure(list(
            gcm_slots = gcm_slots,
            W_list = W_list, w_slots = w_slots,
            mu = list(own_result$mu, h0_result$mu),
            sd = list(own_result$sd, h0_result$sd),
            lam = lambdas,
            excess = list(own_result$excess, h0_result$excess),
            valid = list(own_result$valid, h0_result$valid),
            has_groups = TRUE,
            split_scale = split_scale,
            group_idx = group_idx,
            n_genes = n_genes, n_cells = n_cells
        ), class = 'BanksyLazy')
        rm(own_result, h0_result)
    } else {
        # Single-matrix path
        if (verbose) message('Building sparse weight matrix')
        W <- .buildWeightMatrix(knn_list[[1]], n_cells)
        rm(knn_list); gc(verbose = FALSE)

        own_result <- .lazy_own_scaling(data_own, split_scale, group_idx,
                                        groups = NULL, ugroups = NULL,
                                        n_genes, n_cells, scale_max, verbose)
        h0_result <- .lazy_h0_scaling(data_own, W, split_scale, group_idx,
                                      n_genes, n_cells, scale_max, verbose)

        if (verbose) {
            n_own_clip <- if (!is.null(own_result$excess)) length(own_result$excess@x) else 0
            n_h0_clip <- if (!is.null(h0_result$excess)) length(h0_result$excess@x) else 0
            message('Clipping corrections: own=', n_own_clip,
                    ' H0=', n_h0_clip, ' entries')
        }

        banksy_op <- structure(list(
            gcm = data_own, W = W,
            mu = list(own_result$mu, h0_result$mu),
            sd = list(own_result$sd, h0_result$sd),
            lam = lambdas,
            excess = list(own_result$excess, h0_result$excess),
            valid = list(own_result$valid, h0_result$valid),
            has_groups = FALSE,
            split_scale = split_scale,
            group_idx = group_idx,
            n_genes = n_genes, n_cells = n_cells
        ), class = 'BanksyLazy')
        rm(own_result, h0_result)
    }

    # Reclaim memory before iterative solve
    gc(verbose = FALSE)

    # SVD solve
    pca_backend <- match.arg(pca_backend, c("cpp", "r"))
    if (pca_backend == "cpp") {
        # irlba algorithm with C++ Q buffer management.
        # Ported from irlba::irlba R source. The implicit restart rotates
        # V via lanczos_extract (C++), avoiding the memory spike from
        # R's copy-on-modify on V %*% Bsvd$v.
        # Peak memory = Q (n_cells * work) + extraction buffer (n_cells * k).
        work <- as.integer(min(npcs + 7L, n_cells - 1L))
        maxit <- 1000L
        tol <- 1e-5
        svtol <- tol
        if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via ',
                             'C++ irlba (work=', work, ')')

        m <- 2L * n_genes
        n <- n_cells
        k <- npcs
        eps <- .Machine$double.eps
        eps23 <- eps^(2/3)
        sqrteps <- sqrt(eps)

        # C++ managed V buffer (n x work) -- never copied by R
        V_buf <- numeric(as.double(n) * work)
        W_mat <- matrix(0.0, m, work)
        F_vec <- numeric(n)
        B <- NULL

        set.seed(42)
        q1 <- rnorm(n); q1 <- q1 / lanczos_norm(q1)
        lanczos_set_col(V_buf, 0L, q1, n)
        rm(q1)

        t0 <- proc.time()["elapsed"]
        mprod <- 0L
        iter <- 0L
        Smax <- 1.0
        Smin <- NULL
        lastsv <- numeric(0)
        converged <- FALSE
        restart <- FALSE

        while (iter < maxit) {
            iter <- iter + 1L
            j <- if (iter == 1L && !restart) 1L else k + 1L

            # -- Lanczos bidiagonalization: steps j..work --
            # First step of each cycle
            VJ <- lanczos_get_col(V_buf, j - 1L, n)
            avj <- as.numeric(.banksy_lazy_mult(
                banksy_op, matrix(VJ, ncol = 1)))
            rm(VJ)
            W_mat[, j] <- avj
            mprod <- mprod + 1L

            if (iter > 1L && j > 1L)
                W_mat[, j] <- W_mat[, j] - as.numeric(
                    W_mat[, 1:(j-1), drop=FALSE] %*%
                    crossprod(W_mat[, 1:(j-1), drop=FALSE], W_mat[, j]))

            S <- sqrt(sum(W_mat[, j]^2))
            if (S < eps23) {
                W_mat[, j] <- rnorm(m)
                if (j > 1L) W_mat[, j] <- W_mat[, j] - as.numeric(
                    W_mat[, 1:(j-1), drop=FALSE] %*%
                    crossprod(W_mat[, 1:(j-1), drop=FALSE], W_mat[, j]))
                W_mat[, j] <- W_mat[, j] / sqrt(sum(W_mat[, j]^2))
                S <- 0
            } else {
                W_mat[, j] <- W_mat[, j] / S
            }

            # Inner Lanczos loop
            while (j <= work) {
                # Adjoint
                F_vec <- as.numeric(.banksy_lazy_mult(
                    banksy_op, matrix(W_mat[, j], ncol = 1), transpose = TRUE))
                mprod <- mprod + 1L
                F_vec <- F_vec - S * lanczos_get_col(V_buf, j - 1L, n)
                F_vec <- lanczos_reorth(V_buf, j, F_vec, n)

                if (j + 1L <= work) {
                    R <- lanczos_norm(F_vec)
                    if (R < eps23) {
                        F_vec <- rnorm(n)
                        F_vec <- lanczos_reorth(V_buf, j, F_vec, n)
                        lanczos_set_col(V_buf, j,
                            F_vec / lanczos_norm(F_vec), n)
                        R <- 0
                    } else {
                        lanczos_set_col(V_buf, j, F_vec / R, n)
                    }

                    # Build B
                    if (is.null(B)) {
                        B <- matrix(c(S, R), nrow = 1)
                    } else {
                        B <- rbind(cbind(B, 0),
                                   c(rep(0, ncol(B) - 1), S, R))
                    }

                    # Forward for next step
                    VJP1 <- lanczos_get_col(V_buf, j, n)
                    W_mat[, j + 1L] <- as.numeric(.banksy_lazy_mult(
                        banksy_op, matrix(VJP1, ncol = 1)))
                    rm(VJP1)
                    mprod <- mprod + 1L
                    W_mat[, j + 1L] <- W_mat[, j + 1L] - R * W_mat[, j]
                    if (j > 1L)
                        W_mat[, j + 1L] <- W_mat[, j + 1L] - as.numeric(
                            W_mat[, 1:j, drop=FALSE] %*%
                            crossprod(W_mat[, 1:j, drop=FALSE], W_mat[, j + 1L]))
                    S <- sqrt(sum(W_mat[, j + 1L]^2))
                    if (S < eps23) {
                        W_mat[, j + 1L] <- rnorm(m)
                        if (j > 0L) W_mat[, j + 1L] <- W_mat[, j + 1L] - as.numeric(
                            W_mat[, 1:j, drop=FALSE] %*%
                            crossprod(W_mat[, 1:j, drop=FALSE], W_mat[, j + 1L]))
                        W_mat[, j + 1L] <- W_mat[, j + 1L] / sqrt(sum(W_mat[, j + 1L]^2))
                        S <- 0
                    } else {
                        W_mat[, j + 1L] <- W_mat[, j + 1L] / S
                    }
                } else {
                    B <- rbind(B, c(rep(0, j - 1), S))
                }
                j <- j + 1L
            }

            # -- Convergence test --
            Bsz <- nrow(B)
            R_F <- lanczos_norm(F_vec)
            F_vec <- F_vec / R_F
            Bsvd <- svd(B)
            Smax <- max(Smax, Bsvd$d[1])

            R_vals <- R_F * Bsvd$u[Bsz, , drop = FALSE]
            ct_conv <- all(abs(R_vals[1:k]) < tol * Smax)
            sv_conv <- length(lastsv) >= k &&
                all(abs(Bsvd$d[1:k] - lastsv[1:k]) < svtol * Bsvd$d[1:k])
            lastsv <- Bsvd$d

            if (verbose && (iter <= 2L || iter %% 5 == 0 || ct_conv || sv_conv)) {
                rss <- tryCatch({
                    l <- readLines("/proc/self/status", warn = FALSE)
                    v <- grep("^VmRSS:", l, value = TRUE)
                    as.numeric(gsub("[^0-9]", "", v)) / 1024^2
                }, error = function(e) NA)
                elapsed <- proc.time()["elapsed"] - t0
                message(sprintf(
                    "  iter=%d  mprod=%d  sv[%d]=%.4e  RSS=%.1fGB  t=%.0fs",
                    iter, mprod, k, Bsvd$d[k], rss, elapsed))
            }

            if (ct_conv || sv_conv) {
                converged <- TRUE
                if (verbose) message(sprintf(
                    "  Converged: iter=%d, mprod=%d", iter, mprod))
                break
            }

            # -- Implicit restart --
            V_rot <- lanczos_extract(V_buf, Bsvd$v[, 1:k, drop = FALSE],
                                     n, Bsz, k)
            for (i in seq_len(k))
                lanczos_set_col(V_buf, i - 1L, V_rot[, i], n)
            rm(V_rot)
            lanczos_set_col(V_buf, k, F_vec, n)

            B <- cbind(diag(Bsvd$d[1:k], nrow = k),
                       R_F * Bsvd$u[Bsz, 1:k])

            W_mat[, 1:k] <- W_mat[, 1:Bsz, drop = FALSE] %*%
                Bsvd$u[, 1:k, drop = FALSE]

            restart <- TRUE
        }

        if (!converged)
            warning("C++ irlba did not converge; try increasing maxit or work")

        # Final extraction (C++ managed, no spike)
        Bsvd <- svd(B)
        V <- lanczos_extract(V_buf, Bsvd$v[, 1:npcs, drop = FALSE],
                             n, nrow(B), npcs)
        rm(V_buf); gc(verbose = FALSE)

        U <- W_mat[, 1:nrow(B), drop = FALSE] %*%
             Bsvd$u[, 1:npcs, drop = FALSE]
        rm(W_mat)

        d <- Bsvd$d[1:npcs]
        rm(Bsvd, B, F_vec)
        pca <- list(u = U, v = V, d = d)

    } else {
        if (verbose) message('Computing BANKSY PCA (', npcs, ' PCs) via irlba')

        .irlba_iter <- 0L
        .irlba_t0 <- proc.time()["elapsed"]
        .instrumented_mult <- function(A, x, transpose = FALSE) {
            .irlba_iter <<- .irlba_iter + 1L
            if (.irlba_iter %% 20 == 0) gc(verbose = FALSE)
            if (.irlba_iter %% 20 == 1) {
                rss <- tryCatch({
                    l <- readLines("/proc/self/status", warn = FALSE)
                    v <- grep("^VmRSS:", l, value = TRUE)
                    as.numeric(gsub("[^0-9]", "", v)) / 1024^2
                }, error = function(e) NA)
                elapsed <- proc.time()["elapsed"] - .irlba_t0
                message(sprintf("  irlba iter=%d  t=%.0fs  RSS=%.1fGB",
                                .irlba_iter, elapsed, rss))
            }
            .banksy_lazy_mult(A, x, transpose)
        }
        pca <- irlba::irlba(banksy_op, nv = npcs, mult = .instrumented_mult)
    }

    # Cell embeddings: V * D
    embeddings <- sweep(pca$v, 2, pca$d, `*`)
    rownames(embeddings) <- cell_names
    colnames(embeddings) <- paste0('PC', seq_len(npcs))

    # Feature loadings
    loadings <- pca$u
    feat_names <- c(gene_names, paste0(gene_names, '.m0'))
    rownames(loadings) <- feat_names
    colnames(loadings) <- paste0('PC', seq_len(npcs))

    # Percent variance
    total_var <- sum(pca$d^2)
    stdev <- pca$d / sqrt(max(1, n_cells - 1))

    list(embeddings = embeddings, loadings = loadings,
         stdev = stdev, total_var = total_var)
}

# Own expression scaling params + clipping excess
.lazy_own_scaling <- function(data_own, split_scale, group_idx, groups,
                              ugroups, n_genes, n_cells, scale_max, verbose,
                              gcm_list = NULL) {
    if (verbose) message('Computing scaling parameters for own expression')

    if (split_scale) {
        mu <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        sd <- matrix(0, nrow = n_genes, ncol = length(group_idx))
        for (gr in seq_along(group_idx)) {
            grp <- if (!is.null(gcm_list)) gcm_list[[gr]]
                   else data_own[, group_idx[[gr]], drop = FALSE]
            n_c <- as.double(ncol(grp))
            mu[, gr] <- Matrix::rowMeans(grp)
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (Matrix::rowMeans(grp^2) - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
    } else {
        if (!is.null(gcm_list)) {
            # Aggregate global stats from per-group data
            row_sums <- numeric(n_genes)
            row_sq_sums <- numeric(n_genes)
            for (gr in seq_along(gcm_list)) {
                row_sums <- row_sums + Matrix::rowSums(gcm_list[[gr]])
                row_sq_sums <- row_sq_sums + Matrix::rowSums(gcm_list[[gr]]^2)
            }
            n_c <- as.double(n_cells)
            mu <- row_sums / n_c
            sd <- sqrt(pmax(n_c / (n_c - 1) * (row_sq_sums / n_c - mu^2), 0))
        } else {
            n_c <- as.double(n_cells)
            if (inherits(data_own, 'sparseMatrix')) {
                mu <- Matrix::rowMeans(data_own)
                sd <- sqrt(pmax(
                    n_c / (n_c - 1) * (Matrix::rowMeans(data_own^2) - mu^2), 0
                ))
            } else {
                mu <- rowMeans(data_own)
                sd <- sqrt(pmax(
                    n_c / (n_c - 1) * (rowMeans(data_own^2) - mu^2), 0
                ))
            }
        }
        sd[sd == 0] <- 1
        valid <- NULL
    }

    # Compute clipping excess via C++ (no R temporaries)
    if (verbose) message('Computing clipping excess for own expression')
    if (!is.null(gcm_list)) {
        mu_mat <- if (is.matrix(mu)) mu else matrix(0, 0, 0)
        sd_mat <- if (is.matrix(sd)) sd else matrix(0, 0, 0)
        mu_vec <- if (!is.matrix(mu)) mu else numeric(0)
        sd_vec <- if (!is.matrix(sd)) sd else numeric(0)
        exc <- own_excess_cpp(
            lapply(gcm_list, slot, 'i'), lapply(gcm_list, slot, 'p'),
            lapply(gcm_list, slot, 'x'), group_idx,
            split_scale, mu_mat, sd_mat, mu_vec, sd_vec,
            valid, scale_max, n_genes, n_cells)
        excess <- if (length(exc$i) > 0L)
            Matrix::sparseMatrix(i = exc$i, j = exc$j, x = exc$x,
                                 dims = c(n_genes, n_cells))
        else NULL
    } else {
        excess <- .lazy_own_excess(data_own, mu, sd, valid, split_scale,
                                   group_idx, groups, ugroups,
                                   n_genes, n_cells, scale_max)
    }

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

.lazy_own_excess <- function(data_own, mu, sd, valid, split_scale, group_idx,
                             groups, ugroups, n_genes, n_cells, scale_max,
                             gcm_list = NULL) {
    if (!is.null(gcm_list)) {
        # Per-group excess from materialized dgCMatrix list
        exc_i <- integer(0); exc_j <- integer(0); exc_x <- numeric(0)
        for (gr in seq_along(gcm_list)) {
            g <- gcm_list[[gr]]
            cid <- group_idx[[gr]]
            gi <- g@i + 1L
            gj <- rep(cid, diff(g@p))
            if (split_scale) {
                g_mu <- mu[cbind(gi, gr)]
                g_sd <- sd[cbind(gi, gr)]
            } else {
                g_mu <- mu[gi]
                g_sd <- sd[gi]
            }
            exceed <- g@x > g_mu + scale_max * g_sd
            if (split_scale) exceed <- exceed & valid[gi]
            if (any(exceed)) {
                z_vals <- (g@x[exceed] - g_mu[exceed]) / g_sd[exceed]
                exc_i <- c(exc_i, gi[exceed])
                exc_j <- c(exc_j, gj[exceed])
                exc_x <- c(exc_x, z_vals - scale_max)
            }
        }
        if (length(exc_i) > 0)
            return(Matrix::sparseMatrix(i = exc_i, j = exc_j, x = exc_x,
                                        dims = c(n_genes, n_cells)))
        else
            return(NULL)
    }
    if (inherits(data_own, 'sparseMatrix')) {
        own_i <- data_own@i + 1L
        own_j <- rep(seq_len(n_cells), diff(data_own@p))
        if (split_scale) {
            own_gr <- match(groups[own_j], ugroups)
            own_mu <- mu[cbind(own_i, own_gr)]
            own_sd <- sd[cbind(own_i, own_gr)]
        } else {
            own_mu <- mu[own_i]
            own_sd <- sd[own_i]
        }
        exceed <- data_own@x > own_mu + scale_max * own_sd
        if (split_scale) exceed <- exceed & valid[own_i]
        if (any(exceed)) {
            ei <- own_i[exceed]
            z_vals <- (data_own@x[exceed] - own_mu[exceed]) / own_sd[exceed]
            Matrix::sparseMatrix(
                i = ei, j = own_j[exceed],
                x = z_vals - scale_max,
                dims = c(n_genes, n_cells))
        } else {
            NULL
        }
    } else {
        if (split_scale) {
            excess <- Matrix::sparseMatrix(
                i = integer(0), j = integer(0),
                dims = c(n_genes, n_cells))
            for (gr in seq_along(group_idx)) {
                cid <- group_idx[[gr]]
                z_own <- (data_own[, cid, drop = FALSE] - mu[, gr]) / sd[, gr]
                z_own[!valid, ] <- 0
                exceed_mask <- z_own > scale_max
                if (any(exceed_mask)) {
                    excess[, cid] <- Matrix::Matrix(
                        (z_own - scale_max) * exceed_mask, sparse = TRUE)
                }
            }
            if (length(excess@x) == 0) NULL else excess
        } else {
            z_own <- (data_own - mu) / sd
            exceed_mask <- z_own > scale_max
            if (any(exceed_mask)) {
                Matrix::Matrix((z_own - scale_max) * exceed_mask, sparse = TRUE)
            } else {
                NULL
            }
        }
    }
}

# H0 scaling params + two-pass clipping excess
.lazy_h0_scaling <- function(data_own, W, split_scale, group_idx,
                             n_genes, n_cells, scale_max, verbose) {
    if (verbose) message('Computing scaling params and clipping for H0')
    chunk_sz <- 100L

    # Pass 1: compute mu, ss, and per-gene max of H0 = data_own %*% W
    if (split_scale) {
        n_groups <- length(group_idx)
        mu <- matrix(0, nrow = n_genes, ncol = n_groups)
        ss <- matrix(0, nrow = n_genes, ncol = n_groups)
        max_h0 <- matrix(-Inf, nrow = n_genes, ncol = n_groups)
        for (ch_start in seq(1L, n_genes, by = chunk_sz)) {
            ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
            ri <- ch_start:ch_end
            chunk <- as.matrix(data_own[ri, , drop = FALSE] %*% W)
            for (gr in seq_along(group_idx)) {
                cid <- group_idx[[gr]]
                chunk_group <- chunk[, cid, drop = FALSE]
                mu[ri, gr] <- rowMeans(chunk_group)
                ss[ri, gr] <- rowSums(chunk_group * chunk_group)
                max_h0[ri, gr] <- apply(chunk_group, 1, max)
            }
        }
        sd <- matrix(0, nrow = n_genes, ncol = n_groups)
        for (gr in seq_along(group_idx)) {
            n_c <- as.double(length(group_idx[[gr]]))
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (ss[, gr] / n_c - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
        thresh <- mu + scale_max * sd
    } else {
        mu <- numeric(n_genes)
        ss <- numeric(n_genes)
        max_h0 <- rep(-Inf, n_genes)
        for (ch_start in seq(1L, n_genes, by = chunk_sz)) {
            ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
            ri <- ch_start:ch_end
            chunk <- as.matrix(data_own[ri, , drop = FALSE] %*% W)
            mu[ri] <- rowMeans(chunk)
            ss[ri] <- rowSums(chunk * chunk)
            max_h0[ri] <- apply(chunk, 1, max)
        }
        n_c <- as.double(n_cells)
        sd <- sqrt(pmax(n_c / (n_c - 1) * (ss / n_c - mu * mu), 0))
        sd[sd == 0] <- 1
        thresh <- mu + scale_max * sd
        valid <- NULL
    }

    # Identify genes that need clipping
    if (split_scale) {
        clip_genes <- which(valid & rowSums(max_h0 > thresh) > 0)
    } else {
        clip_genes <- which(max_h0 > thresh)
    }
    if (verbose) message('H0 genes requiring clipping: ', length(clip_genes),
                         ' / ', n_genes)

    # Pass 2: compute excess for clipped genes only
    excess <- .lazy_h0_excess(data_own, W, mu, sd, split_scale, group_idx,
                              clip_genes, chunk_sz, n_genes, n_cells, scale_max)

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

.lazy_h0_excess <- function(data_own, W, mu, sd, split_scale, group_idx,
                            clip_genes, chunk_sz, n_genes, n_cells, scale_max) {
    if (length(clip_genes) == 0) return(NULL)

    exc_cap <- max(1024L, length(clip_genes) * 10L)
    exc_i <- integer(exc_cap)
    exc_j <- integer(exc_cap)
    exc_x <- numeric(exc_cap)
    exc_n <- 0L

    clip_chunks <- unique((clip_genes - 1L) %/% chunk_sz)
    for (ch_idx in clip_chunks) {
        ch_start <- ch_idx * chunk_sz + 1L
        ch_end <- min(ch_start + chunk_sz - 1L, n_genes)
        ri <- ch_start:ch_end
        ri_clip <- ri[ri %in% clip_genes]
        chunk <- as.matrix(data_own[ri_clip, , drop = FALSE] %*% W)
        scan_idx <- if (split_scale) group_idx else list(seq_len(n_cells))
        for (gr in seq_along(scan_idx)) {
            cid <- scan_idx[[gr]]
            if (split_scale) {
                z_chunk <- (chunk[, cid, drop = FALSE] -
                    mu[ri_clip, gr]) / sd[ri_clip, gr]
            } else {
                z_chunk <- (chunk - mu[ri_clip]) / sd[ri_clip]
            }
            wh <- which(z_chunk > scale_max, arr.ind = TRUE)
            if (nrow(wh) > 0) {
                new_n <- nrow(wh)
                while (exc_n + new_n > length(exc_i)) {
                    exc_cap <- exc_cap * 2L
                    length(exc_i) <- exc_cap
                    length(exc_j) <- exc_cap
                    length(exc_x) <- exc_cap
                }
                idx <- seq(exc_n + 1L, exc_n + new_n)
                exc_i[idx] <- ri_clip[wh[, 1]]
                exc_j[idx] <- cid[wh[, 2]]
                exc_x[idx] <- z_chunk[wh] - scale_max
                exc_n <- exc_n + new_n
            }
        }
    }

    if (exc_n == 0L) return(NULL)
    Matrix::sparseMatrix(
        i = exc_i[1:exc_n], j = exc_j[1:exc_n],
        x = exc_x[1:exc_n], dims = c(n_genes, n_cells))
}

# Per-group H0 scaling via C++ column-wise sweep (no intermediate allocation)
.lazy_h0_scaling_grouped <- function(data_own, W_list, split_scale, group_idx,
                                     n_genes, n_cells, scale_max, verbose,
                                     gcm_list = NULL) {
    if (verbose) message('Computing scaling params and clipping for H0 (per-group)')
    n_groups <- length(group_idx)

    # Extract gcm CSC slots (always per-group)
    gcm_i <- lapply(gcm_list, slot, 'i')
    gcm_p <- lapply(gcm_list, slot, 'p')
    gcm_x <- lapply(gcm_list, slot, 'x')

    # Extract W CSC slots
    w_csc <- lapply(W_list, function(w)
        if (is(w, 'dgCMatrix')) w else as(w, 'dgCMatrix'))
    w_i <- lapply(w_csc, slot, 'i')
    w_p <- lapply(w_csc, slot, 'p')
    w_x <- lapply(w_csc, slot, 'x')
    rm(w_csc)

    if (verbose) {
        for (gr in seq_along(group_idx))
            message('  Group ', gr, ': ', length(group_idx[[gr]]), ' cells')
    }

    # Pass 1: per-group stats via C++
    stats <- h0_group_stats_cpp(
        gcm_i, gcm_p, gcm_x, nrow(gcm_list[[1]]),
        w_i, w_p, w_x, group_idx)

    mu_g <- stats$mu
    ss_g <- stats$ss
    max_g <- stats$max_val
    rm(stats)

    # Derive mu, sd, valid
    if (split_scale) {
        mu <- mu_g
        sd <- matrix(0, nrow = n_genes, ncol = n_groups)
        for (gr in seq_along(group_idx)) {
            n_c <- as.double(length(group_idx[[gr]]))
            sd[, gr] <- sqrt(pmax(
                n_c / (n_c - 1) * (ss_g[, gr] / n_c - mu[, gr]^2), 0
            ))
        }
        valid <- rowSums(sd == 0) == 0
        sd[sd == 0] <- 1
    } else {
        n_g_vec <- vapply(group_idx, length, integer(1))
        mu <- rowSums(sweep(mu_g, 2, n_g_vec, `*`)) / n_cells
        ss <- rowSums(ss_g)
        n_c <- as.double(n_cells)
        sd <- sqrt(pmax(n_c / (n_c - 1) * (ss / n_c - mu * mu), 0))
        sd[sd == 0] <- 1
        valid <- NULL
    }

    # Identify clip genes
    if (split_scale) {
        thresh <- mu + scale_max * sd
        clip_genes <- which(valid & rowSums(max_g > thresh) > 0)
    } else {
        thresh <- mu + scale_max * sd
        max_h0 <- apply(max_g, 1, max)
        clip_genes <- which(max_h0 > thresh)
    }
    rm(mu_g, ss_g, max_g)
    if (verbose) message('H0 genes requiring clipping: ', length(clip_genes),
                         ' / ', n_genes)

    # Pass 2: excess via C++ (only for clip genes, no dense intermediates)
    if (length(clip_genes) == 0) {
        excess <- NULL
    } else {
        mu_mat <- if (is.matrix(mu)) mu else matrix(mu, ncol = 1)
        sd_mat <- if (is.matrix(sd)) sd else matrix(sd, ncol = 1)
        exc <- h0_group_excess_cpp(
            gcm_i, gcm_p, gcm_x, nrow(gcm_list[[1]]),
            w_i, w_p, w_x, group_idx,
            as.integer(clip_genes), split_scale, mu_mat, sd_mat, scale_max)
        excess <- if (length(exc$i) > 0L)
            Matrix::sparseMatrix(i = exc$i, j = exc$j, x = exc$x,
                                 dims = c(n_genes, n_cells))
        else NULL
    }

    list(mu = mu, sd = sd, valid = valid, excess = excess)
}

# Lazy operator multiply (forward + adjoint)
#' S3 dim method for BanksyLazy operator
#' @param x A BanksyLazy object
#' @return Integer vector of length 2: c(2*n_genes, n_cells)
#' @export
dim.BanksyLazy <- function(x) c(x$n_genes * 2L, x$n_cells)

.as_base <- function(x) {
    if (inherits(x, 'Matrix')) x <- as.matrix(x)
    if (!is.matrix(x)) x <- as.matrix(x)
    storage.mode(x) <- 'double'
    x
}

# --- Lazy operator dispatch ---

.banksy_lazy_mult <- function(A, x, transpose = FALSE) {
    if (inherits(x, 'BanksyLazy')) {
        tmp <- A; A <- x; x <- tmp; transpose <- TRUE
    }
    x <- .as_base(x)
    k <- ncol(x)
    if (!transpose) .banksy_forward(A, x, k)
    else .banksy_adjoint(A, x, k)
}

.banksy_forward <- function(A, x, k) {
    if (isTRUE(A$has_groups)) {
        ws <- A$w_slots
        gs_i <- A$gcm_slots$i
        gs_p <- A$gcm_slots$p
        gs_x <- A$gcm_slots$x

        if (A$split_scale) {
            result <- banksy_forward_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx, x,
                TRUE, A$mu[[1]], A$sd[[1]], A$mu[[2]], A$sd[[2]],
                A$lam, A$valid[[1]], A$valid[[2]])

            if (!is.null(A$excess[[1]]) || !is.null(A$excess[[2]])) {
                ng <- A$n_genes
                own <- result[1:ng, , drop = FALSE]
                h0  <- result[(ng+1):(2*ng), , drop = FALSE]
                if (!is.null(A$excess[[1]]))
                    own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
                if (!is.null(A$excess[[2]]))
                    h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
                result <- rbind(own, h0)
            }
        } else {
            result <- banksy_forward_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx, x,
                FALSE,
                matrix(0, 0, 0), matrix(0, 0, 0),
                matrix(0, 0, 0), matrix(0, 0, 0),
                A$lam, NULL, NULL)

            ng <- A$n_genes
            own <- result[1:ng, , drop = FALSE]
            h0  <- result[(ng+1):(2*ng), , drop = FALSE]
            cs <- colSums(x)
            own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
            h0  <- A$lam[2] * (h0  - outer(A$mu[[2]], cs)) / A$sd[[2]]
            if (!is.null(A$excess[[1]]))
                own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
            if (!is.null(A$excess[[2]]))
                h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
            result <- rbind(own, h0)
        }
        result
    } else {
        cs <- colSums(x)
        own <- .as_base(A$gcm %*% x)
        own <- A$lam[1] * (own - outer(A$mu[[1]], cs)) / A$sd[[1]]
        if (!is.null(A$excess[[1]]))
            own <- own - A$lam[1] * .as_base(A$excess[[1]] %*% x)
        Wx <- .as_base(A$W %*% x)
        h0 <- .as_base(A$gcm %*% Wx)
        h0 <- A$lam[2] * (h0 - outer(A$mu[[2]], cs)) / A$sd[[2]]
        if (!is.null(A$excess[[2]]))
            h0 <- h0 - A$lam[2] * .as_base(A$excess[[2]] %*% x)
        rbind(own, h0)
    }
}

.banksy_adjoint <- function(A, x, k) {
    ng <- A$n_genes
    xo <- x[1:ng, , drop = FALSE]
    xh <- x[(ng+1):(2*ng), , drop = FALSE]

    if (isTRUE(A$has_groups)) {
        ws <- A$w_slots
        gs_i <- A$gcm_slots$i
        gs_p <- A$gcm_slots$p
        gs_x <- A$gcm_slots$x

        if (A$split_scale) {
            if (!is.null(A$valid[[1]])) xo[!A$valid[[1]], ] <- 0
            if (!is.null(A$valid[[2]])) xh[!A$valid[[2]], ] <- 0

            r <- banksy_adjoint_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx,
                xo, xh,
                numeric(k), numeric(k),
                TRUE, A$mu[[1]], A$sd[[1]], A$mu[[2]], A$sd[[2]],
                A$lam)

            if (!is.null(A$excess[[1]]))
                r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
            if (!is.null(A$excess[[2]]))
                r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        } else {
            xo_s <- xo / A$sd[[1]]
            xh_s <- xh / A$sd[[2]]
            adj_o <- colSums(A$mu[[1]] * xo_s)
            adj_h <- colSums(A$mu[[2]] * xh_s)

            r <- banksy_adjoint_cpp(
                gs_i, gs_p, gs_x, A$n_genes, A$n_cells,
                ws$i, ws$p, ws$x, ws$ncol,
                A$group_idx,
                xo_s, xh_s, adj_o, adj_h,
                FALSE,
                matrix(0, 0, 0), matrix(0, 0, 0),
                matrix(0, 0, 0), matrix(0, 0, 0),
                A$lam)

            if (!is.null(A$excess[[1]]))
                r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
            if (!is.null(A$excess[[2]]))
                r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        }
        r
    } else {
        xo_s <- xo / A$sd[[1]]
        xh_s <- xh / A$sd[[2]]
        adj_o <- colSums(A$mu[[1]] * xo_s)
        r <- A$lam[1] * (.as_base(crossprod(A$gcm, xo_s)) -
             matrix(adj_o, A$n_cells, k, byrow = TRUE))
        if (!is.null(A$excess[[1]]))
            r <- r - A$lam[1] * .as_base(crossprod(A$excess[[1]], xo))
        adj_h <- colSums(A$mu[[2]] * xh_s)
        ht <- .as_base(crossprod(A$gcm, xh_s))
        ht <- .as_base(crossprod(A$W, ht))
        ht <- ht - matrix(adj_h, A$n_cells, k, byrow = TRUE)
        r <- r + A$lam[2] * ht
        if (!is.null(A$excess[[2]]))
            r <- r - A$lam[2] * .as_base(crossprod(A$excess[[2]], xh))
        r
    }
}
