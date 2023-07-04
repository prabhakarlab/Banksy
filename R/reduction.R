
#' Run PCA
#'
#' @param bank BanksyObject
#' @param lambda (numeric vector) spatial weighting parameter
#' @param M (numeric vector) run PCA using up to the m-th azimuthal fourier harmonic (default: 1) 
#' @param npcs (numeric) number of principal components to compute
#' @param verbose (logical) print messages
#'
#' @importFrom irlba prcomp_irlba
#'
#' @return BanksyObject with PCA
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.2)
#' 
RunBanksyPCA <- function(bank, lambda = 0.2, M = 1, npcs = 20, verbose = TRUE) {
    
    for (m in M) {
        for (lam in lambda) {
            if (verbose) message('Running PCA for M=', m, ' lambda=', lam)
            
            x <- getBanksyMatrix(bank, lambda = lam, M = m)$expr
            cell.names <- colnames(x)
            
            pca <- prcomp_irlba(t(x), n = npcs)
            pca$var <- pca$sdev ^ 2 / sum(pca$sdev ^ 2) * 100
            rownames(pca$x) <- cell.names
            pca.name <- paste0('pca_M', m, '_lam', lam)
            bank@reduction[[pca.name]] <- pca
            
        }
    }
  return(bank)
}

#' Plots scree plot for PCA
#'
#' @param bank BanksyObject
#' @param lambda (numeric) spatial weighting parameter
#' @param M (numeric) compute up to the k-th azimuthal fourier harmonic (default: 1) 
#'
#' @importFrom ggplot2 ggplot aes geom_point theme element_blank element_line
#'
#' @return Scree plot
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.2)
#' scree <- plotScree(bank, lambda = 0.2)
#' scree
#' 
plotScree <- function(bank, lambda, M = 1) {

  pca.name <- paste0('pca_M', M, '_lam', lambda)
  found <- pca.name %in% names(bank@reduction)
  if (!found) stop('Run PCA with M=', M, ' lambda=', lambda)
  x <- bank@reduction[[pca.name]]$var
  data <- data.frame(PCs = seq_len(length(x)), Variance = x)
  PCs <- Variance <- NULL
  plot <- ggplot(data, aes(x = PCs, y = Variance)) + geom_point() +
    theme_blank(0,0,'none')

  return(plot)
}

#' Run UMAP
#'
#' @param bank BanksyObject
#' @param lambda (numeric vector) spatial weighting parameter
#' @param M (numeric) compute up to the k-th azimuthal fourier harmonic (default: 1) 
#' @param ncomponents (numeric) number of umap components to compute
#' @param pca (logical) run UMAP on PCs (TRUE) or BANKSY matrix (FALSE)
#' @param npcs (numeric) number of principal components to use for umap
#' @param nneighbors (numeric) number of neighbors to use for umap
#' @param spread (numeric) effective scale of embedded points
#' @param mindist (numeric) effective min. dist. between embedded points
#' @param nepochs (numeric) number of epochs to run umap optimization
#' @param verbose (logical) print messages
#' @param ... parameters to pass to uwot::umap
#'
#' @importFrom uwot umap
#'
#' @return BanksyObject with dimensionality reduction
#'
#' @export
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' bank <- NormalizeBanksy(bank)
#' bank <- ScaleBanksy(bank)
#' bank <- ComputeBanksy(bank)
#' bank <- RunBanksyPCA(bank, lambda = 0.2)
#' bank <- RunBanksyUMAP(bank, lambda = 0.2)
#' 
RunBanksyUMAP <- function(bank, lambda = 0.2, M = 1, ncomponents = 2, pca = TRUE, npcs = 20,
                    nneighbors = 30, spread = 3, mindist = 0.1, nepochs = 300,
                    verbose = TRUE, ...) {
    
    for (m in M) {
        for (lam in lambda) {
            if (!pca) {
                if (verbose) message('Computing UMAP on Banksy matrix')
                x <- t(getBanksyMatrix(bank, lambda = lam, M = m)$expr)
                
            } else {
                if (npcs < 2)
                    stop('Require more than 2 PCs')
                pca.name <- paste0('pca_M', m, '_lam', lam)
                found <- pca.name %in% names(bank@reduction)
                if (!found)
                    stop('Compute PCA with M=', m, ' lambda=', lam, ' before running UMAP')
                x <- bank@reduction[[pca.name]]$x
                nc <- min(ncol(x), npcs)
                if (nc < npcs)
                    warning('Only ', nc, ' PCs available')
                message('Computing UMAP with ', nc, ' PCs')
                x <- x[, seq_len(nc)]
            }
            
            if (verbose) message('Running UMAP for M=', m, ' lambda=', lam)
            umap <-
                umap(
                    x,
                    n_neighbors = nneighbors,
                    n_components = ncomponents,
                    spread = spread,
                    min_dist = mindist,
                    n_epochs = nepochs,
                    ...
                )
            rownames(umap) <- rownames(x)
            colnames(umap) <- paste0('UMAP_', seq_len(ncol(umap)))
            umap.name <- paste0('umap_M', m, '_lam', lam)
            bank@reduction[[umap.name]] <- umap
            
        }
    }

  return(bank)
}
