# DimReduction

#' Run PCA
#'
#' @param bank BanksyObject
#' @param lambda numeric vector - weighting parameter
#' @param npcs number of principal components to compute
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
#' bank <- RunPCA(bank, lambda = 0.2)
#' 
RunPCA <- function(bank, lambda, npcs = 20) {

  for (lam in lambda) {

    message('Running PCA for lambda=', lam)

    x <- getBanksyMatrix(bank, lambda = lam)$expr
    cell.names <- colnames(x)

    pca <- prcomp_irlba(t(x), n = npcs)
    pca$var <- pca$sdev^2 / sum(pca$sdev^2) * 100
    rownames(pca$x) <- cell.names
    pca.name <- paste0('pca_', lam)
    bank@reduction[[pca.name]] <- pca

  }
  return(bank)
}

#' Plots scree plot for PCA
#'
#' @param bank BanksyObject
#' @param lambda numeric
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
#' bank <- RunPCA(bank, lambda = 0.2)
#' scree <- plotScree(bank, lambda = 0.2)
#' scree
#' 
plotScree <- function(bank, lambda) {

  pca.name <- paste0('pca_', lambda)
  found <- pca.name %in% names(bank@reduction)
  if (!found) stop('Run PCA with lambda=', lambda)
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
#' @param lambda numeric vector - weighting parameter
#' @param ncomponents number of umap components to compute
#' @param pca if TRUE, run UMAP on pca, if not runs on Banksy matrix
#' @param npcs number of principal components to use for umap
#' @param nneighbors number of neighbors to use for umap
#' @param spread effective scale of embedded points
#' @param mindist effective min. dist. between embedded points
#' @param nepochs number of epochs to run umap optimization
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
#' bank <- RunPCA(bank, lambda = 0.2)
#' bank <- RunUMAP(bank, lambda = 0.2)
#' 
RunUMAP <- function(bank, lambda, ncomponents = 2, pca = TRUE, npcs = 20,
                    nneighbors = 30, spread = 3, mindist = 0.3, nepochs = 300,
                    ...) {

  for (lam in lambda) {

    if (!pca) {

      message('Computing UMAP on Banksy matrix')
      x <- t(getBanksyMatrix(bank, lambda = lambda)$expr)

    } else {

      if (npcs < 2) stop('Require more than 2 PCs')
      pca.name <- paste0('pca_', lam)
      found <- pca.name %in% names(bank@reduction)
      if (!found) stop('Compute PCA with lambda ', lam, ' before running UMAP')
      x <- bank@reduction[[pca.name]]$x
      nc <- min(ncol(x), npcs)
      if (nc < npcs) warning('Only ', nc, ' PCs available')
      message('Computing UMAP with ', nc, ' PCs')
      x <- x[, seq_len(nc)]
    }

    message('Running UMAP for lambda=', lam)
    umap <- umap(x, n_neighbors = nneighbors, n_components = ncomponents,
                 spread = spread, min_dist = mindist, n_epochs = nepochs, ...)
    rownames(umap) <- rownames(x)
    colnames(umap) <- paste0('UMAP_', seq_len(ncol(umap)))
    umap.name <- paste0('umap_', lam)
    bank@reduction[[umap.name]] <- umap

  }

  return(bank)
}
