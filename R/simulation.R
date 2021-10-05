
# Simulate dataset for testing
#' @importFrom stats rnorm
simulateDataset <- function(n_clusters = 5,
                            n_cells = 100,
                            n_genes = 100,
                            seed = 1234) {
  set.seed(seed)

  ## Get count for each cluster
  count_clusters <- vector(mode = 'numeric', length = n_clusters)
  for (i in seq_len(n_clusters-1)) {
    count_clusters[i] <- floor(n_cells / n_clusters)
  }
  count_clusters[n_clusters] <- n_cells - sum(count_clusters)
  count_clusters <- count_clusters[sample(n_clusters)]
  memberships <- rep(seq_len(n_clusters), count_clusters)

  gcm <- matrix(data = rnorm(n_genes * n_cells, mean = 10, sd = 2),
                nrow = n_genes)
  for (i in seq_len(n_clusters)) {
    gene_idx <- sample(n_genes, floor(n_genes / n_clusters / 2))
    cell_idx <- which(memberships == i)
    gcm[gene_idx, cell_idx] <- gcm[gene_idx, cell_idx] +
      rnorm(length(gene_idx)*length(cell_idx), mean = 5)
  }

  ## Set hyperpriors
  mu_x0 <- rnorm(n_clusters, mean = 1000, sd = 250)
  mu_y0 <- rnorm(n_clusters, mean = 1000, sd = 250)
  sigma <- rnorm(n_clusters, mean = 100, sd = 10)

  ## Generate Locs
  locs <- Map(function(mu_x, mu_y, sigma, num_cells){
    data.frame(sdimx = rnorm(num_cells, mean = mu_x, sd = sigma),
               sdimy = rnorm(num_cells, mean = mu_y, sd = sigma))},
    mu_x0, mu_y0, sigma, count_clusters)
  locs <- do.call(rbind.data.frame, locs)

  permute <- sample(n_cells)
  gcm <- gcm[, permute]
  locs <- locs[permute, ]
  memberships <- memberships[permute]

  rownames(gcm) <- paste0('gene_', seq_len(n_genes))
  colnames(gcm) <- paste0('cell_', seq_len(n_cells))
  rownames(locs) <- colnames(gcm)

  gcm <- as.matrix(gcm)
  gcm[gcm < 0] <- 0

  return(list(gcm, locs, memberships))
}

