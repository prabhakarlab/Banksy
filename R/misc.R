#' Simulate an unrealistic spatial omics dataset.
#'
#' @param n_cells An integer scalar specifying the approximate number of cells.
#' @param n_genes An integer scalar specifying the number of genes.
#' @param n_rings An integer scalar specifying the number of spatial rings.
#' @param rate A numeric scalar specifying the Poisson rate parameter for
#'  simulating counts.
#'
#'
#' @importFrom stats rpois runif
#' @importFrom SpatialExperiment SpatialExperiment
#'
#' @return A SpatialExperiment object.
#'
#' @export
#'
#' @examples
#' set.seed(2023)
#' rings <- simulateDataset(n_cells = 5e3, n_genes = 50, n_rings = 8)
#' rings
#' table(rings$cluster)
#' ggspavis::plotSpots(rings, annotate = 'cluster', size = 2)
#'
simulateDataset <- function(n_cells = 300,
                            n_genes = 30,
                            n_rings = 3,
                            rate = 10) {
    # Generate spatial locations
    n_cells <- ceiling(n_cells / (0.25 * pi))
    x <- runif(n_cells, min = -1, max = 1)
    y <- runif(n_cells, min = -1, max = 1)
    r <- sqrt(x^2 + y^2)
    points <- data.frame(x, y, r)
    points <- points[order(r), ]
    r0 <- 1 / sqrt(n_rings)
    radii_thres <- c(0, sqrt(seq(n_rings)) * r0, Inf)
    points$group <-
        factor(cut(points$r, radii_thres), labels = seq(n_rings + 1))
    points <- points[points$group != n_rings + 1, ]
    points$group <- droplevels(points$group)
    n_cells <- nrow(points)

    # Generate expression
    counts <- matrix(rpois(n_genes * n_cells, lambda = rate),
        nrow = n_genes,
        ncol = n_cells
    )
    n_pring <- floor(n_genes / n_rings)
    gene_idx <-
        split(seq(n_pring * n_rings), rep(seq(n_rings), each = n_pring))
    for (i in seq(n_rings)) {
        curr_gid <- gene_idx[[i]]
        curr_cid <- which(points$group == i)
        counts[curr_gid, curr_cid] <-
            counts[curr_gid, curr_cid] +
            rpois(length(curr_gid) * length(curr_cid), lambda = rate / 2)
    }
    locs <- as.matrix(points[, c("x", "y")])
    se <- SpatialExperiment(
        assays = list(counts = counts),
        spatialCoords = locs,
        colData = data.frame(cluster = points$group, in_tissue = TRUE)
    )
    rownames(se) <- paste0("gene_", seq(nrow(se)))
    colnames(se) <- paste0("cell_", seq(ncol(se)))
    se
}
