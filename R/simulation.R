
#' Simulate a dataset consisting of expression and location information 
#' with concentric circular structure
#' 
#' @param n_buds Number of rings
#' @param n_genes Number of genes
#' @param n_deg Number of de genes per bud
#' @param de_shift Mean shift for de genes per bud
#' @param de_dropout Fraction of cells where gene is de per bud
#' @param n_cells Number of cells per bud
#' @param n_circles Number of circles per bud
#' @param n_thickness Thickness of each circle per bud
#' @param radius_init Initial radius of each bud
#' @param radius_increment Radius increment per bud
#' @param x_stretch Amount of stretching in the x-direction per bud
#' @param y_stretch Amount of stretching in the y-direction per bud
#' @param offset_x Offset to shift x center of bud
#' @param offset_y Offset to shift y center of bud
#' @param skew Rotation of bud
#' @param seed Seed for sampling
#' 
#' @importFrom stats rnorm
#'  
#' @return List of expression, locations, and metadata
#'  
#' @examples 
#' data <- simulateDataset() 
#'  
#' @export
simulateDataset <- function(
    n_buds = 3, n_genes = 100,
    n_deg = c(10,5,7), de_shift = c(3,4,2), de_dropout = c(0.8,0.9,0.7),
    n_cells = c(250, 150, 100), n_circles = c(3,3,3),
    n_thickness = list(c(3,1,3), c(2,1,4), c(2,2,2)),
    radius_init = c(5,3,3), radius_increment = c(4,3,3),
    x_stretch = c(1.1,1,1), y_stretch = c(0.9,1.1,1),
    offset_x = c(0,50,20), offset_y = c(0,40,35),
    skew = c(0,0,0.5), seed = 1000) {
    
    locs <- generateLocs(n_buds = n_buds, 
                         n_circles = n_circles,
                         n_cells = n_cells,
                         n_thickness = n_thickness,
                         radius_init = radius_init,
                         radius_increment = radius_increment,
                         x_stretch = x_stretch,
                         y_stretch = y_stretch,
                         offset_x = offset_x,
                         offset_y = offset_y,
                         skew = skew)
    
    n_cells <- nrow(locs)
    n_clusters <- length(unique(locs$Label))
    
    gcm <- matrix(data = rnorm(n_genes * n_cells, mean = 15, sd = 2), 
                  nrow = n_genes)
    de_genes <- sample(seq_len(nrow(gcm)), sum(n_deg))
    
    start <- c(1, cumsum(n_deg)[-length(n_deg)]+1)
    end <- cumsum(n_deg)
    
    for (i in seq_len(n_clusters)) {
        
        curr_de <- de_genes[start[i]:end[i]]
        cell_idx <- which(locs$Label == i)
        cell_idx <- sample(cell_idx, size = length(cell_idx) * de_dropout[i])
        gcm[curr_de, cell_idx] <- gcm[curr_de, cell_idx] +
            rnorm(length(curr_de)*length(cell_idx), mean = de_shift[i])
    }
    
    rownames(gcm) <- paste0('gene_', seq_len(nrow(gcm)))
    colnames(gcm) <- rownames(locs) <- paste0('cell_', seq_len(ncol(gcm)))
    meta <- locs[,3,drop=FALSE]
    return(list(gcm=gcm, locs=locs[,seq_len(2)], meta=meta))
}


#' @importFrom stats runif rnorm
generateBud <- function(n_circles, n_thickness, n_cells, 
                         radius_init, radius_increment,
                         x_stretch, y_stretch,
                         skew) {
      
    # First, get the number of cells per radius
    cells_per_radius <- round(n_cells / sum(n_thickness))
    n_cells <- cells_per_radius * sum(n_thickness)
    x <- rep(0, length(n_cells))
    y <- rep(0, length(n_cells))
    total_circles <- sum(n_thickness)
    lab <- rep(seq_len(n_circles), n_thickness * cells_per_radius)
    
    radius <- radius_init
    for (i in seq_len(total_circles)) {
        
        start <- (i-1)*cells_per_radius + 1
        end <- i*cells_per_radius 
        theta <- runif(cells_per_radius, min = 0, max = 2 * pi)
        x[start:end] <- radius * x_stretch * sin(theta + skew) + 
            rnorm(cells_per_radius)
        y[start:end] <- radius * y_stretch * cos(theta) + 
            rnorm(cells_per_radius)
        radius <- radius + radius_increment
    }
    
    list(x=x,y=y,Label=lab)
    
} 


generateLocs <- function(n_buds, n_circles, n_thickness, n_cells, 
                         radius_init, radius_increment, x_stretch, y_stretch,
                         offset_x, offset_y, skew) {
    
    locs <- vector(mode = 'list', length = n_buds)
    
    for (i in seq_len(n_buds)) {
        
        curr <- generateBud(n_circles = n_circles[i], 
                            n_thickness = n_thickness[[i]], 
                            n_cells = n_cells[i], 
                            radius_init = radius_init[i], 
                            radius_increment = radius_increment[i],
                            x_stretch = x_stretch[i], 
                            y_stretch = y_stretch[i],
                            skew = skew[i])
        locs[[i]] <- curr
    }
    
    lens <- vapply(locs, function(x) length(x$x), FUN.VALUE = numeric(1))
    res <- do.call(rbind.data.frame, locs)
    off_x <- rep(offset_x, lens)
    off_y <- rep(offset_y, lens)
    res$x <- res$x + off_x
    res$y <- res$y + off_y
    res$x <- res$x - mean(res$x)
    res$y <- res$y - mean(res$y)
    return(res)
    
}

#' @importFrom stats median aggregate
#' @importFrom data.table data.table setkey key 
binerize <- function(gcm, loc, mdist = NULL, normalize=TRUE) {
    
    if (is.null(mdist)) {
        knn <- kNN(loc, k = 1)
        mdist <- median(knn$dist) 
        message('Bin distance: ', mdist)
    }
    
    grp.x <- cut(loc$sdimx, 
                 breaks = c(min(loc$sdimx)-mdist, 
                            seq(min(loc$sdimx), max(loc$sdimx), by = mdist), 
                            max(loc$sdimx)+mdist), 
                 dig.lab = 10)
    grp.y <- cut(loc$sdimy, 
                 breaks = c(min(loc$sdimy)-mdist, 
                            seq(min(loc$sdimy), max(loc$sdimy), by = mdist), 
                            max(loc$sdimy)+mdist), 
                 dig.lab = 10)
    dloc <- data.table(loc,grp.x,grp.y)
    dloc$id <- seq_len(nrow(dloc))
    setkey(dloc, grp.x, grp.y)
    
    i <- .GRP <- id <- NULL
    dloc[, i := .GRP, by = key(dloc)]
    dloc$grid_x <- as.numeric(gsub('\\(|,.*', '', dloc$grp.x))
    dloc$grid_y <- as.numeric(gsub('\\(|,.*', '', dloc$grp.y))
    dloc <- dloc[order(id),]
    
    message('Getting bin locations')
    if (!normalize) mdist <- 1
    sloc <- data.frame(sdimx = dloc$grid_x/mdist, sdimy = dloc$grid_y/mdist)
    sloc <- sloc[!duplicated(sloc),]
    rownames(sloc) <- paste0('spot_', unique(dloc$i))
    
    message('Aggregating bin features')
    gcm <- t(gcm)
    gcm <- data.frame(cbind(grp=dloc$i, gcm))
    gsm <- aggregate(gcm[,2:ncol(gcm)], list(gcm$grp), mean)
    gsm <- t(gsm) 
    colnames(gsm) <- paste0('spot_', gsm[1,])
    gsm <- gsm[-1,]
    
    return(list(expression=gsm, locations=sloc, binstats=dloc))
    
}
