
#' Subset a BanksyObject
#' @param x a BanksyObject
#' @param cells cells to filter by
#' @param dims dimensions to filter by - must correspond to valid columns in
#'  cell.locs
#' @param features genes to filter by
#' @param metadata metadata to filter by - must be valid column in meta.data
#' @param dataset dataset to subset dimensions by
#'
#' @importFrom rlang enquo quo_get_expr
#'
#' @return subset BanksyObject
#'
#' @export
#' 
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' # Subset cells with Label = 1
#' bank <- SubsetBanksy(bank, metadata = Label == 1)
#' 
SubsetBanksy <- function(
  x, cells=NULL, features=NULL, dims=TRUE, metadata=TRUE, dataset=NULL) {

  nfeaturesBef <- nrow(x@own.expr)
  ncellsBef <- ncol(x@own.expr)

  ## Enquote all input conditions
  dims <- rlang::enquo(dims)
  metadata <- rlang::enquo(metadata)

  ## Subset dimensions
  x <- subsetLocations(x, dims, dataset)

  ## Subset cells
  x <- subsetCells(x, cells)

  ## Subset metadata
  x@meta.data <- subset(x@meta.data, subset = eval(rlang::quo_get_expr(metadata)))
  ## Subset features
  x <- subsetFeatures(x, features)

  ## Consistent Cell filtering
  survivingCells <- subsetConsistent(x)
  x <- subsetCells(x, survivingCells)

  ## Filter cell locations, metadata and dim reductions based on surviving cells
  x@meta.data <- x@meta.data[x@meta.data$cell_ID %in% survivingCells,,drop=FALSE]
  r1 <- lapply(x@reduction[grepl('pca', names(x@reduction))], function(red) {
    red$x <- red$x[rownames(red$x) %in% survivingCells,,drop=FALSE]
    red
  })
  r2 <- lapply(x@reduction[!grepl('pca', names(x@reduction))], function(red) {
    red <- red[rownames(red) %in% survivingCells,,drop=FALSE]
    red
  })
  x@reduction <- c(r1, r2)
  if (is.list(x@own.expr)) {
    x@cell.locs <- lapply(x@cell.locs, function(x) {
      keep <- which(rownames(x) %in% survivingCells)
      x <- x[keep,,drop=FALSE]
    })
  } else {
    keep <- which(rownames(x@cell.locs) %in% survivingCells)
    x@cell.locs <- x@cell.locs[keep,,drop=FALSE]
  }

  ## For multiple datasets - remove empty datasets
  if (is.list(x@own.expr)) {
    x <- cleanSubset(x)
  }

  return(x)
}

#' Split a BanksyObject by metadata column
#'
#' @param bank BanksyObject
#' @param by metadata column
#' @param names names for new datasets
#'
#' @return BanksyObject
#'
#' @export
#' 
#' @examples 
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' # Split BanksyObject into multiple datasets based on label
#' bank <- SplitBanksy(bank, by = 'Label')
#' bank
#'
SplitBanksy <- function(bank, by, names = NULL) {

  if (is.list(bank@own.expr)) {
    warning('SplitBanksy only operates on BanksyObjects with a single dataset')
    return(bank)
  }
  if (!(by %in% names(bank@meta.data)))  stop(by, ' not in metadata.')
  groups <- unique(bank@meta.data[[by]])
  nGroups <- length(groups)

  if (is.null(names)) {
    names <- paste(by, groups, sep = "_")
  } else if (length(names) != nGroups) {
    stop('Invalid number of names.')
  }

  own.expr <- nbr.expr <- cell.locs <- vector('list', length = nGroups)
  if (is.null(bank@nbr.expr)) nbr.expr <- NULL

  for (i in seq_len(nGroups)) {

    currGroup <- groups[i]
    message('Processing ', by, ' ', currGroup)
    cells <- bank@meta.data$cell_ID[bank@meta.data[[by]] == currGroup]
    gcm <- bank@own.expr[, cells]
    ncm <- bank@nbr.expr[, cells]
    loc <- bank@cell.locs[cells, ]
    colnames(gcm) <- paste0(names[i], '_', cells)
    if (!is.null(ncm)) colnames(ncm) <- paste0(names[i], '_', cells)
    rownames(loc) <- paste0(names[i], '_', cells)
    own.expr[[i]] <- gcm
    nbr.expr[[i]] <- ncm
    cell.locs[[i]] <- loc
  }

  names(own.expr)  <- names(cell.locs) <- names
  if (!is.null(nbr.expr)) names(nbr.expr) <- names

  bank@own.expr <- own.expr
  bank@nbr.expr <- nbr.expr
  bank@cell.locs <- cell.locs
  bank@reduction <- bank@reduction
  bank@meta.data <- bank@meta.data
  bank@meta.data$cell_ID <- unlist(lapply(cell.locs, row.names))

  return(bank)
}

# Index of helper functions and where they are called
# subsetDims ---------------- SubsetBanksy
# subsetLocations ----------- SubsetBanksy
# subsetCells --------------- SubsetBanksy
# subsetFeatures ------------ SubsetBanksy
# subsetConsistent ---------- SubsetBanksy

#' @importFrom rlang quo_get_expr
subsetDims <- function(locs, dims) {
  locs <- subset(data.frame(locs), subset = eval(rlang::quo_get_expr(dims)))
  return(locs)
}

subsetLocations <- function(x, dims=TRUE, dataset=NULL) {

  if (is.list(x@own.expr)) {

    if (is.null(dataset)) {
      message('Dataset not specified. Subsetting all data sets by dimension.')
      x@cell.locs <- lapply(x@cell.locs, function(x) {
        subsetDims(x, dims)
      })
    } else {
      if (!(dataset %in% names(x@cell.locs))) {
        stop(paste0('Dataset ', dataset, ' not found. One of ',
                    paste(names(x@cell.locs), collapse = ' ')))
      }
      assayid <- which(names(x@cell.locs) == dataset)
      assay <- x@cell.locs[[assayid]]
      x@cell.locs[[assayid]] <- subsetDims(assay, dims)
    }

  } else {
    x@cell.locs <- subsetDims(x@cell.locs, dims)
  }

  return(x)
}

subsetCells <- function(x, cells) {

  if (is.null(cells)) {
    return(x)
  }

  if (is.list(x@own.expr)) {
    x@own.expr <- lapply(x@own.expr, function(x) {
      keep <- which(colnames(x) %in% cells)
      x <- x[,keep, drop = FALSE]
      x
    })
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- lapply(x@nbr.expr, function(x) {
        keep <- which(colnames(x) %in% cells)
        x <- x[,keep, drop = FALSE]
        x
      })
    }
  } else {
    keep <- which(colnames(x@own.expr) %in% cells)
    x@own.expr <- x@own.expr[, keep, drop = FALSE]
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- x@nbr.expr[, keep, drop = FALSE]
    }
  }
  return(x)
}

subsetFeatures <- function(x, features) {

  if (is.null(features)) {
    return(x)
  }

  if (is.list(x@own.expr)) {
    keep <- which(rownames(x@own.expr[[1]]) %in% features )
    x@own.expr <- lapply(x@own.expr, function(x) {
      x <- x[keep,,drop = FALSE]
      x
    })
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- lapply(x@nbr.expr, function(x) {
        x <- x[keep,,drop=FALSE]
        x
      })
    }

  } else {
    keep <- which(rownames(x@own.expr) %in% features)
    x@own.expr <- x@own.expr[keep,,drop=FALSE]
    if (!is.null(x@nbr.expr)) {
      x@nbr.expr <- x@nbr.expr[keep,,drop=FALSE]
    }
  }
  return(x)
}

subsetConsistent <- function(x) {

  cellMeta <- x@meta.data$cell_ID

  if (is.list(x@own.expr)) {
    cellCell <- as.character(unlist(lapply(x@own.expr, colnames)))
    cellLocs <- as.character(unlist(lapply(x@cell.locs, rownames)))
  } else {
    cellCell <- colnames(x@own.expr)
    cellLocs <- rownames(x@cell.locs)
  }

  cellOut <- intersect(cellMeta, intersect(cellCell, cellLocs))
  return(cellOut)
}

cleanSubset <- function(bank) {

  keep <- unlist(lapply(bank@own.expr, function(x) all(dim(x) != 0)))
  bank@own.expr <- bank@own.expr[keep]
  bank@nbr.expr <- bank@nbr.expr[keep]
  bank@cell.locs <- bank@cell.locs[keep]

  cells <- unlist(lapply(bank@own.expr, colnames))

  bank@meta.data <- bank@meta.data[bank@meta.data$cell_ID %in% cells,,drop=FALSE]
  r1 <- lapply(bank@reduction[grepl('pca', names(bank@reduction))], 
               function(red) {
                   red$x <- red$x[rownames(red$x) %in% cells,,drop=FALSE]
                   red})
  r2 <- lapply(bank@reduction[!grepl('pca', names(bank@reduction))], 
               function(red) {
                   red <- red[rownames(red) %in% cells,,drop=FALSE]
                   red})
  bank@reduction <- c(r1, r2)
  
  
  return(bank)
}
