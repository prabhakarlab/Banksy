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
  bank@dim.reduction <- lapply(bank@dim.reduction, function(dim.red) {
    dim.red <- dim.red[rownames(dim.red) %in% cells,,drop=FALSE]
    dim.red
  })

  return(bank)
}
