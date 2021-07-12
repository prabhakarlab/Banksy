#' S4 class BanksyObject
#'
#' @slot own.expr own expression
#' @slot nbr.expr neighbour expression
#' @slot custom.expr custom expression
#' @slot cell.locs cell locations
#' @slot meta.data metadata
#' @slot dim.reduction dimension reductions
#'
# #' @importClassesFrom Matrix dgCMatrix
#'
#' @name BanksyObject-class
#' @rdname BanksyObject-class
#'
#' @exportClass BanksyObject
setClassUnion('assay', c('data.frame', 'matrix', 'list', 'NULL'))
setClass('BanksyObject',
         slots = list(own.expr = 'assay',
                      nbr.expr = 'assay',
                      custom.expr = 'assay',
                      cell.locs = 'assay',
                      meta.data = 'data.frame',
                      dim.reduction = 'assay'
         ))

## Constructor -----------------------------------------------------------------
#' @param own.expr own expression
#' @param nbr.expr neighbour expression
#' @param custom.expr custom expression
#' @param cell.locs cell locations
#' @param meta.data metadata
#' @param dim.reduction dimension reductions
#' @param genes.filter intersect or union
#' @param min.cells.expressed min number of cells a gene is expressed in
#'   (to be used with intersect)
#'
#' @name BanksyObject
#' @rdname BanksyObject-class
#'
#' @importFrom methods new
#'
#' @return BanksyObject
#'
#' @export
BanksyObject <-
  function(own.expr = NULL,
           nbr.expr = NULL,
           custom.expr = NULL,
           cell.locs = NULL,
           meta.data = NULL,
           dim.reduction = NULL,
           genes.filter = 'intersect',
           min.cells.expressed = -1) {

    dimnames <- paste0('sdim', c('x', 'y','z'))
    object <- new('BanksyObject')

    # Multiple dataset case
    if (is.list(own.expr)) {

      if (length(own.expr) != length(cell.locs)) {
        stop('Require equal number of expression and location assays')
      }

      # Names
      nassays <- length(own.expr)
      if (is.null(names(own.expr))) {
        names(own.expr) <- paste0('d', seq_len(nassays))
        names(cell.locs) <- paste0('d', seq_len(nassays))

      }

      # Filter genes
      own.expr <- geneFilter(own.expr,
                             genes.filter,
                             min.cells.expressed)

      for (i in seq_len(nassays)) {
        gcm <- own.expr[[i]]
        locs <- cell.locs[[i]]
        locs <- locs[match(colnames(gcm), rownames(locs)), ]
        if (!identical(colnames(gcm), rownames(locs))) {
          stop(paste0('Ensure cell locations correspond to gene-cell matrix for ', names(own.expr)[i]))
        }
        gcm <- as.matrix(gcm)
        colnames(gcm) <- paste0(names(own.expr)[i], '_', colnames(gcm))
        names(locs) <- dimnames[seq_len(ncol(locs))]
        rownames(locs) <- colnames(gcm)
        own.expr[[i]] <- gcm
        cell.locs[[i]] <- locs
      }

      # Metadata
      cells <- unlist(lapply(own.expr, colnames))
      ncells <- sapply(own.expr, ncol)
      dnames <- rep(names(own.expr), ncells)
      mdata <- data.frame(cell_ID = cells,
                          dataset = dnames)
    } else if (!is.null(own.expr)) {

      if (is.null(cell.locs)) {
        stop('Provide cell location assay')
      }

      own.expr <- as.matrix(own.expr)
      own.expr <- own.expr[rowSums(own.expr > 0) >= min.cells.expressed, ]

      names(cell.locs) <- dimnames[seq_len(ncol(cell.locs))]
      cell.locs <- cell.locs[match(colnames(own.expr), rownames(cell.locs)), ]

      if (!identical(colnames(own.expr), rownames(cell.locs))) {
        stop('Ensure cell locations correspond to gene-cell matrix.')
      }
      mdata <- data.frame(cell_ID = colnames(own.expr))

    } else {
      return(object)
    }

    if (!is.null(meta.data)) {
      if (nrow(meta.data) != nrow(mdata)) {
        stop('Wrong number of rows in meta.data.
              Ensure nrows == number of cells across all datasets.')
      } else {
        mdata <- cbind(mdata, meta.data)
      }
    }

    object@own.expr <- own.expr
    object@meta.data <- mdata
    object@cell.locs <- cell.locs
    return(object)
  }

## Validity
setValidity('BanksyObject', function(object) {
  check <- TRUE
  if (!is(object@own.expr, 'assay')) check <- FALSE
  if (!is(object@nbr.expr, 'assay')) check <- FALSE
  if (!is(object@custom.expr, 'assay')) check <- FALSE
  if (!is(object@cell.locs, 'assay')) check <- FALSE
  if (!is(object@meta.data, 'data.frame')) check <- FALSE
  if (!is(object@dim.reduction, 'assay')) check <- FALSE
  return(check)
})

## Getters ---------------------------------------------------------------------

#' @param object BanksyObject
#'
#' @describeIn BanksyObject-class getter own.expr
#' @exportMethod own.expr
setGeneric('own.expr', function(object) standardGeneric('own.expr'))
#' @describeIn BanksyObject-class getter own.expr
setMethod('own.expr', signature(object = 'BanksyObject'),
          function(object) return(object@own.expr))

#' @describeIn BanksyObject-class getter nbr.expr
#' @exportMethod nbr.expr
setGeneric('nbr.expr', function(object) standardGeneric('nbr.expr'))
#' @describeIn BanksyObject-class getter nbr.expr
setMethod('nbr.expr', signature(object = 'BanksyObject'),
          function(object) return(object@nbr.expr))

#' @describeIn BanksyObject-class getter custom.expr
#' @exportMethod custom.expr
setGeneric('custom.expr', function(object) standardGeneric('custom.expr'))
#' @describeIn BanksyObject-class getter custom.expr
setMethod('custom.expr', signature(object = 'BanksyObject'),
          function(object) return(object@custom.expr))

#' @describeIn BanksyObject-class getter cell.locs
#' @exportMethod cell.locs
setGeneric('cell.locs', function(object) standardGeneric('cell.locs'))
#' @describeIn BanksyObject-class getter cell.locs
setMethod('cell.locs', signature(object = 'BanksyObject'),
          function(object) return(object@cell.locs))

#' @describeIn BanksyObject-class getter meta.data
#' @exportMethod meta.data
setGeneric('meta.data', function(object) standardGeneric('meta.data'))
#' @describeIn BanksyObject-class getter meta.data
setMethod('meta.data', signature(object = 'BanksyObject'),
          function(object) return(object@meta.data))

#' @describeIn BanksyObject-class getter dim.reduction
#' @exportMethod dim.reduction
setGeneric('dim.reduction', function(object) standardGeneric('dim.reduction'))
#' @describeIn BanksyObject-class getter dim.reduction
setMethod('dim.reduction', signature(object = 'BanksyObject'),
          function(object) return(object@dim.reduction))

## Setters ---------------------------------------------------------------------

#' @param object BanksyObject
#' @param value to be assigned
#'
#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for own.expr
#' @exportMethod own.expr<-
setGeneric('own.expr<-',
           function(object, value) standardGeneric('own.expr<-'))
#' @describeIn BanksyObject-class setter for own.expr
setReplaceMethod('own.expr', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@own.expr <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for nbr.expr
#' @exportMethod nbr.expr<-
setGeneric('nbr.expr<-',
           function(object, value) standardGeneric('nbr.expr<-'))
#' @describeIn BanksyObject-class setter for nbr.expr
setReplaceMethod('nbr.expr', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@nbr.expr <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for custom.expr
#' @exportMethod custom.expr<-
setGeneric('custom.expr<-',
           function(object, value) standardGeneric('custom.expr<-'))
#' @describeIn BanksyObject-class setter for custom.expr
setReplaceMethod('custom.expr', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@custom.expr <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for cell.locs
#' @exportMethod cell.locs<-
setGeneric('cell.locs<-',
           function(object, value) standardGeneric('cell.locs<-'))
#' @describeIn BanksyObject-class setter for cell.locs
setReplaceMethod('cell.locs', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@cell.locs <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for meta.data
#' @exportMethod meta.data<-
setGeneric('meta.data<-',
           function(object, value) standardGeneric('meta.data<-'))
#' @describeIn BanksyObject-class setter for meta.data
setReplaceMethod('meta.data', signature(object = 'BanksyObject', value = 'data.frame'),
                 function(object, value) {
                   object@meta.data <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for dim.reduction
#' @exportMethod dim.reduction<-
setGeneric('dim.reduction<-',
           function(object, value) standardGeneric('dim.reduction<-'))
#' @describeIn BanksyObject-class setter for dim.reduction
setReplaceMethod('dim.reduction', signature(object = 'BanksyObject', value = 'list'),
                 function(object, value) {
                   object@dim.reduction <- value
                   validObject(object)
                   return(object)
                 })

## Display ---------------------------------------------------------------------

setMethod('show', signature(object = 'BanksyObject'),
          function(object){

            cat('Object of class', class(object), '\n')
            nassays <- length(object@own.expr)

            ## Raw expression
            if (is.null(object@own.expr)) {
              cat('Assays: NULL\n')
            } else {
              if (is.list(object@own.expr)) {
                cat('Number of assays: ', nassays, '\n', sep='')
                dnames <- names(object@own.expr)
                for (i in seq_len(nassays)) {
                  cat(dnames[i], ': ',
                      ncol(object@own.expr[[i]]), ' cells ',
                      nrow(object@own.expr[[i]]), ' features\n',
                      sep='')
                }
              } else {
                cat('Assay with', ncol(object@own.expr), 'cells',
                    nrow(object@own.expr), 'features\n')
              }
            }

            if (is.null(object@cell.locs)) {
              cat('Spatial dimensions: NULL\n')
            } else if (is.list(object@own.expr)) {
              cat('Spatial dimensions:\n')
              dnames <- names(object@cell.locs)
              for (i in seq_len(nassays)) {
                cat(dnames[i], ': ',
                    paste(names(object@cell.locs[[i]]), collapse = ' '),
                    '\n', sep='')
              }
            } else {
              cat('Spatial dimensions:', names(object@cell.locs), '\n')
            }

            cat('Metadata names:', colnames(object@meta.data), '\n')

            cat('Dimension reductions:', names(object@dim.reduction), '\n')

          })

#' @param x BanksyObject
#' @aliases head,BanksyObject-class
#' @exportMethod head
#' @describeIn BanksyObject-class head for banksy object
setGeneric('head', function(x) standardGeneric('head'))
setMethod('head', signature(x = 'BanksyObject'),
          function(x, n=5) {

            if (is.list(x@own.expr)) {
              cat('own expression:\n')
              rd <- min(nrow(x@own.expr[[1]]), n)
              cd <- min(ncol(x@own.expr[[1]]), n)
              print(x@own.expr[[1]][seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\nneighbour expression:\n')
              rd <- min(nrow(x@nbr.expr[[1]]), n)
              cd <- min(ncol(x@nbr.expr[[1]]), n)
              print(x@nbr.expr[[1]][seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\ncell locations:\n')
              rd <- min(nrow(x@cell.locs[[1]]), n)
              cd <- min(ncol(x@cell.locs[[1]]), n)
              print(x@cell.locs[[1]][seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\nmetadata:\n')
              rd <- min(nrow(x@meta.data), n)
              cd <- min(ncol(x@meta.data), n)
              print(x@meta.data[seq_len(rd),seq_len(cd),drop=FALSE])
            } else {
              cat('own expression:\n')
              rd <- min(nrow(x@own.expr), n)
              cd <- min(ncol(x@own.expr), n)
              print(x@own.expr[seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\nneighbour expression:\n')
              rd <- min(nrow(x@nbr.expr), n)
              cd <- min(ncol(x@nbr.expr), n)
              print(x@nbr.expr[seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\ncell locations:\n')
              rd <- min(nrow(x@cell.locs), n)
              cd <- min(ncol(x@cell.locs), n)
              print(x@cell.locs[seq_len(rd),seq_len(cd),drop=FALSE])

              cat('\nmetadata:\n')
              rd <- min(nrow(x@meta.data), n)
              cd <- min(ncol(x@meta.data), n)
              print(x@meta.data[seq_len(rd),seq_len(cd),drop=FALSE])
            }
          })

