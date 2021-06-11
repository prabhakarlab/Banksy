#' S4 class BanksyObject
#'
#' @slot own.expr own expression
#' @slot nbr.expr neighbour expression
#' @slot own.norm.scaled.expr own norm/scaled expression
#' @slot nbr.norm.scaled.expr neighbour norm/scaled expression
#' @slot custom.expr custom expression
#' @slot cell.locs cell locations
#' @slot meta.data metadata
#' @slot dim.reduction dimension reductions
#'
#' @importClassesFrom Matrix dgCMatrix
#'
#' @name BanksyObject-class
#' @rdname BanksyObject-class
#'
#' @exportClass BanksyObject
setClassUnion('assay', c('data.frame', 'dgCMatrix', 'NULL'))
setClassUnion('listN', c('list', 'NULL'))
setClass('BanksyObject',
         slots = list(own.expr = 'assay',
                      nbr.expr = 'assay',
                      own.norm.scaled.expr = 'assay',
                      nbr.norm.scaled.expr = 'assay',
                      custom.expr = 'assay',
                      cell.locs = 'assay',
                      meta.data = 'data.frame',
                      dim.reduction = 'listN'
         ))

## Constructor -----------------------------------------------------------------
#' @param own.expr own expression
#' @param nbr.expr neighbour expression
#' @param own.norm.scaled.expr own norm/scaled expression
#' @param nbr.norm.scaled.expr neighbour norm/scaled expression
#' @param custom.expr custom expression
#' @param cell.locs cell locations
#' @param meta.data metadata
#' @param dim.reduction dimension reductions
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
           own.norm.scaled.expr = NULL,
           nbr.norm.scaled.expr = NULL,
           custom.expr = NULL,
           cell.locs = NULL,
           meta.data = data.frame(),
           dim.reduction = NULL) {
    object <- new('BanksyObject',
                  own.expr = own.expr,
                  nbr.expr = nbr.expr,
                  own.norm.scaled.expr = own.norm.scaled.expr,
                  nbr.norm.scaled.expr = nbr.norm.scaled.expr,
                  custom.expr = custom.expr,
                  cell.locs = cell.locs,
                  meta.data = meta.data,
                  dim.reduction = dim.reduction)
    if (!is.null(own.expr)) {
      meta.data(object) <- data.frame(cell_ID = colnames(own.expr))
    }
    return(object)
  }

## Validity
setValidity('BanksyObject', function(object) {
  check <- TRUE
  if (!is(object@own.expr, 'assay')) check <- FALSE
  if (!is(object@nbr.expr, 'assay')) check <- FALSE
  if (!is(object@own.norm.scaled.expr, 'assay')) check <- FALSE
  if (!is(object@nbr.norm.scaled.expr, 'assay')) check <- FALSE
  if (!is(object@custom.expr, 'assay')) check <- FALSE
  if (!is(object@cell.locs, 'assay')) check <- FALSE
  if (!is(object@meta.data, 'data.frame')) check <- FALSE
  if (!is(object@dim.reduction, 'listN')) check <- FALSE
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

#' @describeIn BanksyObject-class getter own.norm.scaled.expr
#' @exportMethod own.norm.scaled.expr
setGeneric('own.norm.scaled.expr', function(object) standardGeneric('own.norm.scaled.expr'))
#' @describeIn BanksyObject-class getter own.norm.scaled.expr
setMethod('own.norm.scaled.expr', signature(object = 'BanksyObject'),
          function(object) return(object@own.norm.scaled.expr))

#' @describeIn BanksyObject-class getter nbr.norm.scaled.expr
#' @exportMethod nbr.norm.scaled.expr
setGeneric('nbr.norm.scaled.expr', function(object) standardGeneric('nbr.norm.scaled.expr'))
#' @describeIn BanksyObject-class getter nbr.norm.scaled.expr
setMethod('nbr.norm.scaled.expr', signature(object = 'BanksyObject'),
          function(object) return(object@nbr.norm.scaled.expr))

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
#' @describeIn BanksyObject-class setter for own.norm.scaled.expr
#' @exportMethod own.norm.scaled.expr<-
setGeneric('own.norm.scaled.expr<-',
           function(object, value) standardGeneric('own.norm.scaled.expr<-'))
#' @describeIn BanksyObject-class setter for own.norm.scaled.expr
setReplaceMethod('own.norm.scaled.expr', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@own.norm.scaled.expr <- value
                   validObject(object)
                   return(object)
                 })

#' @importFrom methods validObject
#' @describeIn BanksyObject-class setter for nbr.norm.scaled.expr
#' @exportMethod nbr.norm.scaled.expr<-
setGeneric('nbr.norm.scaled.expr<-',
           function(object, value) standardGeneric('nbr.norm.scaled.expr<-'))
#' @describeIn BanksyObject-class setter for nbr.norm.scaled.expr
setReplaceMethod('nbr.norm.scaled.expr', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                   object@nbr.norm.scaled.expr <- value
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
            cat('class:', class(object), '\n')

            ## Raw expression
            if (is.null(object@own.expr)) cat('own expression: NULL\n')
            else cat('own expression:', ncol(object@own.expr), 'cells', nrow(object@own.expr), 'features\n')

            if (is.null(object@nbr.expr)) cat('neighbour expression: NULL\n')
            else cat('neighbour expression:', ncol(object@nbr.expr), 'cells', nrow(object@nbr.expr), 'features\n')

            if (is.null(object@own.norm.scaled.expr)) cat('own norm. scaled expression: NULL\n')
            else cat('own norm. scaled expression:', ncol(object@own.norm.scaled.expr), 'cells', nrow(object@own.norm.scaled.expr), 'features\n')

            if (is.null(object@nbr.norm.scaled.expr)) cat('neighbour norm. scaled expression: NULL\n')
            else cat('neighbour norm. scaled expression:', ncol(object@nbr.norm.scaled.expr), 'cells', nrow(object@nbr.norm.scaled.expr), 'features\n')

            if (is.null(object@cell.locs)) cat('dimensions: NULL\n')
            else cat('spatial dimensions:', colnames(object@cell.locs), '\n')

            cat('metadata names:', colnames(object@meta.data), '\n')

            cat('dimension reductions:', names(object@dim.reduction), '\n')
          })

#' @param x BanksyObject
#' @aliases head,BanksyObject-class
#' @exportMethod head
#' @describeIn BanksyObject-class head for banksy object
setGeneric('head', function(x) standardGeneric('head'))
setMethod('head', signature(x = 'BanksyObject'),
          function(x, n=5) {
            cat('own expression:\n')
            rd <- min(nrow(x@own.expr), n)
            cd <- min(ncol(x@own.expr), n)
            print(x@own.expr[seq_len(rd),seq_len(cd),drop=FALSE])

            cat('\nneighbour expression:\n')
            rd <- min(nrow(x@nbr.expr), n)
            cd <- min(ncol(x@nbr.expr), n)
            print(x@nbr.expr[seq_len(rd),seq_len(cd),drop=FALSE])

            cat('\nown normalized scaled expression:\n')
            rd <- min(nrow(x@own.norm.scaled.expr), n)
            cd <- min(ncol(x@own.norm.scaled.expr), n)
            print(x@own.norm.scaled.expr[seq_len(rd),seq_len(cd),drop=FALSE])

            cat('\nneighbour normalized scaled expression:\n')
            rd <- min(nrow(x@nbr.norm.scaled.expr), n)
            cd <- min(ncol(x@nbr.norm.scaled.expr), n)
            print(x@nbr.norm.scaled.expr[seq_len(rd),seq_len(cd),drop=FALSE])

            cat('\ncell locations:\n')
            rd <- min(nrow(x@cell.locs), n)
            cd <- min(ncol(x@cell.locs), n)
            print(x@cell.locs[seq_len(rd),seq_len(cd),drop=FALSE])

            cat('\nmetadata:\n')
            rd <- min(nrow(x@meta.data), n)
            cd <- min(ncol(x@meta.data), n)
            print(x@meta.data[seq_len(rd),seq_len(cd),drop=FALSE])

          })

