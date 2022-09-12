#' S4 class BanksyObject
#'
#' @slot own.expr own expression
#' @slot nbr.expr neighbour expression
#' @slot harmonics azimuthal fourier harmonics
#' @slot custom.expr custom expression
#' @slot cell.locs cell locations
#' @slot meta.data metadata
#' @slot reduction dimension reductions
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
                        harmonics = 'assay',
                        custom.expr = 'assay',
                        cell.locs = 'assay',
                        meta.data = 'data.frame',
                        reduction = 'list'
        ))

## Constructor -----------------------------------------------------------------
#' @param own.expr own expression
#' @param nbr.expr neighbour expression
#' @param harmonics azimuthal fourier harmonics
#' @param custom.expr custom expression
#' @param cell.locs cell locations
#' @param meta.data metadata
#' @param reduction dimension reductions
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
#' 
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(
#'     own.expr = d$gcm, 
#'     cell.locs = d$locs, 
#'     meta.data = d$meta
#'     )
BanksyObject <- function(own.expr = NULL, nbr.expr = NULL, harmonics = NULL, 
                         custom.expr = NULL, cell.locs = NULL, meta.data = NULL,
                         reduction = NULL, genes.filter = 'intersect', 
                         min.cells.expressed = -1) {

    dimnames <- paste0('sdim', c('x', 'y','z'))
    object <- new('BanksyObject')

    if (is.null(own.expr)) stop('Provide expression matrix.')

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
        own.expr <- geneFilter(own.expr, genes.filter, min.cells.expressed)

        for (i in seq_len(nassays)) {
            gcm <- own.expr[[i]]
            locs <- cell.locs[[i]]
            locs <- locs[match(colnames(gcm), rownames(locs)), ]
            if (!identical(colnames(gcm), rownames(locs))) {
            stop('Ensure cell locations correspond to gene-cell matrix for ', 
                    names(own.expr)[i])
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
        ncells <- vapply(own.expr, ncol, FUN.VALUE = numeric(1))
        nCount <- unlist(lapply(own.expr, colSums))
        NODG <- unlist(lapply(own.expr, function(x) colSums(x > 0)))
        dnames <- rep(names(own.expr), ncells)
        mdata <- data.frame(cell_ID = cells, dataset = dnames, nCount = nCount, NODG = NODG)
        
    } else if (!is.null(own.expr)) {
        if (is.null(cell.locs)) stop('Provide cell location assay')
    
        own.expr <- as.matrix(own.expr)
        own.expr <- own.expr[rowSums(own.expr > 0) >= min.cells.expressed, ]
        
        names(cell.locs) <- dimnames[seq_len(ncol(cell.locs))]
        cell.locs <- cell.locs[match(colnames(own.expr), rownames(cell.locs)), ]
        
        if (!identical(colnames(own.expr), rownames(cell.locs))) {
            stop('Ensure cell locations correspond to gene-cell matrix.')
        }
        mdata <- data.frame(cell_ID = colnames(own.expr),
                            nCount = colSums(own.expr),
                            NODG = colSums(own.expr > 0))

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
    object@meta.data <- mdata[!duplicated(as.list(mdata))]
    object@cell.locs <- cell.locs
    return(object)
}

## Validity
setValidity('BanksyObject', function(object) {
    check <- TRUE
    if (!is(object@own.expr, 'assay')) check <- FALSE
    if (!is(object@nbr.expr, 'assay')) check <- FALSE
    if (!is(object@harmonics, 'assay')) check <- FALSE
    if (!is(object@custom.expr, 'assay')) check <- FALSE
    if (!is(object@cell.locs, 'assay')) check <- FALSE
    if (!is(object@meta.data, 'data.frame')) check <- FALSE
    if (!is(object@reduction, 'list')) check <- FALSE
    return(check)
})

## Getters ---------------------------------------------------------------------

#' @param object BanksyObject
#'
#' @describeIn BanksyObject-class getter own.expr
#' @exportMethod own.expr
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' mat <- own.expr(bank)
setGeneric('own.expr', function(object) standardGeneric('own.expr'))
#' @describeIn BanksyObject-class getter own.expr
setMethod('own.expr', signature(object = 'BanksyObject'),
            function(object) return(object@own.expr))

#' @describeIn BanksyObject-class getter nbr.expr
#' @exportMethod nbr.expr
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' mat <- own.expr(bank)
setGeneric('nbr.expr', function(object) standardGeneric('nbr.expr'))
#' @describeIn BanksyObject-class getter nbr.expr
setMethod('nbr.expr', signature(object = 'BanksyObject'),
            function(object) return(object@nbr.expr))

#' @describeIn BanksyObject-class getter harmonics
#' @exportMethod harmonics
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' lst <- harmonics(bank)
setGeneric('harmonics', function(object) standardGeneric('harmonics'))
#' @describeIn BanksyObject-class getter harmonics
setMethod('harmonics', signature(object = 'BanksyObject'),
          function(object) return(object@harmonics))

#' @describeIn BanksyObject-class getter custom.expr
#' @exportMethod custom.expr
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' mat <- nbr.expr(bank)
setGeneric('custom.expr', function(object) standardGeneric('custom.expr'))
#' @describeIn BanksyObject-class getter custom.expr
setMethod('custom.expr', signature(object = 'BanksyObject'),
            function(object) return(object@custom.expr))

#' @describeIn BanksyObject-class getter cell.locs
#' @exportMethod cell.locs
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' df <- cell.locs(bank)
setGeneric('cell.locs', function(object) standardGeneric('cell.locs'))
#' @describeIn BanksyObject-class getter cell.locs
setMethod('cell.locs', signature(object = 'BanksyObject'),
            function(object) return(object@cell.locs))

#' @describeIn BanksyObject-class getter meta.data
#' @exportMethod meta.data
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' df <- meta.data(bank)
setGeneric('meta.data', function(object) standardGeneric('meta.data'))
#' @describeIn BanksyObject-class getter meta.data
setMethod('meta.data', signature(object = 'BanksyObject'),
            function(object) return(object@meta.data))

#' @describeIn BanksyObject-class getter reduction
#' @exportMethod reduction
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' lst <- reduction(bank)
setGeneric('reduction', function(object) standardGeneric('reduction'))
#' @describeIn BanksyObject-class getter reduction
setMethod('reduction', signature(object = 'BanksyObject'),
            function(object) return(object@reduction))

#' @describeIn BanksyObject-class getter clust.names
#' @exportMethod clust.names
#' @examples 
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' vec <- clust.names(bank)
setGeneric('clust.names', function(object) standardGeneric('clust.names'))
#' @describeIn BanksyObject-class getter clust.names
setMethod('clust.names', signature(object = 'BanksyObject'),
            function(object) return(names(object@meta.data)[
                                    grepl('^clust', names(object@meta.data))]))

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
#' @describeIn BanksyObject-class setter for harmonics
#' @exportMethod harmonics<-
setGeneric('harmonics<-',
           function(object, value) standardGeneric('harmonics<-'))
#' @describeIn BanksyObject-class setter for harmonics
setReplaceMethod('harmonics', signature(object = 'BanksyObject', value = 'assay'),
                 function(object, value) {
                     object@harmonics <- value
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
#' @describeIn BanksyObject-class setter for reduction
#' @exportMethod reduction<-
setGeneric('reduction<-',
            function(object, value) standardGeneric('reduction<-'))
#' @describeIn BanksyObject-class setter for reduction
setReplaceMethod('reduction', signature(object = 'BanksyObject', value = 'list'),
                function(object, value) {
                    object@reduction <- value
                    validObject(object)
                    return(object)
                })

## Display ---------------------------------------------------------------------

#' @importFrom methods show
#' @describeIn BanksyObject-class show BanksyObject
#' @exportMethod show
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
        
        cat('Dimension reductions:', names(object@reduction), '\n')
        
        })

#' @importFrom utils head
#' @describeIn BanksyObject-class header for BanksyObject
#' @param x BanksyObject
#' @param n Number of records to display
#' @exportMethod head
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
