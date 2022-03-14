
#' Convert a SingleCellExperiment to BanksyObject
#' 
#' @param se object of class SingleCellExperiment or SpatialExperiment
#' @param expr.assay expression assay to extract
#' @param coord.colnames colnames of expression in metadata from SCE 
#' @param features features to extract from expression assay
#' 
#' @importFrom SummarizedExperiment assays colData 
#' 
#' @return BanksyObject
#' 
#' @export
#' 
#' @examples 
#' data <- simulateDataset()
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'     assays = list(counts = data$gcm), colData = data$locs)
#' bank <- asBanksyObject(sce, expr.assay = 'counts', 
#'                             coord.colnames = c('x', 'y'))
asBanksyObject <- function(se, expr.assay = 'counts', features = NULL, 
                           coord.colnames = c('sdimx', 'sdimy')) {
    cx <- class(se)
    
    # Own expression
    if (!(expr.assay %in% names(assays(se)))) stop(
        'Assay ', expr.assay, ' not found')
    gcm <- assays(se)[[expr.assay]]
    
    # Filter features
    if (!is.null(features)) {
        found <- features %in% rownames(gcm)
        if (!all(found)) {
            miss <- paste(features[!found], collapse = ' ')
            warning('Dropping features ', miss)
            features <- features[found]
        }
        gcm <- gcm[features, ]
    }
    
    # Meta data
    meta <- data.frame(colData(se))
    if ('cell_ID' %in% colnames(meta)) warning(
        'Reassigning cell_ID column in meta.data to colnames of SCE')
    meta$cell_ID <- colnames(se)
    
    # Coordinates 
    if (cx == 'SingleCellExperiment') {
        # Coord colnames
        if (!(length(coord.colnames) %in% 2:3)) stop(
            'Only 2 or 3 spatial dimensions allowed.')
        if (all(coord.colnames %in% colnames(meta))) {
            locs <- meta[, coord.colnames, drop = FALSE]
        }
    }
    
    bank <- BanksyObject(own.expr = gcm, cell.locs = locs, meta.data = meta)
    return(bank)
}
