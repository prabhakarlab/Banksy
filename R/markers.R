
#' Wilcoxon two-sided tests to find markers (based on Seurat)
#'
#' @param bank BanksyObject with unscaled expression
#' @param id1 identity class 1
#' @param id2 identity class 2
#' @param assay assay to find markers with
#' \itemize{
#'   \item{own.expr: own expression (default)}
#'   \item{nbr.expr: neighbor expression}
#' }
#' @param by metadata column to define identity classes (default uses latest clustering)
#' @param q.val.thres q-value threshold to filter out markers (default: 0.01)
#' @param express.thres count threshold for gene to be considered expressed in cell (default: 1)
#' @param only.pos return markers with positive log2 fold change (i.e. for identity class 1)
#' @param pseudocount.use pseudocount for fold change computation
#'
#' @importFrom qvalue qvalue
#'
#' @return dataframe of markers
#' @details returns six columns:
#' \itemize{
#'   \item{p.val: pvalue returned by the Wilcoxon rank-sum test}
#'   \item{pct.1: percentage of cells in id1 expressing gene (based on express.thres)}
#'   \item{pct.2: percentage of cells in id2 expressing gene (based on express.thres)}
#'   \item{log.FC: log2 fold change between mean of id1 and id2 with pseudocount}
#'   \item{q.val: qvalue returned by correction of pvalues}
#'   \item{gene: marker name}
#' }
#'
#' @export
#' 
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' markers <- RunWilcox(bank, 1, 2, by = 'Label')
#'
RunWilcox <- function(bank,
                      id1,
                      id2,
                      assay = 'own.expr',
                      by = NULL,
                      q.val.thres = 0.01,
                      express.thres = 1,
                      only.pos = FALSE,
                      pseudocount.use = 0.1) {
    if (is.null(by))
        by = clust.names(bank)[length(clust.names(bank))]
    
    id1.cells = bank@meta.data$cell_ID[bank@meta.data[[by]] %in% id1]
    id2.cells = bank@meta.data$cell_ID[bank@meta.data[[by]] %in% id2]
    
    id1.mat = slot(bank, assay)[, id1.cells]
    id2.mat = slot(bank, assay)[, id2.cells]
    
    results = lapply(seq_len(nrow(id1.mat)), function(i) {
        out = wilcox.test(id1.mat[i, ], id2.mat[i, ], alternative = 'two.sided')
        pct.1 = sum(id1.mat[i, ] > express.thres) / length(id1.mat[i, ])
        pct.2 = sum(id2.mat[i, ] > express.thres) / length(id2.mat[i, ])
        list(p.val = out$p.value,
             pct.1 = pct.1,
             pct.2 = pct.2)
    })
    
    data.1 <-
        apply(id1.mat, 1, function(x)
            log2(mean(x + pseudocount.use)))
    data.2 <-
        apply(id2.mat, 1, function(x)
            log2(mean(x + pseudocount.use)))
    total.diff <- data.1 - data.2
    
    de = do.call(rbind.data.frame, results)
    de$log.FC = total.diff
    qq = tryCatch(
        qvalue::qvalue(de$p.val),
        error = function(cond)
            qvalue::qvalue(de$p.val, lambda = 0)
    )
    de$q.val = qq$qvalues
    de$gene = rownames(id1.mat)
    if (!is.null(q.val.thres))
        de = de[de$q.val < q.val.thres, ]
    de = de[order(de$log.FC, decreasing = TRUE), ]
    if (only.pos)
        de = de[de$log.FC > 0, ]
    rownames(de) = NULL
    de
    
}

#' Wilcoxon two-sided tests to find markers for all identity class (based on Seurat)
#'
#' @param bank BanksyObject with unscaled expression
#' @param assay assay to find markers with
#' \itemize{
#'   \item{own.expr: own expression (default)}
#'   \item{nbr.expr: neighbor expression}
#' }
#' @param by metadata column to define identity classes (default uses latest clustering)
#' @param q.val.thres q-value threshold to filter out markers (default: 0.01)
#' @param express.thres count threshold for gene to be considered expressed in cell (default: 1)
#' @param only.pos return markers with positive log2 fold change (i.e. for identity class 1)
#' @param pseudocount.use pseudocount for fold change computation
#' @param verbose verbose
#'
#' @return dataframe of markers
#' @details returns seven columns:
#' \itemize{
#'   \item{p.val: pvalue returned by the Wilcoxon rank-sum test}
#'   \item{pct.1: percentage of cells in id1 expressing gene (based on express.thres)}
#'   \item{pct.2: percentage of cells in id2 expressing gene (based on express.thres)}
#'   \item{log.FC: log2 fold change between mean of id1 and id2 with pseudocount}
#'   \item{q.val: qvalue returned by correction of pvalues}
#'   \item{gene: marker name}
#'   \item{cluster: cluster for marker}
#' }
#' 
#' @export
#' 
#' @examples
#' # Generate a simulated dataset
#' d <- simulateDataset()
#' bank <- BanksyObject(own.expr = d$gcm, cell.locs = d$locs, meta.data = d$meta)
#' markers <- RunAllWilcox(bank, by = 'Label')
#'
RunAllWilcox <- function(bank,
                         assay = 'own.expr',
                         by = NULL,
                         q.val.thres = 0.01,
                         express.thres = 1,
                         only.pos = TRUE,
                         pseudocount.use = 0.1,
                         verbose = TRUE) {
    if (is.null(by))
        by = clust.names(bank)[length(clust.names(bank))]
    
    ids = sort(unique(bank@meta.data[[by]]))
    
    de.lst = lapply(ids, function(i) {
        if (verbose)
            message('Find markers for cluster ', i)
        RunWilcox(
            bank,
            assay = assay,
            by = by,
            id1 = i,
            id2 = setdiff(ids, i),
            q.val.thres = q.val.thres,
            express.thres = express.thres,
            only.pos = only.pos,
            pseudocount.use = pseudocount.use
        )
    })
    
    de = do.call(rbind.data.frame, de.lst)
    de$cluster = rep(ids, sapply(de.lst, nrow))
    de
    
}
