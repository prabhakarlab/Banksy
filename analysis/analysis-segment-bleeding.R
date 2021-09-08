
##### Libraries ----------------------------------------------------------------

rm(list=ls())
graphics.off()

library(Banksy)
library(tictoc)
library(gridExtra)
library(ggplot2)
library(data.table)
library(ComplexHeatmap)
library(Seurat)
library(dbscan)
library(plyr)
library(circlize)
library(qvalue)
library(schelpers)

##### Banksy -------------------------------------------------------------------

# bank <- readRDS('moffitt-banksy-object-animal1.rds')
# bank <- ComputeBanksy(bank)
# bank <- ScaleBanksy(bank)
# bank <- ClusterBanksy(bank, resolution = 0.8, kneighbours = 40)
# bank <- saveRDS(bank, 'moffitt-banksy-object-animal1-processed.rds')
bank <- readRDS('moffitt-banksy-object-animal1-processed.rds')

bankODMoffitt <- SubsetBanksy(bank,
                              metadata = grepl('OD Mature 2|Immature', Cell_class))

pdf('plots/heatmap-od.pdf', height = 8, width = 10)
plotHeatmap(bankODMoffitt,
            annotate = TRUE,
            annotate.by = c('Cell_class', 'res0.8_lam0.25_k40'))
dev.off()

bank <- ClusterBanksy(bank, lambda = 0)

p1 <- plotUMAP(bank, by = 'res0.8_lam0_k40', reduction = 'umap_0')
p2 <- plotUMAP(bank, by = 'Cell_class', reduction = 'umap_0')
grid.arrange(p1, p2, ncol = 2)

pdf('plots/heatmap-moffitt-class-barplot.pdf', height = 15, width = 20)
plotHeatmap(bank, annotate = TRUE, annotate.by = c('Cell_class',
                                                   'res0.8_lam0.25_k40'),
            barplot.by = c('genes_detected', 'total_expression'),
            annotation.size = 15,
            cluster.column = TRUE)
dev.off()

param <- 'Cell_class'
s <- plotSpatialDims(bank, by = param, pt.size = 0.05) + facet_wrap( ~ Cell_class)
# u <- plotUMAP(bank, by = param, reduction = 'umap_0.25')
# hm <- plotHeatmap(bank, assay = 'banksy', lambda = 0.25,
#                   annotate = TRUE, annotate.by = param,
#                   cluster.column = TRUE)

pdf('plots/spatdim-cellclass-wrap.pdf', height = 7, width = 8)
s
dev.off()

##### Find Cell class markers --------------------------------------------------

# sClass <- plotSpatialDims(bank, by = 'Cell_class') + facet_wrap( ~ Cell_class)

cell.types <- mapvalues(bank@meta.data$Cell_class,
                        from = unique(bank@meta.data$Cell_class),
                        to  = c('AS', 'IH', 'OD', 'EC', 'UK',
                                'PE', 'EC', 'OD', 'OD', 'EX',
                                'MG', 'EC', 'OD', 'OD', 'OD', 'EP'))

seu <- CreateSeuratObject(counts = bank@own.expr)
seu <- ScaleData(seu)
Idents(seu) <- cell.types
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers <- data.table(markers)
top.markers <- markers[, head(.SD, 15), by = cluster]

pdf('plots/cell-class-markers.pdf', height = 5, width = 15)
DotPlot(seu, features = unique(top.markers$gene),
        cols = 'RdBu', cluster.idents = T) +
    theme(axis.text.x = element_text(angle = 45, size = 8,
                                     vjust = 0.8, hjust = 0.8))
dev.off()

##### Subcluster the sparse OD based on neighborhood ---------------------------

clusters <- 6
bankOD <- SubsetBanksy(bank,
                       metadata = res0.8_lam0.25_k40 %in% clusters &
                                  grepl('OD', Cell_class))

set.seed(1234)
bankOD <- ClusterBanksy(bankOD, lambda = 1, resolution = 1.2, kneighbours = 30)
table(bankOD@meta.data$res1.2_lam1_k30)

subcluster <- 'res1.2_lam1_k30'
sOD <- plotSpatialDims(bankOD, by = subcluster)
uOD <- plotUMAP(bankOD, by = subcluster, reduction = 'umap_1')
hmOD <- plotHeatmap(bankOD, assay = 'banksy', lambda = 0.25,
            annotate = TRUE, annotate.by = subcluster, cluster.column = TRUE)


##### Find NBR markers in sparse ODs -------------------------------------------

seuODnbr <- CreateSeuratObject(counts = bankOD@nbr.expr)
Idents(seuODnbr) <- bankOD@meta.data$res1.2_lam1_k30
seuODnbr <- ScaleData(seuODnbr)
markersODnbr <- FindAllMarkers(seuODnbr, only.pos = TRUE)
markersODnbr <- data.table(markersODnbr)
top.nbr.markers <- markersODnbr[p_val_adj < 0.05, head(.SD, 15), by = cluster]

pdf('plots/od-nbr-markers.pdf', height = 4, width = 12)
DotPlot(seuODnbr, features = unique(top.nbr.markers$gene),
        cols = 'RdBu', cluster.idents = T) +
    theme(axis.text.x = element_text(angle = 45, size = 8,
                                     vjust = 0.8, hjust = 0.8))
dev.off()

##### Find OWN markers in sparse ODs -------------------------------------------

seuODown <- CreateSeuratObject(counts = bankOD@own.expr)
Idents(seuODown) <- bankOD@meta.data$res1.2_lam1_k30
seuODown <- ScaleData(seuODown)
markersODown <- FindAllMarkers(seuODown, only.pos = TRUE)
markersODown <- data.table(markersODown)
# top.own.markers <- markersODown[p_val_adj < 0.05, head(.SD, 10), by = cluster]

pdf('plots/od-own-and-nbr-markers.pdf', height = 4, width = 12)
col <- rep(1, 63)
col[1:9] <- '4'
col[1:2] <- '2'
DotPlot(seuODown, features = unique(
    c(unique(markersODown$gene),
      gsub('.nbr', '', unique(top.nbr.markers$gene)))),
        cols = 'RdBu', cluster.idents = T, idents = c(1,2,3,4,5)) +
    theme(axis.text.x = element_text(angle = 45, size = 8,
                                     vjust = 0.8, hjust = 0.8,
                                     colour = col))
dev.off()

##### Match OD OWN and NBR to MOFFITT markers ----------------------------------

ref <- markers[p_val_adj < 0.05]
ref <- ref[cluster != 'OD']

## OWN markers
markersODown <- markersODown[p_val_adj < 0.05]
hit <- lapply(markersODown$gene, function(x) {
    id <- which(ref$gene == x)
    x <- ref$cluster[id]
    if (length(x) == 0) return(NA)
    return(as.character(x))
})
clust <- rep(markersODown$cluster, unlist(lapply(hit, length)))
hit <- unlist(hit)
table(clust, hit)

## NBR markers
markersODnbr$gene <- gsub('.nbr', '', markersODnbr$gene)
markersODnbr <- markersODnbr[p_val_adj < 0.05]
hit <- lapply(markersODnbr$gene, function(x) {
    id <- which(ref$gene == x)
    x <- ref$cluster[id]
    if (length(x) == 0) return(NA)
    return(as.character(x))
    })
clust <- rep(markersODnbr$cluster, unlist(lapply(hit, length)))
hit <- unlist(hit)
table(clust, hit)

##### Get neighbors ------------------------------------------------------------

cell.types <- mapvalues(bank@meta.data$Cell_class,
                        from = unique(bank@meta.data$Cell_class),
                        to  = c('AS', 'IH', 'OD', 'EC', 'UK',
                                'PE', 'EC', 'OD', 'OD', 'EX',
                                'MG', 'EC', 'OD', 'OD', 'OD', 'EP'))
dims <- c('sdimx', 'sdimy')
K <- 6
neighbors <- kNN(bank@cell.locs[,dims], K)$id
neighbors <- apply(neighbors, 2, function(x) as.character(cell.types[x]))
rownames(neighbors) <- bank@meta.data$cell_ID
OD <- rownames(neighbors) %in% bankOD@meta.data$cell_ID
neighborsOD <- neighbors[OD,]

##### CELL X CLASS fractions ---------------------------------------------------

getFractions <- function(neighbors) {
    count <- table(neighbors)
    id <- match(names(count), classes)
    fractions <- rep(0, length(classes))
    fractions[id] <- count
    fractions <- fractions / sum(fractions)
    return(fractions)
}

classes <- unique(as.vector(neighborsOD))
classes <- sort(classes)
fractions <- t(apply(neighborsOD, 1, getFractions))
colnames(fractions) <- classes

#### Aggregate mean fractions for each cluster ---------------------------------

fractions <- data.table(fractions)
fractions$cluster <- bankOD@meta.data$res1.2_lam1_k30

aggr.fractions <- fractions[, lapply(.SD, mean) , by = cluster]
aggr.mat <- as.matrix(aggr.fractions[,-1])
rownames(aggr.mat) <- aggr.fractions$cluster

pal <- Banksy:::getPalette(9)
maxlim <- max(aggr.mat)
minlim <- min(aggr.mat)
names(pal) <- classes
hmAggr <- Heatmap(aggr.mat, name = 'Fractions',
                  col = colorRamp2(c(0, maxlim/2, maxlim),
                                   c('#FFF5F0', '#FB6A4A', '#A50F15'))) %v%
    HeatmapAnnotation(Celltypes = classes, col = list(Celltypes = pal))

pdf('plots/aggregated-fractions.pdf', width = 6, height = 5)
hmAggr
dev.off()

##### Diff prop ----------------------------------------------------------------

expected <- table(cell.types)
expected <- expected / sum(expected)
aggr.mat <- sweep(aggr.mat, 2, expected, '-')

pal <- Banksy:::getPalette(9)
maxlim <- max(aggr.mat)
minlim <- min(aggr.mat)
names(pal) <- classes
hmAggr <- Heatmap(aggr.mat, name = 'Fractions',
                  col = colorRamp2(c(minlim, 0, maxlim),
                                   c("#67A9CF", 'white', '#EF8A62'))) %v%
    HeatmapAnnotation(Celltypes = classes, col = list(Celltypes = pal))

pdf('plots/aggregated-fractions-expected.pdf', width = 6, height = 5)
hmAggr
dev.off()

##### DE between tight ODs and sparse ODs --------------------------------------

od <- SubsetBanksy(bank, metadata = Cell_class == 'OD Mature 2' &
                       res0.8_lam0.25_k40 %in% c(3,6))

seu <- CreateSeuratObject(od@own.expr)
Idents(seu) <- od@meta.data$res0.8_lam0.25_k40
markersODDE <- FindAllMarkers(seu)
markersODDE <- data.table(markersODDE)
markersODDE <- markersODDE[cluster == 6]
markersODDE <- markersODDE[order(avg_log2FC, decreasing = T)]

odSparse <- SubsetBanksy(bankOD, metadata = Cell_class == 'OD Mature 2')

od@meta.data$subcluster <- od@meta.data$res0.8_lam0.25_k40
id <- match(odSparse@meta.data$cell_ID, od@meta.data$cell_ID)
od@meta.data$subcluster[id] <- paste0(od@meta.data$subcluster[id],
                                              '.',
                                              odSparse@meta.data$res1.2_lam1_k30)

pdf('plots/od-heatmap.pdf', height = 9, width = 9)
plotHeatmap(od, annotate = TRUE, annotate.by = c('subcluster', 'res0.8_lam0.25_k40'))
dev.off()

pdf('plots/od-spatdim.pdf', height = 5, width = 12)
p1 <- plotSpatialDims(od, by = 'subcluster')
p2 <- plotSpatialDims(od, by = 'res0.8_lam0.25_k40')
grid.arrange(p1, p2, ncol=2)
dev.off()

seu <- CreateSeuratObject(counts = od@own.expr)
Idents(seu) <- od@meta.data$subcluster
seu <- ScaleData(seu)
markerSubcluster <- FindMarkersBaseline(seu, baseline = '3')

pdf('plots/od-own-6.6.pdf', height = 4, width = 15)
features <- unique(
    c(
        markerSubcluster$markers$`cluster 6.6`$gene,
        unique(gsub('.nbr', '', top.nbr.markers$gene))
    )
)
col <- rep(1, length(features))
col[1:length(markerSubcluster$markers$`cluster 6.6`$gene)] <- 2
DotPlot(seu, features = features,idents = c(3, 6.6), cols = 'RdBu') +
    theme(axis.text.x = element_text(angle = 45, size = 8,
                                         vjust = 0.8, hjust = 0.8,
                                         colour = col))
dev.off()


pdf('plots/od-dense-sparse.pdf', height = 7, width = 16)
do.call(grid.arrange, c(marker$plots, ncol = 3))
dev.off()


markerSubcluster$markers

table(ref$cluster[match(markerSubcluster$markers$`cluster 6.1`$sym, ref$gene)])
table(ref$cluster[match(markerSubcluster$markers$`cluster 6.2`$sym, ref$gene)])
table(ref$cluster[match(markerSubcluster$markers$`cluster 6.3`$sym, ref$gene)])
table(ref$cluster[match(markerSubcluster$markers$`cluster 6.4`$sym, ref$gene)])
table(ref$cluster[match(markerSubcluster$markers$`cluster 6.5`$sym, ref$gene)])
table(ref$cluster[match(markerSubcluster$markers$`cluster 6.6`$sym, ref$gene)])

save(markers, markersODnbr, markersODown, markersODDE, file = 'markerList.rda')


##### Plotting specific genes

Plin3 <- od@own.expr['Mbp',od@meta.data$cell_ID]
plotData <- od@cell.locs
names(plotData)[seq_len(2)] <- c("dimx", "dimy")
plotData <- cbind(plotData,
                  data.frame(Plin3 = Plin3,
                             Cluster = od@meta.data$subcluster))

ggplot(plotData, aes(x = dimx, y = dimy, color = Plin3)) +
    xlab("x coordinates") + ylab("y coordinates") +
    theme_minimal() + geom_point(size = 2) +
    scale_color_gradient2(low = 'blue', mid = 'white', high = 'red') +
    theme(legend.text = element_text(size = 5),
          legend.title = element_blank(), axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(), panel.background = element_blank())


tmp <- readRDS('moffitt-banksy-object-animal1.rds')
tmp <- tmp@own.expr[,od@meta.data$cell_ID]
saveRDS(tmp, file = 'od-gcm.rds')
saveRDS(od@cell.locs, file = 'od-locs.rds')
saveRDS(od@meta.data, file = 'od-metadata.rds')
