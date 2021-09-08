
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
library(qvalue)
library(schelpers)

###### Find DE between sparse and dense OD -------------------------------------

bank <- readRDS('moffitt-banksy-object-animal1-processed.rds')
bank <- SubsetBanksy(bank, metadata = grepl('OD Mature 2', Cell_class) & 
                                      res0.8_lam0.25_k40 %in% c(3,6))
bank@meta.data$OD_cluster <- mapvalues(bank@meta.data$res0.8_lam0.25_k40,
                                       from = c(3,6), to = c('Dense', 'Sparse'))
plotSpatialDims(bank, by = 'OD_cluster') 

seu <- CreateSeuratObject(counts = bank@own.expr)
Idents(seu) <- bank@meta.data$OD_cluster
markers <- FindAllMarkers(seu, only.pos = TRUE)
markers <- data.table(markers)
top.markers <- markers[p_val_adj < 0.01]
od.markers <- top.markers$gene
    
cor.genes <- lapply(od.markers, function(gene) {
    
    message(gene)
    cor <- rep(NA, nrow(bank@own.expr))
    names(cor) <- rownames(bank@own.expr)
    for (i in 1:nrow(bank@own.expr)) {
        if (grepl('^Blank', rownames(bank@own.expr)[i])) {
            next
        }
        cor[i] <- cor(bank@own.expr[gene,],
                      bank@own.expr[i,])
    }
    cor <- sort(cor, decreasing = TRUE)[2:11]
    return(names(cor))
})
od.set <- unique(unlist(cor.genes), od.markers)
saveRDS(od.set, file = 'od-markers.rds')
od.set <- readRDS('od-markers.rds')

###### Extract ODs -------------------------------------------------------------

# counts <- Read10X('scrnaseq/')
# saveRDS(counts, file = 'scrnaseq/counts.rds')
# counts <- readRDS('scrnaseq/counts.rds')

# seu <- CreateSeuratObject(counts = counts)
# saveRDS(seu, file = 'scrnaseq/seu.rds')
# seu <- readRDS('scrnaseq/seu.rds')
# 
# meta <- openxlsx::read.xlsx('tables/aau5324_moffitt_table-s1.xlsx')
# names(meta) <- c('Cell_name', 'Sex', 'Rep', 'Cell_class', 
#                  'Non_neuronal_cluster', 'Neuronal_cluster')
# seu@meta.data <- cbind(seu@meta.data, meta[match(meta$Cell_name, rownames(seu@meta.data)),])
# Idents(seu) <- seu@meta.data$Cell_class
# seu <- subset(seu, idents = c('Immature oligodendrocyte',
#                               'Mature oligodendrocyte',
#                               'Newly formed oligodendrocyte'))
# saveRDS(seu, 'scrnaseq/seu-od.rds')
# seu <- subset(seu, idents = c('Mature oligodendrocyte'))
# saveRDS(seu, 'scrnaseq/seu-odm.rds')

###### All ODs -----------------------------------------------------------------

seu <- readRDS('scrnaseq/seu-od.rds')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu)

seu <- RunTSNE(seu, dims = 1:20, perplexity = 100)
tsne <- DimPlot(seu, reduction = 'tsne', group.by = 'Cell_class', cols = scpalette())
seu <- RunUMAP(seu, dims = 1:20, min.dist = 0.3)
umap <- DimPlot(seu, reduction = 'umap', group.by = 'Cell_class', cols = scpalette())

pdf('plots/rna-od-red.pdf', height = 6, width = 15)
grid.arrange(umap, tsne, ncol = 2)
dev.off()

seu <- FindNeighbors(seu)
seu <- FindClusters(seu, resolution = 0.1)
DimPlot(seu, reduction = 'umap', group.by = 'seurat_clusters', cols = scpalette())
DimPlot(seu, reduction = 'tsne', group.by = 'seurat_clusters', cols = scpalette())

markers <- FindAllMarkers(seu)
top.markers <- data.table(markers)
top.markers <- top.markers[, head(.SD, 10), by = cluster]

pdf('plots/rna-feature.pdf', height = 25, width = 26)
FeaturePlot(seu, reduction = 'tsne', features = c('Neu1', 'Pdgfa', 'Cspg4',
                                                  'Foxo4', 'Galc',
                                                  'Mbp', 'Mag', 'Mog', 'Omg', 'Cnp',
                                                  'S100b', 'Plin3', 'Lpar1'))
dev.off()

pdf('plots/rna-feature-plin3.pdf', height = 25, width = 26)
plin3 <- readRDS('plin3-markers.rds')
FeaturePlot(seu, reduction = 'tsne', features = plin3)
dev.off()

mat <- GetAssayData(seu, 'scale.data')
plin3 <- mat[which(rownames(mat) == 'Plin3'),]
cor <- rep(NA, nrow(mat))
names(cor) <- rownames(mat)
for (i in 1:nrow(mat)) {
    if (i %% 10 == 0) print(paste0(round(i/nrow(mat)*100,2), ' % complete'))
    cor[i] <- cor(plin3, mat[i,])
}
cor <- sort(cor, decreasing = TRUE)

###### Mature ODs --------------------------------------------------------------

seu <- readRDS('scrnaseq/seu-odm.rds')
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)
ElbowPlot(seu)

seu <- RunTSNE(seu, dims = 1:20, perplexity = 100)
DimPlot(seu, reduction = 'tsne', group.by = 'Non_neuronal_cluster', cols = scpalette())
seu <- RunUMAP(seu, dims = 1:20, min.dist = 0.3)
DimPlot(seu, reduction = 'umap', group.by = 'Non_neuronal_cluster', cols = scpalette())

seu <- FindNeighbors(seu)
seu <- FindClusters(seu, resolution = 0.2)
umap <- DimPlot(seu, reduction = 'umap', group.by = 'seurat_clusters', cols = scpalette())
tsne <- DimPlot(seu, reduction = 'tsne', group.by = 'seurat_clusters', cols = scpalette())

pdf('plots/rna-odm-red.pdf', height = 6, width = 12)
grid.arrange(umap, tsne, ncol = 2)
dev.off()

markers <- FindAllMarkers(seu, only.pos = TRUE)
markers.filtered <- data.table(markers)[p_val_adj < 0.01]
top.markers <- markers.filtered[, head(.SD, 10), by = cluster]

pdf('plots/rna-feature-odm.pdf', height = 20, width = 25)
feat <- c('Mbp', 'Lpar1', 'Mag', 'Mog',
          'S100b', 'Klk6', 'Sepp1', 'Plin3', 'Pmp22', 'Enpp6')
FeaturePlot(seu, reduction = 'umap', features = feat)
dev.off()

markers.plin3 <- markers.filtered[cluster == 2][pct.1 - pct.2 > 0.35]
markers.plin3
saveRDS(markers.plin3$gene, 'plin3-markers.rds')

##### Sub-cluster scrna-seq manually -------------------------------------------

mat <- as.matrix(GetAssayData(seusub, 'counts'))
mat <- mat[od.set[od.set %in% rownames(mat)], ]

rm(seusub)

mat <- Banksy:::normalizer(mat, 100)
mat <- Banksy:::scaler(mat)

## Run PCA
pca <- prcomp(t(mat))
plot(pca$x[,1], pca$x[,2], pch = 19, col = rgb(0,0,0,0.2), cex = 0.5)

## Elbow plot
vv <- pca$sdev^2 / sum(pca$sdev^2) * 100
ndims <- 50
plot(1:ndims, vv[1:ndims], pch = 19, cex = 0.5)

ndims <- 20

## UMAP with pca
umap <- uwot::umap(pca$x[,1:ndims], init = 'spectral', min_dist = 0.2)
plot(umap[,1],umap[,2], pch = 19, col = rgb(0,0,0,0.2), cex = 0.5)

## Densmap
densmap <- densvis::densmap(pca$x[,1:ndims], dens_frac = 0.5, dens_lambda = 0.5, min_dist = 0.5)
plot(densmap[,1], densmap[,2], pch = 19, col = rgb(0,0,0,0.2), cex = 0.5)

## consolidate
data <- data.frame(PC1 = pca$x[,1],
                   PC2 = pca$x[,2],
                   UMAP1 = umap[,1],
                   UMAP2 = umap[,2],
                   DMAP1 = densmap[,1],
                   DMAP2 = densmap[,2],
                   t(mat))

ggplot(data, aes(x=UMAP1, y=UMAP2, col = Mbp)) +
    geom_point(size = 1, alpha = 0.9) + theme_minimal() + 
    theme(panel.background = element_rect(fill = 'black'),
          plot.background = element_rect(fill = 'black'),
          legend.text = element_text(color = 'white'),
          axis.text = element_text(color = 'white'),
          axis.line = element_line(color = 'white'),
          panel.grid = element_blank(),
          plot.title = element_text(color = 'white', hjust = 0.5)) +
    scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
    ggtitle('Mbp')

genePlot <- function(gene) {
    
    p <- ggplot(data, aes(x=UMAP1, y=UMAP2, col = .data[[gene]])) +
        geom_point(size = 1, alpha = 0.9) + theme_minimal() + 
        theme(panel.background = element_rect(fill = 'black'),
              plot.background = element_rect(fill = 'black'),
              legend.text = element_text(color = 'white'),
              axis.text = element_text(color = 'white'),
              axis.line = element_line(color = 'white'),
              panel.grid = element_blank(),
              plot.title = element_text(color = 'white', hjust = 0.5)) +
        scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0) +
        ggtitle(gene)
    
    return(p)
}

p1<-genePlot('Mbp')
p2<-genePlot('Lpar1')
p3<-genePlot('Plin3')
grid.arrange(p1, p3, ncol=2)

gene = 'Lpar1'
ggplot(data, aes(x=UMAP1, y=UMAP2, col = .data[[gene]])) +
    geom_point(size = 1, alpha = 0.9) + theme_minimal() + 
    theme(panel.background = element_rect(fill = 'black'),
          plot.background = element_rect(fill = 'black'),
          legend.text = element_text(color = 'white'),
          axis.text = element_text(color = 'white'),
          axis.line = element_line(color = 'white'),
          panel.grid = element_blank(),
          plot.title = element_text(color = 'white', hjust = 0.5)) +
    scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 2) +
    ggtitle(gene)
