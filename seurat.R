#=====================================================================================
#
#  Code chunk 0 加载包和数据，设置工作路径
#
#=====================================================================================

library(dplyr)
library(Seurat)
library(patchwork)

setwd("E:/scp259.cssmilie/script")
source('analysis.r')


setwd("E:/scp259.cssmilie/")
cell_subsets = read.table('input/cell_subsets.txt', sep='\t', header=F, stringsAsFactors=F)
meta = read.table('input/all.meta2.txt', sep='\t', header=T, row.names=1, stringsAsFactors=F)
epi.batch = setNames(meta[colnames(epi.counts), 'Sample'], colnames(epi.counts))


epi.counts = readMM('input/gene_sorted-Epi.matrix.mtx')
rownames(epi.counts) = readLines('input/Epi.genes.tsv')
colnames(epi.counts) = readLines('input/Epi.barcodes2.tsv')

epi.batch = setNames(meta[colnames(epi.counts), 'Sample'], colnames(epi.counts))

scp259.seu <- CreateSeuratObject(
    counts = epi.counts, 
    project = "scp259", 
    min.cells = 3, 
    min.features = 200
)
rm(epi.counts)
gc()
saveRDS(scp259.seu, "data/scp259.seu0.RDS")
save(cell_subsets, meta, epi.batch, file = "data/meta0.RData")

#=====================================================================================
#
#  Code chunk 0 加入metadata
#
#=====================================================================================

rm(list = ls())
gc()
scp259.seu <- readRDS("data/scp259.seu.RDS")
load("data/meta0.RData")

# meta <- subset(meta, meta$Location=="Epi")
meta <- meta[-1,]
meta <- meta[,-c(1,6,7)]
colnames(meta)[3] <- "geo.accession"
colnames(meta)[4] <- "inflammation"

meta[which(meta[4]=="Healthy"),4] <- "Control"

meta <- data.frame(
    meta,
    meta[4])
colnames(meta)[5] <- "subject status"
meta[which(meta[5]=="Healthy"),5] <- "Healthy Control"
meta[which(meta[5]=="Non-inflamed"),5] <- "Ulcerative Colitis"
meta[which(meta[5]=="Inflamed"),5] <- "Ulcerative Colitis"
meta <- data.frame(
    meta,
    "epithelial cells")
colnames(meta)[6] <- "cell.type"

scp259.seu <- AddMetaData(object = scp259.seu, metadata = meta)

saveRDS(scp259.seu, "data/scp259.seu.RDS")

#=====================================================================================
#
#  Code chunk 0 数据清洗
#
#=====================================================================================

scp259.seu <- readRDS("data/scp259.seu.RDS")

# [["orig.ident"]] <- scp259.seu[["geo.accession"]]

scp259.seu[["percent.mt"]] <- PercentageFeatureSet(scp259.seu, pattern = "^MT-")
scp259.seu[["percent.rb"]] <- PercentageFeatureSet(scp259.seu, pattern = "^RP[SL]")

VlnPlot(
    scp259.seu, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
    #    split.by = "geo.accession",
    pt.size = 0,
    ncol = 4
)


scp259.seu <- subset(
    scp259.seu, 
    subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 75 & percent.rb < 50 & nCount_RNA < 40000
)
scp259.seu1

# saveRDS(scp259.seu, "data/scp259.seu.beforeclean.RDS")


VlnPlot(
    scp259.seu, 
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rb"), 
    pt.size = 0,
    ncol = 4
)
saveRDS(scp259.seu, "data/scp259.seu1.RDS")







scp259.seu <- NormalizeData(scp259.seu)

scp259.seu <- FindVariableFeatures(
    scp259.seu, selection.method = "vst", 
    nfeatures = 2000
)


scp259.seu <- ScaleData(scp259.seu)

scp259.seu <- RunPCA(
    scp259.seu, 
    features = scp259.seu@var.genes
    #VariableFeatures(object = scp259.seu)
)
gc()

scp259.har <- RunHarmony(scp259.seu, group.by.vars = "geo.accession")


#降维聚类
scp259.har <- RunUMAP(scp259.har, reduction = "harmony", dims = 1:30)
scp259.har <- FindNeighbors(scp259.har, reduction = "harmony", dims = 1:30) %>% FindClusters(resolution = 0.3)
##作图
#group_by_cluster
plot1 = DimPlot(scp259.har, reduction = "umap", label=T) 
#group_by_sample
plot2 = 
    DimPlot(scp259.har, reduction = "umap", group.by='geo.accession', split.by = 'geo.accession', ncol = 6) 
#combinate
plot1+plot2


plot2


saveRDS(scp259.har, "data/scp259.har.RDS")




DimPlot(scp259.har, reduction = "umap", label=T) 

scp259.har <- RunTSNE(scp259.har, reduction = "harmony", dims = 1:30)

DimPlot(scp259.har, reduction = "tsne", label=T) 



marker.all <- FindAllMarkers(
    scp259.har, 
    min.pct = 0.25
)
gc()

# 先找到EEC是哪一群
VlnPlot(
    scp259.har, 
    features = c("EPCAM", "KRT8", "KRT18", "TPH1", "CHGA", "CHGB", "PCSK1N"),
    pt.size = 0,
    ncol = 2
)
FeaturePlot(
    scp259.har, 
    features = c("EPCAM", "KRT8", "KRT18", "TPH1", "CHGA", "CHGB", "PCSK1N")
)
DotPlot(
    scp259.har, 
    features = c("EPCAM", "KRT8", "KRT18", "TPH1", "CHGA", "CHGB", "PCSK1N")
)


scp259.har1 <- RenameIdents(
    scp259.har,
    `22` = "EEC"
)


# 识别保守marker
marker.ec <- FindConservedMarkers(
    scp259.har, 
    ident.1 = 22,
    grouping.var = "subject.status", 
    verbose = FALSE
)


#==============================================================================
#
#  Code chunk 5 EEC细分群，把EC挑出来
#
#==============================================================================

eec.har <- readRDS("data/eec.har.RDS")

eec.har1 <- NormalizeData(eec.har)

eec.har1 <- FindVariableFeatures(
    eec.har1, selection.method = "vst", 
    nfeatures = 2000
)

eec.har1 <- ScaleData(eec.har1)
eec.har1 <- RunPCA(
    eec.har1, 
    features = VariableFeatures(object = eec.har1)
)
eec.har1 <- RunHarmony(eec.har1, group.by.vars = "geo.accession")


ElbowPlot(eec.har1)

eec.har1 <- RunUMAP(eec.har1, reduction = "harmony", dims = 1:13)

eec.har1 <- FindNeighbors(eec.har1, reduction = "harmony", dims = 1:13) 
eec.har1 <- FindClusters(eec.har1, resolution = 0.2)

CairoPDF(
    "output/harmony后整体eec umap2.pdf", 6, 5)
# 5.5*5
DimPlot(
    eec.har1, 
    reduction = "umap", 
    label=TRUE,
    repel=TRUE
)
dev.off()



FeaturePlot(
    eec.har1, 
    features = c("TPH1"),
    label = F
)


eec.har1
# An object of class Seurat 
# 20028 features across 604 samples within 1 assay 
# Active assay: RNA (20028 features, 2000 variable features)
# 4 dimensional reductions calculated: pca, harmony, umap, tsne


table(eec.har1@active.ident)
# 0   1   2 
# 271 167 166 






ggsave(
    filename = "harmony后后EEC的umap.png",
    path = "output/",
    width = 2000  ,
    height = 1500  ,
    units = "px"
)



eec.har1 <- RenameIdents(
    eec.har1,
    `0` = "EC",
    `1` = "EEC1",
    `2` = "EEC2"
)


saveRDS(eec.har1, "data/eec.har2.RDS")

eec.har1 <- readRDS("data/eec.har2.RDS")
# EEC marker
CairoPDF(
    "output/eec marker featureplot.pdf", 10, 8)
FeaturePlot(
    eec.har1, 
    features = c("TPH1", "CHGA", "CHGB", "PYY", "GCG", "SST",  "SCT", "PCSK1N", "SCG5"),
    combine = TRUE,
    label = F
) +
    theme(legend.position = "none",
          axis.title.y = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()
          
    )  







dev.off()





scp259.har1 <- subset(
    scp259.har,
    subset = seurat_clusters == "1" |
        seurat_clusters == "2" |
        seurat_clusters == "3" |
        seurat_clusters == "4" |
        seurat_clusters == "5" |
        seurat_clusters == "6" |
        seurat_clusters == "7" |
        seurat_clusters == "8" |
        seurat_clusters == "9" |
        seurat_clusters == "10" |
        seurat_clusters == "11" |
        seurat_clusters == "12" |
        seurat_clusters == "13" |
        seurat_clusters == "14" |
        seurat_clusters == "15" |
        seurat_clusters == "16" |
        seurat_clusters == "17" |
        seurat_clusters == "18" |
        seurat_clusters == "19" |
        seurat_clusters == "20" |
        seurat_clusters == "21" |
        seurat_clusters == "23" |
        seurat_clusters == "24" |
        seurat_clusters == "25" )


gc()

scp259.har1 <- merge(eec.har1, scp259.har1, idents = TRUE)

saveRDS(scp259.har1, "data/scp259.har.ecdefined.RDS")


table(scp259.har1@active.ident)


ec0 <- readRDS("data/ec.RDS")

DimPlot(ec, reduction = "umap", split.by = "subject.status1")


mtx.ec <- GetAssayData(ec, slot = "counts")

# 删除稀有基因，< 10%
keep_gene = rowSums(mtx.ec != 0) > ncol(mtx.ec)/10
mtx.ec1 = mtx.ec[keep_gene, ]
mtx.ec1

ec <- CreateSeuratObject(mtx.ec1)
ec <- AddMetaData(ec, ec0@meta.data)
Idents(ec) <- "subject.status1"


inf2con <- FindMarkers(
    ec, 
    ident.1 = "EC_Inflamed", 
    ident.2 = "EC_Control", 
    verbose = FALSE
)
unif2con <- FindMarkers(
    ec, 
    ident.1 = "EC_Non-inflamed", 
    ident.2 = "EC_Control", 
    verbose = FALSE
)
inf2unif <- FindMarkers(
    ec, 
    ident.1 = "EC_Inflamed", 
    ident.2 = "EC_Non-inflamed", 
    verbose = FALSE
)





















