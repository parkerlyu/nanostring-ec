#=====================================================================================
#
#  Code chunk 1 加载包、路径和数据
#
#=====================================================================================

rm(list=ls())

library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(UpSetR)
library(pheatmap)
library(Cairo)




#=====================================================================================
#
#  Code chunk 2 差异基因分析1（UC EC vs Con EC）
#
#=====================================================================================

sample <- SAMPLE[,-c(6,7,14,15)]
colsample <- colSAMPLE[-c(6,7,14,15),1:2]


dds <- DESeqDataSetFromMatrix(
    sample, 
    colsample,
    design= ~ subject
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("subject", "UC", "Con"))

res$padj[is.na(res$padj)] <- 1
res <- res[order(res$log2FoldChange),]

write.table(res, file="UCEC2ConEC.res.csv", row.names = T, col.names = T, quote = FALSE, sep=',')

## 提取差异表达基因
deseq2.sig<- data.frame(subset(
    res, 
    padj < 0.05
    & abs(log2FoldChange) > 1
))

write.table(deseq2.sig, file = "UCEC2ConEC.sig.csv", row.names = T, col.names = T, quote = FALSE,sep=',') 

#=====================================================================================
#
#  Code chunk 3 绘制差异基因火山图1（UC EC vs Con EC）
#
#=====================================================================================

res <- data.frame(res)

res$threshold <- 
    factor(
        ifelse(res$padj < 0.05 & abs(res$log2FoldChange) >= 1, 
               ifelse(res$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig$threshold <- 
    factor(
        ifelse(deseq2.sig$padj < 0.05 & abs(deseq2.sig$log2FoldChange) >= 1, 
               ifelse(deseq2.sig$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig <- data.frame(deseq2.sig)
# res$padj <- res$padj + 1*10^(-5)
# deseq2.sig$padj <- deseq2.sig$padj + 1*10^(-5)

CairoPDF(file = "VP.UCEC2ConEC.pdf",3.8,3.8)
ggplot(
    res,
    aes(
        x = log2FoldChange,
        y = -log10(padj)
        ,
        color = threshold
    ) 
)+ 
    geom_point(size = 0.6, pch = 20)+
    scale_color_manual(values=c("#CD2626",
                                #"#000080",
                                "#808080"))+
    ylab('-log10padj')+
    xlab('log2FoldChange')+
    theme_bw()+
    theme(
        legend.position="none",
        legend.title = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)
    )+
    geom_vline(xintercept=-1, linetype=2, colour="gray30") +
    geom_vline(xintercept=1, linetype=2, colour="gray30") +
    geom_hline(yintercept=-log10(0.05), linetype=2, colour="gray30") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    #  scale_x_continuous(limit = c(-5.1,5.1)) + 
    
    geom_text_repel(
        data = deseq2.sig,
        aes(
            x = log2FoldChange,
            y = -log10(padj),
            label = rownames(deseq2.sig)),
        size = 2.5,
        segment.color = "black", 
        show.legend = T,
        min.segment.length = 0.3,
        max.overlaps= 10,
        force = 3
    )
dev.off()


#=====================================================================================
#
#  Code chunk 4 差异基因分析2（UC EC vs UC Epi）
#
#=====================================================================================

sample <- SAMPLE[,8:15]
colsample <- colSAMPLE[8:15,1:2]



dds <- DESeqDataSetFromMatrix(
    sample, 
    colsample,
    design= ~ cell
)
dds <- DESeq(dds)
res1 <- results(dds, contrast = c("cell", "EC", "Epithelium"))

## padj为NA的替换为1
res1$padj[is.na(res1$padj)] <- 1
res1 <- res1[order(res1$log2FoldChange),]

write.table(res1, file="UCEC2UCEpi.res.csv", row.names = T, col.names = T, quote = FALSE, sep=',')


## 提取差异表达基因，padj<0.1！！！
deseq2.sig1<- data.frame(subset(
    res1, 
    padj < 0.05
    & abs(log2FoldChange) > 1
))

write.table(deseq2.sig1, file = "UCEC2UCEpi.sig.csv", row.names = T, col.names = T, quote = FALSE,sep=',') 


#=====================================================================================
#
#  Code chunk 5 绘制差异基因火山图2（UC EC vs UC Epi）
#
#=====================================================================================

res1 <- data.frame(res1)

res1$threshold <- 
    factor(
        ifelse(res1$padj < 0.05 & abs(res1$log2FoldChange) >= 1, 
               ifelse(res1$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig1$threshold <- 
    factor(
        ifelse(deseq2.sig1$padj < 0.05 & abs(deseq2.sig1$log2FoldChange) >= 1, 
               ifelse(deseq2.sig1$log2FoldChange>= 1 ,'Up','Down'),'No Significance'),
        levels=c('Up','Down','No Significance')
    )

deseq2.sig1 <- data.frame(deseq2.sig1)
# res1$padj <- res1$padj + 1*10^(-5)
# deseq2.sig1$padj <- deseq2.sig1$padj + 1*10^(-5)

CairoPDF(file = "VP.UCEC2UCEpi.pdf",3.8,3.8)
ggplot(
    res1,
    aes(
        x = log2FoldChange,
        y = -log10(padj)
        ,
        color = threshold
    ) 
)+ 
    geom_point(size = 0.6, pch = 20)+
    scale_color_manual(values=c("#CD2626","#000080","#808080"))+
    ylab('-log10padj')+
    xlab('log2FoldChange')+
    theme_bw()+
    theme(
        legend.position="none",
        legend.title = element_blank(),
        axis.text=element_text(size=10),
        axis.title=element_text(size=10)
    )+
    geom_vline(xintercept=-1, linetype=2, colour="gray30") +
    geom_vline(xintercept=1, linetype=2, colour="gray30") +
    geom_hline(yintercept=-log10(0.05), linetype=2, colour="gray30") +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank()
    ) +
    scale_x_continuous(limit = c(-6,6)) + 
    
    geom_text_repel(
        data = deseq2.sig1,
        aes(
            x = log2FoldChange,
            y = -log10(padj),
            label = rownames(deseq2.sig1)),
        size = 2.5,
        segment.color = "black", 
        show.legend = T,
        min.segment.length = 0.5,
        max.overlaps= 10,
        force = 10
    )
dev.off()


#=====================================================================================
#
#  Code chunk 6 保存结果用于富集分析
#
#=====================================================================================

save(res0, res, res1, deseq2.sig0, deseq2.sig, deseq2.sig1, file = "DEG_all.RData")


#=====================================================================================
#
#  Code chunk 7 UpsetR
#
#=====================================================================================

sig0.up <- rownames(subset(deseq2.sig0, threshold == "Up"))
sig0.down <- rownames(subset(deseq2.sig0, threshold == "Down"))

sig.up <- rownames(subset(deseq2.sig, threshold == "Up"))
sig.down <- rownames(subset(deseq2.sig, threshold == "Down"))

sig1.up <- rownames(subset(deseq2.sig1, threshold == "Up"))
sig1.down <- rownames(subset(deseq2.sig1, threshold == "Down"))


dataForUpSetPlot <- list(
    #  UC.EC_vs_Con.EC1=sig.down,  
    UC.EC_vs_UC.Epi1=sig1.down,
    Con.EC_vs_Con.Epi1=sig0.down,
    
    UC.EC_vs_Con.EC=sig.up,
    UC.EC_vs_UC.Epi=sig1.up,
    Con.EC_vs_Con.Epi=sig0.up  
)

#setsBarColors <-c('#EA4335', '#FBBC05', '#34A853', '#4285F4')

CairoPDF("USP_updown.pdf",10,5)
upset(
    fromList(dataForUpSetPlot),
    nsets=length(dataForUpSetPlot),
    sets = names(dataForUpSetPlot),
    keep.order = TRUE,
    point.size = 3,
    line.size = 1,
    number.angles = 0,
    text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 2), # ytitle, ylabel, xtitle, xlabel, sets, number
    order.by="degree",
    matrix.color="black",
    main.bar.color = 'black',
    mainbar.y.label = 'Number of Intersection Genes',
)
dev.off()

save.image("all.RData")


#=====================================================================================
#
#  Code chunk 7 绘制差异基因热图（总）
#
#=====================================================================================

load("DEG.RData")

SAMPLE <- read.csv(
    "E:/nanostring.remaster/Data_Input/RawCounts_Rearranged_Clean_LI.csv",
    header=T,row.names=1)

dds <- DESeqDataSetFromMatrix(
    SAMPLE, 
    colSAMPLE,
    design= ~ subject
)
dds <- DESeq(dds)

normSAMPLE <- counts(dds, normalized=TRUE)
normSAMPLE.mad <- apply(normSAMPLE, 1, mad)
normSAMPLE <- normSAMPLE[order(normSAMPLE.mad, decreasing=T), ]









SAMPLE.sig.all <- unique(c(rownames(deseq2.sig0), rownames(deseq2.sig), rownames(deseq2.sig1)))




SAMPLE.sig.all <- unique(c(rownames(deseq2.sig0), rownames(deseq2.sig)))
SAMPLE.sig.all <- normSAMPLE[SAMPLE.sig.all,]

annotation_col <- colSAMPLE[,c("cell","subject")]

ann_colors = list(
    colSAMPLE$cell,
    colSAMPLE$subject)


CairoPDF("HM_EC_all.pdf", 8, 8)

pheatmap(
    SAMPLE.sig.all,
    scale = "row" , 
    cluster_cols = F,
    cluster_rows = T,
    clustering_distance_rows = "correlation",
    
    show_rownames=T,
    show_colnames=F,
    
    border = F,
    
    color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    
    treeheight_row = 15, 
    treeheight_col = 10,
    
    legend = T,
    
    annotation_col = annotation_col,
    annotation_legend = TRUE,
    annotation_colors = ann_colors,
    
    fontsize = 6,
    cellwidth = 15
    #,
    #cellheight = 6
) 

dev.off()










