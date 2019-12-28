#!/usr/bin/env R

# 20180723
# Nocturnin: GST or NOCT (del 2-15) overexpression RNA-seq
# aligned by STAR, quantified by HTSeq
# with customized genome and annotation (ref human + plasmid + ERCC)

# rm(list = ls())

deseqDirName <- '~/data/noct/humanOE/05-deseq2_check_results_20191219/'
source('~/projects/styles/R/general.R')

##### MANUAL HTSEQ TABLE BINDING #####
htseqDirName <- '~/data/noct/humanOE/04-htseq'
htseqFileNameSet <- list.files(path = htseqDirName, pattern = '.htseq.count', all.files = FALSE, full.names = TRUE)
print(htseqFileNameSet)

stripSampleName <- function(filename){
  sampleName <- strsplit(x = basename(filename), split = '.htseq.count', fixed = TRUE)[[1]][1]
  return(sampleName)
}

stripSampleCondition <- function(filename){
  sampleCondition <- strsplit(x = basename(filename), split = '_', fixed = TRUE)[[1]][1]
  if(sampleCondition == 'Del'){
    return('NOCT')
  } else {
    return(sampleCondition)    
  }

}

##### VIGNETTE #####
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#htseq

library('DESeq2')

########################################################
# what if I source the patch for lfcSE after importing #
########################################################
# would it replace the orignial functon?
source('~/projects/noct/humanOE/DESeq2.1.61.1.patch.R')

###### htseq-count input ######

directory <- htseqDirName

sampleFileNameSet <- htseqFileNameSet
sampleNameSet <- sapply(X = sampleFileNameSet, FUN = stripSampleName)
sampleConditionSet <- sapply(X = sampleFileNameSet, FUN = stripSampleCondition)

sampleTable <- data.frame(sampleName = sampleNameSet,
                          fileName = basename(sampleFileNameSet),
                          condition = sampleConditionSet)
# build the DESeqDataSet 
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = directory,
                                       design= ~ condition)
ddsHTSeq
# > ddsHTSeq
# class: DESeqDataSet 
# dim: 58489 6 
# metadata(1): version
# assays(1): counts
# rownames(58489): ENSG00000000003 ENSG00000000005 ... pFC3F_GST pFC3F_NOCT_del2-15
# rowData names(0):
#   colnames(6): Del_3F_11_S5.stranded Del_3F_12_S6.stranded ... GST_3F_8_S2.stranded
# GST_3F_9_S3.stranded
# colData names(1): condition
print(head(ddsHTSeq))
head(ddsHTSeq@assays$data$counts)
rawCountTbl <- ddsHTSeq@assays$data$counts
print(head(rownames(ddsHTSeq)))
# write.table(x = rawCountTbl, 
#             file = getOutFilename(name = 'deseq_binded_raw_counts', saveTo = deseqDirName, 
#                                   proj = 'noct', suffix = 'csv'), 
#             col.names = NA,row.names = TRUE, sep = ',', quote = TRUE)

# to match the variable names used in vignette
dds <- ddsHTSeq

###### pre-filtering ######
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]
dds
# from numbers of genes (see "rowname" counts), low counts filter 58489 to 23417
# > dds
# class: DESeqDataSet 
# dim: 23417 6 
# metadata(1): version
# assays(1): counts
# rownames(23417): ENSG00000000003 ENSG00000000419 ... pFC3F_GST pFC3F_NOCT_del2-15
# rowData names(0):
#   colnames(6): Del_3F_11_S5.stranded Del_3F_12_S6.stranded ... GST_3F_8_S2.stranded
# GST_3F_9_S3.stranded
# colData names(1): condition

###### set reference level ######
print(dds$condition)
# dds$condition <- factor(dds$condition, levels = c('GST', 'NOCT'))
dds$condition <- relevel(dds$condition, ref = 'GST')
print(dds$condition)

###### collapse TECHNICAL replicates ######
# do not collapse biological replicates using function `collapseReplicates()`

###### differential expression analysis ######
dds <- DESeq(dds)
res <- results(dds)
print(res)
# specified the coefficient or contrast IF was not done in 'setting reference level'

###### log fold change shrinkage for visualization and ranking ######
resultsNames(dds)
# > resultsNames(dds)
# [1] "Intercept"            "condition_Del_vs_GST"
# DESeq2 v 1.20
# resLFC <- lfcShrink(dds, coef = 'condition_Del_vs_GST', type = 'normal')
# DESeq2 v 1.16
resLFC <- lfcShrink_patched(dds, contrast = c('condition', 'NOCT', 'GST'), type = 'normal', res = res)
write.table(x = resLFC, 
            file = getOutFilename(name = 'patched.new_deseq2_shrinkage', proj = 'noct',
                                  saveTo = deseqDirName, suffix = '.csv'), 
            col.names = NA, row.names = TRUE, sep = ',', quote = TRUE)

###### p-values and adjusted p-values ###### 

resLfcOrdered <- resLFC[order(resLFC$pvalue), ]
summary(resLfcOrdered)
print(head(resLfcOrdered))
print(tail(resLfcOrdered))
write.table(x = resLfcOrdered, 
            file = getOutFilename(name = 'patched.deseq2_lfcShrink_ordered_table', proj = 'noct',
                                  saveTo = deseqDirName, suffix = '.csv'),
            col.names = NA, row.names = TRUE, quote = TRUE, sep = ',')
# out of 23417 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 454, 1.9% 
# LFC < 0 (down)   : 368, 1.6% 
# outliers [1]     : 3, 0.013% 
# low counts [2]   : 6356, 27% 
# (mean count < 10)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

# ###### independent hypothesis weighting ######
# library('IHW')
# resIHW <- results(dds, filterFun = ihw)
# summary(resIHW)
# write.table(x = as.data.frame(resIHW), 
#             file = paste(deseqDirName, '20190315-noct-deseq2_ihw_table.csv', sep = ''), 
#             col.names = NA, row.names = TRUE, quote = TRUE, sep = ',')
# # out of 23417 with nonzero total read count
# # adjusted p-value < 0.1
# # LFC > 0 (up)     : 467, 2% 
# # LFC < 0 (down)   : 386, 1.6% 
# # outliers [1]     : 3, 0.013% 
# # [1] see 'cooksCutoff' argument of ?results
# # [2] see metadata(res)$ihwResult on hypothesis weighting
# 
# sum(resIHW$padj < 0.1, na.rm=TRUE)
# sum(resIHW$padj < 0.05, na.rm=TRUE)
# metadata(resIHW)$ihwResult

###### exploring and exporting results ######

# it all uses results from shrunken LFC test from this point on

library('ggplot2')

plotSavePath <- '~/data/noct/humanOE/05-deseq2_rerun_20190315/'

####### volcano plot of sig genes ####### 

plotData <- as.data.frame(resLFC)
# alphaThreshold <- 0.1
divideLevels <- function(x){
  if (abs(x) > 2){
    return(2)
  }else if (abs(x) > 1){
    return(1)
  }else{
    return(0)
  }
}
plotData$highChange <- factor(sapply(X = plotData$log2FoldChange, FUN = divideLevels))
print(head(plotData$highChange))
plotData$negLog10Adjp <- -log10(plotData$padj)
print(head(plotData))

shrunkenVolcano <- ggplot(data = plotData, aes(x = log2FoldChange, y = negLog10Adjp, color = highChange)) + 
  geom_point(size = 0.75) +
  geom_hline(aes(yintercept = -log10(0.05), linetype = '-log10(0.05)'), color = 'blue') + 
  scale_x_continuous(breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(breaks = seq(0, 55, 5)) + 
  scale_color_manual(values = c('grey50', 'black', 'red'),
                     breaks = c(0, 1, 2),
                     labels = c('<= 2 fold', '2 fold < change <= 4 fold', '> 4 fold')) +
  scale_linetype_manual(name = 'Significance', values = c('-log10(0.05)' = 'dashed')) + 
  labs(title = 'Differential Gene Expression with NOCT ∆2-15 Overexpression (NOCT/GST)',
       subtitle = 'Data analyzed by STAR, HTSeq, and DESeq2.',
       x = 'Shrunken Log2 Fold Change',
       y = 'Adjusted p-value (log 10)',
       caption = 'NA in adjusted p-values by DESeq2 removed',
       color = 'Fold change levels: ',
       linetype = 'Significance:') +
  theme_bw() +
  theme(panel.grid.major = element_line(color = 'grey', linetype = 'solid', size = 0.2),
        # panel.grid.minor = element_line(color = 'black', linetype = 'dotted', size = 0.2),
        axis.text = element_text(color = 'black'),
        legend.position = 'bottom', legend.direction = 'horizontal', legend.box = 'vertical', 
        legend.margin = unit(x = 0.005, units = 'in'))
print(shrunkenVolcano)

# source('~/projects/style/R/general.R')
ggsave(filename = getOutFilename(name = 'volcanoPlot', 
                                 saveTo = plotSavePath, proj = 'noct', 
                                 suffix = 'svg'), plot = shrunkenVolcano, 
       width = 7, height = 6, units = 'in')
ggsave(filename = getOutFilename(name = 'volcanoPlot', 
                                 saveTo = plotSavePath, proj = 'noct', 
                                 suffix = 'pdf'), plot = shrunkenVolcano, 
       width = 7, height = 6, units = 'in')

print(plotData[which.max(plotData$negLog10Adjp), ])
print(plotData[which.min(plotData$log2FoldChange), ])
markedRedPointSet <- plotData[plotData$highChange == 2, ]
print(markedRedPointSet)
dim(markedRedPointSet)

library('biomaRt')
# list available marts
head(biomaRt::listMarts(host = 'www.ensembl.org'), 10)
# list all available filters for dataset: hsapiens_gene_ensembl
head(biomaRt::listFilters(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", 
                                              mart    = useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")
                                              )), 10)
filterList <- biomaRt::listFilters(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", 
                                                       mart = useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")))
filterList[startsWith(x = filterList$name, prefix = 'ensembl'), ]
filterOfInterest <- c('ensembl_gene_id')
attributeList <- biomaRt::listAttributes(biomaRt::useDataset(dataset = "hsapiens_gene_ensembl", 
                                                             mart = useMart("ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org")), page = 'feature_page')
attributeList[endsWith(x = attributeList$name, suffix = 'description'), ]
attributeList[endsWith(x = attributeList$name, suffix = 'name'), ]
head(attributeList)
attributeOfInterest <- c('ensembl_gene_id', 'external_gene_name', 'description') 

# annotating 
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
markedRedPointAnnot <- getBM(attributes = attributeOfInterest, filters = filterOfInterest, values = rownames(markedRedPointSet), mart = ensembl)
markedRedPointAnnot
# combining results and annotations
print(dim(markedRedPointSet))
print(dim(markedRedPointAnnot))
statsOfInterest<- c('log2FoldChange', 'padj', 'negLog10Adjp')
markedPointRes <- markedRedPointSet[, statsOfInterest]
markedRedPointAnnotMatched <- markedRedPointAnnot[base::match(x = rownames(markedRedPointSet), table = markedRedPointAnnot$ensembl_gene_id), ]
markedPointRes$geneName <- markedRedPointAnnotMatched$external_gene_name
markedPointRes$description <- markedRedPointAnnotMatched$description
markedPointRes$ensemblGeneId <- rownames(markedPointRes)
# markedPointRes <- markedPointRes[, c('ensemblGeneId', 'log2FoldChange', 'padj', 'negLog10Adjp', 'geneName', 'description')]
markedPointRes <- markedPointRes[order(markedPointRes$log2FoldChange, decreasing = FALSE), c('ensemblGeneId', 'log2FoldChange', 'padj', 'negLog10Adjp', 'geneName', 'description')]
print(head(markedPointRes))
write.table(x = markedPointRes, file = getOutFilename(name = 'volcanoPlotOver4Fold', suffix = 'csv'), quote = TRUE, row.names = FALSE, sep = ',', na = '', append = FALSE)

####### most significantly changed genes #######
# in NOCT overexpressed samples 
resLFCOrdered <- resLFC
dim(resLFCOrdered)
resLFCOrdered <- resLFCOrdered[!is.na(resLFCOrdered$padj), ]
dim(resLFCOrdered)
resLFCOrdered <- resLFCOrdered[resLFCOrdered$padj < 0.1, ]
dim(resLFCOrdered)
resLFCOrdered <- resLFCOrdered[order(resLFCOrdered$log2FoldChange, decreasing = TRUE), ]
print(head(as.data.frame(resLFCOrdered), 10))
print(tail(as.data.frame(resLFCOrdered), 10))
statsOfInterest <- c('log2FoldChange', 'padj')
sigGeneSet <- as.data.frame(resLFCOrdered)[, statsOfInterest]
dim(sigGeneSet)
print(head(sigGeneSet))
sigGeneAnnot <- getBM(attributes = attributeOfInterest, filters = filterOfInterest, values = rownames(sigGeneSet), mart = ensembl)
print(head(sigGeneAnnot))
dim(sigGeneAnnot)
sigGeneRes <- sigGeneSet
sigGeneAnnotMatched <- sigGeneAnnot[base::match(x = rownames(sigGeneSet), table = sigGeneAnnot$ensembl_gene_id), ]
sigGeneRes$geneName <- sigGeneAnnotMatched$external_gene_name
sigGeneRes$description <- sigGeneAnnotMatched$description
sigGeneRes$ensemblGeneId <- rownames(sigGeneRes)
print(head(sigGeneRes))
print(tail(sigGeneRes))
colnames(sigGeneRes)
sigGeneRes <- sigGeneRes[, c('ensemblGeneId', 'log2FoldChange', 'padj', 
                             'geneName', 'description')]
print(head(sigGeneRes))
write.table(x = sigGeneRes, 
            file = getOutFilename(name = 'sigGenes', suffix = 'csv'), 
            quote = TRUE, row.names = FALSE, sep = ',', na = '', append = FALSE)

