# Project/code: Nocturnin
# Start date: Thu Jul 18 18:11:20 2019
# Objective: plot gene-level differential expression results
# --------------
# Author: diaorch
# Modification date:  Thu Jul 18 18:11:20 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------

library('ggplot2')
library('reshape2')
library('ggrepel')
library('biomaRt')

#### read DESeq2 results ####
deResTbl <- read.table(file = '20190315-noct-deseq2_lfcShrink_ordered_table.csv', 
                       header = TRUE, sep = ',', quote = '"', stringsAsFactors = FALSE)
print(dim(deResTbl))
print(head(deResTbl))

#### plot simple volcano plot ####
plotData <- deResTbl[, c('log2FoldChange', 'padj')]
plotData$ENSG <- deResTbl$X
plotData$negLog10Q <- - log10(plotData$padj)
v <- ggplot(data = plotData) + 
  geom_point(aes(x = log2FoldChange, y = negLog10Q))
print(v)

#### add threshold coloring to plot data ####
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
print(head(plotData))
print(summary(plotData$log2FoldChange))
print(summary(plotData$negLog10Q))

alpha <- 0.05

#### add gene names to plot data ####
# conda install -c bioconda bioconductor-biomart
ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
listAttributes(mart = ensembl, page = 'feature_page', what = 'description')
head(listAttributes(mart = ensembl, page = 'feature_page'), 20)
gnMapping <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
                   filters = 'ensembl_gene_id' , values = plotData$ENSG, 
                   mart = ensembl)
print(head(gnMapping))
plotData <- merge(x = plotData, y = gnMapping, 
                  by.x = 'ENSG', by.y = 'ensembl_gene_id', 
                  all.x = TRUE, all.y = FALSE)
print(head(plotData))
print(plotData[1300:1305, ])
# plotData$hgnc_symbol[plotData$highChange != 2] <- ''

pcrGenes <- c('DLGAP1', 'RORB', 'TMEFF2', 'PTPRZ1', 'FRRS1L')

print(head(plotData))
subset(plotData, hgnc_symbol == 'NOCT')
subset(plotData, log2FoldChange < -2 & !(hgnc_symbol %in% pcrGenes))
# overexpressed NOCT
plotData[grep(pattern = 'NOCT', x = plotData$ENSG, fixed = TRUE), 'hgnc_symbol'] <- '"NOCT "~Delta~"(2-15)-3F"'
plotData[grep(pattern = 'ENSG00000151014', x = plotData$ENSG, fixed = TRUE), 'hgnc_symbol'] <- 'Endogenous NOCT'
plotData[grep(pattern = 'GST', x = plotData$ENSG, fixed = TRUE), 'hgnc_symbol'] <- 'GST-3F'
noctGenes <- c('"NOCT "~Delta~"(2-15)-3F"', 'Endogenous NOCT', 'GST-3F')

#### plot volcano plot ####

shrunkenVolcano <- 
  ggplot(data = plotData, aes(x = log2FoldChange, y = negLog10Q, color = highChange)) + 
  geom_point(size = 1) +
  geom_hline(aes(yintercept = -log10(alpha), linetype = '-log10(alpha)'), color = 'black', 
             show.legend = FALSE) + 
  geom_text(data = data.frame(y = -log10(alpha), 
                              label = paste('-log[10](', as.character(alpha), ')', sep = '')), 
            parse = TRUE, 
            aes(label = label, y = y), x = 2.4, color = 'black', 
            hjust = 0, vjust = -0.1) +
  geom_text_repel(data = subset(plotData, log2FoldChange < -2 & !(hgnc_symbol %in% pcrGenes)), 
                  aes(label = hgnc_symbol), size = 3.5, 
                  segment.color = 'darkgrey', point.padding = 0, 
                  xlim = c(NA, -2.5), # ylim = c(17.5, 52.5),
                  # direction = 'y', 
                  show.legend = FALSE) +
  geom_label_repel(data = subset(plotData, log2FoldChange < -2 & (hgnc_symbol %in% pcrGenes)),
                   aes(label = hgnc_symbol), size = 3.5,
                   segment.color = 'darkgrey', point.padding = 0.1,
                   fill = 'transparent', 
                   # xlim = c(-3, -2.5), ylim = c(0, NA),
                   xlim = c(NA, NA), ylim = c(27, NA),
                   direction = 'x',
                   show.legend = FALSE) +
  geom_text_repel(data = subset(plotData, log2FoldChange > 2), 
                  aes(label = hgnc_symbol), size = 3.5, 
                  segment.color = 'darkgrey', point.padding = 0, 
                  ylim = c(15, 52.5),
                  direction = 'y',
                  show.legend = FALSE) + 
  geom_text_repel(data = subset(plotData, hgnc_symbol == 'GST-3F'),
                   aes(label = hgnc_symbol), size = 3.5, color = 'blue',
                   segment.color = 'darkgrey', point.padding = 0,
                   xlim = c(NA, -2), ylim = c(7, 7),
                   # direction = 'y',
                   show.legend = FALSE) +
  geom_label_repel(data = subset(plotData, hgnc_symbol == 'Endogenous NOCT'),
                   aes(label = hgnc_symbol), size = 3.5, color = 'blue',
                   segment.color = 'darkgrey', point.padding = 0,
                   fill = 'transparent',
                   xlim = c(-0.9, NA), ylim = c(7, NA),
                   # direction = 'y',
                   show.legend = FALSE) +
  geom_label_repel(data = subset(plotData, hgnc_symbol == '"NOCT "~Delta~"(2-15)-3F"'),
                   parse = TRUE, 
                   aes(label = hgnc_symbol), size = 3.5, color = 'blue',
                   segment.color = 'darkgrey', point.padding = 0,
                   fill = 'transparent',
                   xlim = c(1.6, 1.6), ylim = c(8, 8),
                   # direction = 'y',
                   show.legend = FALSE) +
  scale_x_continuous(# name = 'Shrunken log[2](fold change)', 
                     breaks = seq(-3, 3, 1)) + 
  scale_y_continuous(# name = '-log10(adjusted p-value)', 
                     breaks = seq(0, 55, 5)) + 
  scale_color_manual(# name = 'Expression change level: ', 
                     values = c('grey50', 'black', 'red'),
                     breaks = c(0, 1, 2),
                     labels = c('<= 2 fold', '2 < fold change <= 4 fold', '> 4 fold')) +
  scale_linetype_manual(# name = 'Significance', 
                        values = c('-log10(alpha)' = 'dashed')) +
  labs(# title = 'Differential Gene Expression with NOCT ∆2-15 Overexpression (NOCT/GST)',
       # subtitle = 'Data analyzed by STAR, HTSeq, and DESeq2.',
       x = expression(Shrunken~log[2](fold~change)),
       y = expression(-log[10](adjusted~`p-value`))
       # caption = '6359 NA-values in adjusted p-values by DESeq2 are removed',
       # color = 'Expression change levels: ',
       # linetype = 'Significance:'
       ) +
  theme_bw() + 
  theme(panel.grid.major = element_line(color = 'grey', linetype = 'solid', size = 0.2),
        # panel.grid.minor = element_line(color = 'black', linetype = 'dotted', size = 0.2),
        axis.text = element_text(color = 'black'),
        legend.position = 'bottom', legend.direction = 'horizontal', 
        # legend.key.width = unit(x = 5, units = 'pt'), 
        # legend.spacing.x = unit(x = 5, units = 'pt')
        legend.text = element_text(margin = margin(t = 0, l = 13, r = 13))
        # legend.margin = margin(t = 50, r = 50, b = 50, l = 50, unit = 'pt')
  ) +
  guides(color = guide_legend(
    title = 'Expression change level:',
    position = 'bottom', direction = 'horizontal', 
    keywidth = unit(x = 50, units = 'pt'), label.position = 'bottom', title.position = 'left'
  ))
print(shrunkenVolcano)

