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

# save.image('20191001-volcano_plotting.RData')
# 
# ggsave(plot = shrunkenVolcano, 
#        filename = '20191001-noct-overexpVolcano.svg',  
#        device = 'svg'#, width = unit(x = 5.5, units = 'in'), height = unit(4.25, 'in')
#        )
# ggsave(plot = shrunkenVolcano, 
#        filename = '20191001-noct-overexpVolcano.pdf', 
#        device = 'pdf'#, width = unit(x = 5.5, units = 'in'), height = unit(4.25, 'in')
# )



#### plot GO term specific gene expression level ####
goDetected <- read.table(
  file = '~/data/noct/humanOE/07-ipageData/ensemblDB/20181004-noct-deseq2L2FC.exp_bin9_PAGE/pvmatrix.go.tsv', 
  sep = '\t', header = FALSE, stringsAsFactors = FALSE)
print(head(goDetected))
getGOaccn <- function(x){
  return(strsplit(x = x, split = ' ', fixed = TRUE)[[1]][1])
}
goDetected$GO <- sapply(X = goDetected$V1, FUN = getGOaccn)
# goDetected$categoryCount <- c(1:5)
print(head(goDetected))

ipageDB <- file('~/data/noct/humanOE/06-ipageDb/human_ensg/human_ensg_index.txt')
# ipageData <- readLines(ipageDB)
fieldCount <- max(count.fields(ipageDB, sep = '\t'))
ipageEnsgIndex <- 
  read.table(ipageDB, header = FALSE, sep = '\t', 
           col.names = paste0('V', as.character(1:fieldCount)),
           fill = TRUE, stringsAsFactors = FALSE)
print(dim(ipageEnsgIndex))
print(head(ipageEnsgIndex))
close(ipageDB)

goTermTbl <- ipageEnsgIndex[, 2:ncol(ipageEnsgIndex)]
# print(head(goTermTbl))
# print(dim(goTermTbl))
# allGo <- unique(unlist(x = goTermTbl))
# print(length(allGo))
# print(head(allGo))

# only check the GO terms iPAGE detected, input from read.table
goEnsgMapping <- vector(mode = 'character', length = nrow(goDetected))
names(goEnsgMapping) <- goDetected$GO
print(head(goEnsgMapping))

allGoEnsg <- data.frame(
  matrix(data = NA, nrow = 0, ncol = 8, 
         dimnames = list(
           NULL, 
           c('GO', 'ENSG', 'log2FoldChange', 'padj', 'negLog10Q', 'hgnc_symbol', 'description', 'sig'))))

for (g in names(goEnsgMapping)){
  # g <- names(goEnsgMapping)[1]
  # g <- 'GO:0000776'
  # print(g)
  gFlag <- (!is.na(goTermTbl) & goTermTbl == g)
  # print(head(goTermTbl))
  # print(dim(gFlag))
  # print(gFlag[1:5, 1:5])
  gFlagEnsg <- apply(X = gFlag, MARGIN = 1, FUN = any)
  # print(length(gFlagEnsg))
  # print(sum(gFlagEnsg))
  # print(head(ipageEnsgIndex))
  gEnsgIndex <- ipageEnsgIndex[gFlagEnsg, 1]
  # print(g)
  print(length(gEnsgIndex))
  gEnsg <- paste(as.character(gEnsgIndex), collapse = ',')
  # stringr::str_count(string = gEnsg, pattern = 'ENSG')
  # print(gEnsg)
  goEnsgMapping[g] <- gEnsg
  ensgExpr <- plotData[plotData$ENSG %in% gEnsgIndex, ]
  print(nrow(ensgExpr))
  print('===')
  ensgExpr$sig <- ifelse(ensgExpr$padj < 0.05, 'sig', 'non-sig')
  ensgExpr$GO <- g
  # print(head(ensgExpr))
  allGoEnsg <- rbind(allGoEnsg, ensgExpr)
}
# print(head(goEnsgMapping))
print(head(allGoEnsg))

# save(goDetected, allGoEnsg, ipageEnsgIndex, goEnsgMapping, 
#      file = '~/data/noct/humanOE/07-ipageData/20191211-noct-goTermTable.RData')

# load('~/data/noct/humanOE/07-ipageData/20191211-noct-goTermTable.RData')
# print(head(goDetected))
# goDetected[grep(x = goDetected$V1, pattern = 'neurotransmitter secretion'), ]
# goDetected[grep(x = goDetected$V1, pattern = 'cell adhesion'), ]
# goDetected[grep(x = goDetected$V1, pattern = 'tyrosine kinase activity'), ]
# goDetected[grep(x = goDetected$V1, pattern = 'microglial cell activation'), ]

# check interesting GO terms
goCheck <- c('GO:0007204', 'GO:0007187', 'GO:0042446', 'GO:0042692',
             'GO:0070469', 'GO:0007156', 'GO:0005249', 'GO:0045668',
             'GO:0004714', 'GO:0004970', 'GO:0007269' 
             )
goCheckManuscript <- c('GO:0007156', 'GO:0004714', 'GO:0045668', 
                       'GO:0004970', 'GO:0007269', 'GO:0005249',
                       'GO:0001774')
goCheck <- goCheckManuscript

checkGoEnsg <- allGoEnsg[allGoEnsg$GO %in% goCheck, ]
print(head(checkGoEnsg))
print(dim(checkGoEnsg))
print(summary(checkGoEnsg$log2FoldChange))
# print(head(gnMapping))

# checkGoEnsg <- merge(x = checkGoEnsg, y = gnMapping, by.x = 'ENSG', by.y = 'ensembl_gene_id', all.x = TRUE, all.y = FALSE)
print(head(goDetected))

goEnsgPlotData <- checkGoEnsg[!is.na(checkGoEnsg$padj), ]
print(head(goEnsgPlotData))
goEnsgPlotData[goEnsgPlotData$hgnc_symbol == 'Endogenous NOCT', ]
goEnsgPlotData[goEnsgPlotData$hgnc_symbol == 'Endogenous NOCT', 'hgnc_symbol'] <- 'Endo. NOCT'


ggplot(data = goEnsgPlotData) +
  geom_bar(aes(x = reorder(hgnc_symbol, - log2FoldChange), y = log2FoldChange, fill = sig), stat = 'identity') + 
  facet_wrap(~ GO, scales = 'free_x') + 
  # facet_wrap( ~ regulation + representation + categoryCount, scales = 'free_x') + 
  labs(
    title = 'All genes'
  ) + 
  theme(
    axis.text.x = element_text(angle = 90)
  )

# https://stackoverflow.com/a/3935429/10683747
wrapText <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}

goManual <- c('GO:0007204', 'GO:0070469', 'GO:0007156', 'GO:0007269')
goManual <- c()
# pdf('20191211-noct-goTermsOfInterst.mentioned.pdf', height = 5)
pLabel <- c('A', 'B', 'C', 'D', 'E', 'F', 'G')
names(pLabel) <- goCheck
library('grid')
library('ggplot2')

pList <- list()
for(goTerm in goCheck){
  # goTerm <- goCheck[1]
  if (!(goTerm %in% goManual)){
    print(goTerm)
    if (goTerm == 'GO:0007156'){
      t <- wrapText(x = goDetected$V1[goDetected$GO == goTerm], width = 100)
    } else if(goTerm %in% c('GO:0004970', 'GO:0001774')) {
      t <- wrapText(x = goDetected$V1[goDetected$GO == goTerm], width = 30)
    } else {
      t <- wrapText(x = goDetected$V1[goDetected$GO == goTerm], width = 50)
    }
    
    indivPlotData <- goEnsgPlotData[goEnsgPlotData$GO == goTerm, ]
    p <- ggplot(data = indivPlotData) +
      geom_hline(yintercept = 0, color = 'black') +  
      geom_bar(aes(x = reorder(hgnc_symbol, - log2FoldChange), y = log2FoldChange, fill = sig), 
               stat = 'identity', 
               color = 'black', size = 0.25,
               show.legend = FALSE) + 
      scale_y_continuous(limits = c(-2.5, 2.5)) + 
      scale_fill_manual(
        name = 'Significance:',
        labels = c('sig' = 'Significant', 'non-sig' = 'Non-significant'), 
        values = c('sig' = 'black', 'non-sig' = 'white')) + 
      labs(
        title = t, 
        x = 'Gene name', y = 'log2(fold change)', tag = pLabel[[goTerm]]
      ) + 
      theme_bw() + 
      theme(
        title = element_text(size = 6), 
        axis.text.x = element_text(
          angle = 90, 
          # size = min(max(125 / nrow(indivPlotData), 3), 7), 
          size = 4,
          hjust = 1, vjust = 0.5), 
        axis.text.y = element_text(size = 5), 
        plot.tag = element_text(size = 10)
      )
    pList[[goTerm]] <- p
  }
}

library('gridExtra')
pdf('20191226-noct-goTermsOfInterst.mentioned.updated.pdf', height = 5)
# svg('20191212-noct-goTermsOfInterst.mentioned.updated.svg', height = 5)
grid.arrange(grobs = pList, # ncol = 2, 
             layout_matrix = cbind(c(1, 2, 5), 
                                   c(1, 2, 5), 
                                   c(1, 2, 5), 
                                   c(1, 3, 6), 
                                   c(1, 3, 6), 
                                   c(1, 3, 6), 
                                   c(1, 4, 7), 
                                   c(1, 4, 7)))
dev.off()

# install.packages('cowplot')
library('cowplot')
row2 <- plot_grid(pList[[2]], pList[[3]], pList[[4]], rel_widths = c(3, 3, 2), align = 'h', nrow = 1)
row3 <- plot_grid(pList[[5]], pList[[6]], pList[[7]], rel_widths = c(3, 3, 2), align = 'h', nrow = 1)
pdf('20191226-noct-goTermsOfInterst.pdf', height = 5)
plot_grid(pList[[1]], row2, row3, ncol = 1)
dev.off()
svg('20191226-noct-goTermsOfInterst.svg', height = 5)
plot_grid(pList[[1]], row2, row3, ncol = 1)
dev.off()

