# Project/code: NOCT
# Start date: Tue Dec 24 14:17:14 2019
# Objective: plot notched boxplot for mito- or non-mito-related gene expression
# changes
# --------------
# Author: diaorch
# Modification date:  Tue Dec 24 14:17:14 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
#   1. revision (updated from noct-06-notchedBox.R under main script directory)
# --------------

library('ggplot2')
library('ggsignif')
source(paste(Sys.getenv('HOME'), '/projects/styles/R/general.R', sep = ''))

tblFilename <- '06-tpm/20191224-noct-categorizedTable.addUnknown.TPM.DE.csv'
tbl <- read.table(file = tblFilename, 
                  header = TRUE, sep = ',', quote = '"', stringsAsFactors = FALSE)
# print(head(tbl))
# print(table(tbl[, 'category']))
# print(table(tbl[!is.na(tbl[, 'padj']), 'category']))

plotData <- tbl
wilcox.test(x = plotData$log2FoldChange[plotData$category == 'Unknown'], 
            y = plotData$log2FoldChange[plotData$category == 'Mitochondrial-encoded gene'], 
            alternative = 'two.sided', conf.int = TRUE)

compare_pairs <- list(
  c('Unknown', 'Nuclear-encoded with low-confidence mitochondrially localized gene'), 
  c('Nuclear-encoded with low-confidence mitochondrially localized gene', 'Nuclear-encoded with high-confidence mitochondrially localized gene'), 
  c('Nuclear-encoded with high-confidence mitochondrially localized gene', 'Mitochondrial-encoded gene')
)

geneCountTbl <- table(plotData[, 'category'])
geneCountTblSubset <- table(plotData[!is.na(plotData[, 'padj']), 'category'])
geneCountTbl['Mitochondrial-encoded gene']
geneCountTblSubset['Mitochondrial-encoded gene']

formatCategoryLabel <- function(categoryName){
  paste(gsub(x = gsub(x = categoryName, pattern = ' gene$', replacement = '', perl = TRUE), 
             pattern = 'confidence ', replacement = 'confidence\n', fixed = TRUE), 
        '\n(', as.character(geneCountTblSubset[categoryName]), '/', as.character(geneCountTbl[categoryName]), ' genes)', 
        sep = '')
}
formatCategoryLabel('Mitochondrial-encoded gene')

categoryList <- unique(plotData[, 'category'])
categoryList <- sort(categoryList)
print(categoryList)
categoryLabel <- sapply(categoryList, formatCategoryLabel)
# names(categoryLabel) <- categoryList
print(categoryLabel)

pBox <- ggplot(data = plotData, aes(x = category, y = log2FoldChange)) + 
  geom_jitter(size = 0.4, width = 0.25, color = 'black') + 
  geom_boxplot(fill = 'transparent', color = 'red', 
               outlier.alpha = 0, outlier.size = 0, outlier.color = cbPaletteDark[8],
               notch = TRUE, notchwidth = 0.5,
               width = 0.4) + 
  geom_signif(
    comparisons = compare_pairs, 
    map_signif_level = TRUE, step_increase = 0.075,  vjust = 5, hjust = -0.6, 
    test = 'wilcox.test', test.args = list(alternative = 'two.sided', conf.int = TRUE)
  ) +
  scale_x_discrete(label = categoryLabel) + 
  scale_y_continuous(breaks = seq(-3, 3, 1)) +
  coord_flip(clip = 'off') +
  theme_bw() +
  theme(
    text = element_text(family = 'Arial'), 
    axis.text = element_text(color = 'black'),
    axis.title.y = element_blank()
  ) + 
  labs(
    y = 'Log2 fold change', x = 'Gene category'
  )

savePlotFormat <- function(fmt){
  fmtDeviceList <- list('png' = png, 'svg' = svg, 'pdf' = cairo_pdf, 'eps' = cairo_ps)
  ggsave(filename = getOutFilename(name = 'mitoL2FC.notchedBox', saveTo = '~/projects/lab-notebook/figs/', 
                                     proj = 'noct', suffix = fmt), 
         plot = pBox, 
         width = 7, height = 4, dpi = 300, device = fmtDeviceList[[fmt]])
 
}
sapply(X = c('pdf', 'svg', 'eps'), FUN = savePlotFormat)
