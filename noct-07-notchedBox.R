# Project/code: NOCT
# Start date: Tue Dec 24 14:17:14 2019
# Objective: plot notched boxplot for mito- or non-mito-related gene expression
# changes
# --------------
# Author: diaorch
# Modification date:  Tue Dec 24 14:17:14 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------

source('~/projects/styles/R/general.R')
tblFilename <- '06-tpm/20191224-noct-categorizedTable.addUnknown.TPM.DE.csv'

# print(cbPaletteDark)
# scales::show_col(cbPaletteDark)
# print(cbPaletteLight)
# scales::show_col(cbPaletteLight)

library('ggplot2')

tbl <- read.table(file = tblFilename, 
                  header = TRUE, sep = ',', quote = '"', stringsAsFactors = FALSE)
print(head(tbl))
print(table(tbl$category))

plotData <- tbl

table(plotData$category)
# wilcox.test(x = plotData$log2FoldChange[plotData$category == 'Unknown'], mu = 0, alternative = 'two.sided', conf.int = TRUE)
wilcox.test(x = plotData$log2FoldChange[plotData$category == 'Unknown'], 
            y = plotData$log2FoldChange[plotData$category == 'Mitochondrial-encoded gene'], 
            alternative = 'two.sided', conf.int = TRUE)

library('ggsignif')

compare_pairs <- list(
  c('Unknown', 'Nuclear-encoded with low-confidence mitochondrially localized gene'), 
  c('Nuclear-encoded with low-confidence mitochondrially localized gene', 'Nuclear-encoded with high-confidence mitochondrially localized gene'), 
  c('Nuclear-encoded with high-confidence mitochondrially localized gene', 'Mitochondrial-encoded gene')
)

geneCountTbl <- table(plotData$category)
geneCountTbl['Mitochondrial-encoded gene']

category_label <- c(
  paste('Unknown', '\n(', as.character(geneCountTbl['Unknown']), ' genes)', sep = ''),
  paste('Nuclear-encoded\nlow-confidence mito.-localized', '\n(', as.character(geneCountTbl['Nuclear-encoded with low-confidence mitochondrially localized gene']), ' genes)', sep = ''),
  paste('Nuclear-encoded\nhigh-confidence mito.-localized', '\n(', as.character(geneCountTbl['Nuclear-encoded with high-confidence mitochondrially localized gene']), ' genes)', sep = ''),
  paste('Mitochondrially-encoded', '\n(', as.character(geneCountTbl['Mitochondrial-encoded gene']), ' genes)', sep = '')
  )
names(category_label) <- c('Unknown', 'Nuclear-encoded with low-confidence mitochondrially localized gene',
                           'Nuclear-encoded with high-confidence mitochondrially localized gene', 'Mitochondrial-encoded gene')

ggplot(data = plotData, aes(x = category, y = log2FoldChange)) + 
  geom_jitter(size = 0.4, width = 0.25, color = 'black') + 
  geom_boxplot(fill = 'transparent', color = 'red', 
               outlier.alpha = 0, outlier.size = 0, outlier.color = cbPaletteDark[8], # varwidth = TRUE
               notch = TRUE, notchwidth = 0.5,
               width = 0.4) + 
  # geom_text(aes(x = category, label = category), y = 0, hjust = 0.5, vjust = -2, 
  #           color = cbPaletteDark[7], size = 4.5) +
  geom_signif(
    comparisons = compare_pairs, 
    map_signif_level = TRUE, step_increase = 0.075, vjust = 5, hjust = -1, 
    test = 'wilcox.test', test.args = list(alternative = 'two.sided', conf.int = TRUE)
  ) +
  scale_x_discrete(label = category_label) + 
  scale_y_continuous(breaks = seq(-3, 3, 1)) +
  coord_flip() +
  theme_bw() +
  theme(
    axis.text = element_text(color = 'black'),
    # axis.text.x = element_text(size = 4), 
    # axis.text.y = element_blank()
    axis.title.y = element_blank()
  ) + 
  labs(
    # title = 'Log2 fold change distributions for mitochondrial/nuclear encoded, \nhigh/low confidence mitochondrially localized genes', 
    y = 'Log2 fold change', x = 'Gene category'
  )

ggsave(filename = getOutFilename(name = 'mitoL2FC.notchedBox', saveTo = '~/projects/lab-notebook/figs/', 
                                 proj = 'noct', suffix = '.svg'), 
       plot = last_plot(), width = 7, height = 4)
