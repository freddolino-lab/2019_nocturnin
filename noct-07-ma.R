# Project/code: Noct, overexpression RNA-seq data
# Start date: Sat Nov 23 23:04:02 2019
# Objective: plot MA-plot per collaborator request
# --------------
# Author: diaorch
# Modification date:  Sat Nov 23 23:04:02 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------

tbl <- read.table(file = '20191108-TPM.alignedBySTAR.rounded.DESeqStats.csv', 
                  sep = ',', header = TRUE, stringsAsFactors = FALSE)
print(head(tbl))

# significance threshold
alpha <- 0.05

#### basic stats ####
# ratio, as M in "MA plot"
m <- tbl$log2FoldChange
# average log2-pseudo-counts
# red: test samples average TPM
r <- log2(tbl$Del_3F_11_S5_TPM + 1) + log2(tbl$Del_3F_12_S6_TPM + 1) + log2(tbl$Del_3F_4_S4_TPM + 1)
# green: control samples TPM
g <- log2(tbl$GST_3F_2_S1_TPM + 1) + log2(tbl$GST_3F_8_S2_TPM + 1) + log2(tbl$GST_3F_9_S3_TPM + 1)
# average, as A in "MA plot"
a <- (r + g) / 2
print(head(a))

#### plotting ####
library('ggplot2')
print(head(tbl))
plotData <- data.frame(
  gene_id = tbl$gene_id, 
  log2FoldChange = tbl$log2FoldChange, 
  m = m, 
  a = a, 
  sig = (tbl$padj) < alpha, 
  hgnc_symbol = tbl$hgnc_symbol, 
  stringsAsFactors = FALSE
)

pcrGenes <- c('DLGAP1', 'RORB', 'TMEFF2', 'PTPRZ1', 'FRRS1L')

print(head(plotData))
subset(plotData, hgnc_symbol == 'NOCT')
subset(plotData, log2FoldChange < -2 & !(hgnc_symbol %in% pcrGenes))
# overexpressed NOCT
plotData[grep(pattern = 'NOCT', x = plotData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- '"NOCT "~Delta~"(2-15)-3F"'
plotData[grep(pattern = 'ENSG00000151014', x = plotData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- 'Endogenous NOCT'
plotData[grep(pattern = 'GST', x = plotData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- 'GST-3F'
noctGenes <- c('"NOCT "~Delta~"(2-15)-3F"', 'Endogenous NOCT', 'GST-3F')

library('ggrepel')
plotData <- plotData[complete.cases(plotData), ]
print(head(plotData))
labelData <- subset(plotData, abs(log2FoldChange) > 2 | hgnc_symbol %in% noctGenes)
print(labelData)
labelData[grep(pattern = 'NOCT', x = labelData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- '"NOCT "~Delta~"(2-15)-3F"'
labelData[grep(pattern = 'ENSG00000151014', x = labelData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- '"Endogenous "~NOCT'
labelData[grep(pattern = 'GST', x = labelData$gene_id, fixed = TRUE), 'hgnc_symbol'] <- '"GST-3F"'
noctGenes <- c('"NOCT "~Delta~"(2-15)-3F"', '"Endogenous "~NOCT', '"GST-3F"')
labelData$color <- 'black'
labelData$color[labelData$hgnc_symbol %in% noctGenes] <- 'blue'
print(labelData)
labelData$box <- 'box'
labelData$box[labelData$hgnc_symbol %in% c('LIPG', 'KRT19', '"GST-3F"', 'ADGRB3', 'FBXL7')] <- 'noBox'

print(labelData$hgnc_symbol)
labelPos <- data.frame(
  hgnc_symbol = labelData$hgnc_symbol, 
  # mark         x              x                         x     x     x             x      x
  # gene         1      2      3      4      5       6    7     8     9     10     11     12
  labelX  = c( 3.7,  -1.7,  12.5,    16,    25,   -1.7,  20,  3.5,  7.75,    -1,   11, 12.5), 
  labelY  = c( 2.5,  -2.5, -2.25,  -2.7, -1.25, -2.951, 2.8, -3.3,  -2.9, -1.75, -1.9,    2), 
  segEndX = c(3.55,     2,   8.5,  12.1,    20,    2.2,  17,  3.5,  7.75,  -2.2,  7.5,  7.5), 
  segEndY = c(2.35, -2.35, -2.25, -2.65, -1.01,  -2.95, 2.8, -3.2, -2.71, -1.95, -1.9, 1.75)
)
labelData$labelX <- labelPos$labelX
labelData$labelY <- labelPos$labelY
labelData$segEndX <- labelPos$segEndX
labelData$segEndY <- labelPos$segEndY
print(labelData)

ggplot(data = plotData, aes(x = a, y = m)) + 
  geom_hex(binwidth = c(1, 0.2)) + 
  geom_text(data = labelData[labelData$box == 'noBox', ], 
            aes(x = labelX, y = labelY, label = hgnc_symbol, color = color), 
            parse = TRUE, size = 6.5) + 
  geom_label(data = labelData[labelData$box == 'box', ], 
             aes(x = labelX, y = labelY, label = hgnc_symbol, color = color), 
             parse = TRUE, size = 6.5, fill = 'transparent') + 
  geom_segment(data = labelData, 
               aes(x = a, xend = segEndX, y = m, yend = segEndY, color = color), 
               size = 0.7) + 
  geom_point(data = plotData[plotData$sig == TRUE, ], 
             aes(color = sig, alpha = sig), size = 0.75, show.legend = FALSE) + 
  scale_x_continuous(breaks = seq(0, 60, 5), limits = c(-3, NA)) + 
  scale_color_manual(name = 'sig', breaks = c(),
                     values = c('TRUE' = 'red', 'black' = 'black', 'blue' = 'blue'),
                     labels = c('TRUE' = 'Significant'),
                     guide = FALSE) +
  scale_alpha_manual(name = 'sig', breaks = c('TRUE'),
                     values = c('TRUE' = 1, 'FALSE' = 0),
                     labels = c('TRUE' = 'Significant'),
                     guide = FALSE) +
  scale_fill_viridis_c(
    trans = 'log10'
  ) + 
  labs(
    y = expression(Log[2]~(NOCT[OE]/GST[OE])), 
    x = expression(Mean(log[2](TPM)))
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18), 
    panel.grid.major = element_line(color = 'grey', linetype = 'solid', size = 0.2),
    # panel.grid.major = element_line(color = 'black', linetype = 'solid', size = 0.2),
    # panel.grid.minor = element_line(color = 'black', linetype = 'solid', size = 0.2),
    plot.margin = margin(t = 10, b = 10, l = 10, r = 10), 
    axis.text = element_text(size = 18, color = 'black'),
    legend.position = c(1, 0), # legend.direction = 'vertical',
    legend.justification = c(1.1, -0.1), 
    legend.text = element_text(size = 18, margin = margin(t = 0, l = 0, r = 0, b = 0), hjust = 0.5),
    # legend.text.align = 0, 
    legend.background = element_rect(fill = 'white', color = 'black'),
    legend.margin = margin(t = 5, l = 5, r = 5, b = 5)
  ) +
  guides(
    fill = guide_colorbar(
      title = 'Counts in bins: ', 
      direction = 'horizontal', label.position = 'bottom', title.position = 'top', 
      barwidth = unit(x = 150, units = 'pt') 
    )
  )

ggsave(filename = '20191126-noct-ma.svg', 
       width = 7 * 1.5, height = 4 * 1.5) 
ggsave(filename = '20191126-noct-ma.pdf', 
       width = 7 * 1.5, height = 4 * 1.5) 
