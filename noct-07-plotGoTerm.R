# Project/code: Nocturnin
# Start date: Thu Jul 18 18:11:20 2019
# Objective: plot gene-level differential expression results
# --------------
# Author: diaorch
# Modification date:  20200611
# --------------

source('~/projects/styles/R/general.R')
library('ggplot2')
library('grid')
library('gridExtra')
load('~/data/noct/humanOE/revision/20200614-noct-goToGene.RData')
# allGoEnsg: info of all genes associated with all iPAGE-identified GO terms (96 terms, 4744 genes)
# checkGoEnsg: info for genes associated with GO terms plotted for supplementary figures
# # manually list GO terms of interest
goBone  <- c('GO:0045668', 'GO:0005518', 'GO:0002062')
goMito  <- c('GO:0005747', 'GO:0045333', 'GO:0007005', 'GO:0006631', 'GO:0070125')
goNAD   <- c('GO:0008137', 'GO:0006120', 'GO:0051287')
goNeuro <- c('GO:0007218', 'GO:0005267', 'GO:0005244', 'GO:0051402', 'GO:0001580', 
             'GO:0050909', 'GO:0071805', 'GO:0001843', 'GO:0051928', 'GO:0005262', 
             'GO:0042053')
goBoneFat <- c('GO:0070328')

goCheckRevision <- c(goBone, goMito, goNAD, goNeuro, goBoneFat)
goCheckCategory <- c(rep('Bone', length(goBone)), rep('Mito', length(goMito)), rep('NAD', length(goNAD)), rep('Neuro', length(goNeuro)), rep('BoneFat', length(goBoneFat)))
names(goCheckCategory) <- goCheckRevision
# print(goCheckCategory)
goCheck <- goCheckRevision
# print(length(goCheck))

goDescription <- names(goDetectedAccn)
names(goDescription) <- goDetectedAccn
print(goDescription)

checkGoEnsg <- allGoEnsg[allGoEnsg[, 'GO'] %in% goCheck, ]

# # process plot data
print(ls())
goEnsgPlotData <- checkGoEnsg[!is.na(checkGoEnsg$padj), ]
print(head(goEnsgPlotData, n = 3))
# # 1227 genes associated to the 23 GO terms of interest, 1055 of which had non-NA padj

# # plot distributions of gene expression changes (as log2-fold-changes) for GO terms
# # of interest
# # add new line to semi-automatically wrap text too long for plot titles
wrapText <- function(x, ...) 
{
  paste(strwrap(x, ...), collapse = "\n")
}
# # extract common legend
getPlotLegend <- function(p){
  tmp <- ggplot_gtable(ggplot_build(p))
  legIdx <- which(sapply(tmp$grobs, function(x) x$name == 'guide-box'))
  leg <- tmp$grobs[[legIdx]]
  return(leg)
}

# # one gene did not have HGNC symbol hit, manual added Ensembl gene name
print(which(goEnsgPlotData[, 'hgnc_symbol'] == ''))
print(goEnsgPlotData[which(goEnsgPlotData[, 'hgnc_symbol'] == ''), ])
if (goEnsgPlotData[325, 'hgnc_symbol'] == ''){
  goEnsgPlotData[325, 'hgnc_symbol'] <- 'AC006254.1'
}
print(goEnsgPlotData[325, ])

# # arrange plots by gene counts
pageLayoutMtxList <- list(
  'BoneBoneFat' = list(
  # # Bone/BoneFat
    'GO' = c(goBone, goBoneFat), 
    'pageLayout' = matrix(data = c(1, 1, 1, 3, 3, 3, 2, 2, 2, 2, 4, 4), nrow = 6, ncol = 2, byrow = FALSE), 
    'panelTag' = setNames(c('A', 'C', 'B', 'D'), names(goCheckCategory)[goCheckCategory %in% c('Bone', 'BoneFat')])
), 
  'Mito' = list(
  # # Mito
  'GO' = goMito, 
  'pageLayout' = matrix(data = c(rep(4, 6), rep(5, 6), rep(3, 6), rep(1, 4), rep(2, 2)), nrow = 6, ncol = 4, byrow = FALSE), 
  'panelTag' = setNames(c('D', 'E', 'C', 'A', 'B'), names(goCheckCategory)[goCheckCategory == 'Mito'])
), 
  'NAD' = list(
  # # NAD
  'GO' = goNAD, 
  'pageLayout' = matrix(data = c(rep(1, 6), rep(2, 6), rep(3, 6)), nrow = 6, ncol = 3, byrow = FALSE),
  'panelTag' = setNames(LETTERS[1:3], names(goCheckCategory)[goCheckCategory == 'NAD'])
), 
  'Neuro' = list(
  # # Neuro
  'GO' = goNeuro,
  'pageLayout' = matrix(data = c(rep(3, 6), rep(8, 6), rep(7, 4), rep(5, 2), rep(10, 4), rep(9, 2), rep(1, 3), rep(4, 3), rep(6, 2), rep(2, 3), 11), nrow = 6, ncol = 6, byrow = FALSE), 
  'panelTag' = setNames(LETTERS[1:length(goNeuro)], names(goCheckCategory)[goCheckCategory == 'Neuro'][c(3, 8, 7, 5 ,10, 9, 1, 4, 6, 2, 11)])
  )
)

pListLegend <- NULL
for(i in 1:length(pageLayoutMtxList)){
  goCategory <- names(pageLayoutMtxList)[i]
  print(goCategory)
  pList <- list()
  for (goTerm in pageLayoutMtxList[[goCategory]]$GO){
    goLongText <- goDescription[goTerm]

    indivPlotData <- goEnsgPlotData[goEnsgPlotData$GO == goTerm, ]
    indivPlotData[, 'sig'] <- factor(indivPlotData[, 'sig'], levels = c('sig', 'non-sig'))
    titleText <- wrapText(x = goLongText, width = 50)
    print(titleText)
    p <- ggplot(data = indivPlotData) +
      geom_bar(aes(x = reorder(hgnc_symbol, - log2FoldChange), y = log2FoldChange, fill = sig), 
               stat = 'identity', 
               color = 'black', size = 0.25,
               show.legend = FALSE) + 
      scale_y_continuous(limits = c(-2.5, 1.5), breaks = seq(-2, 1, 1)) + # check ranges, for l2fc (-2.5, 1.5)
      scale_fill_manual(
        name = 'Significance (by DESeq2):',
        labels = c('sig' = 'Significant', 'non-sig' = 'Non-significant'), 
        values = c('sig' = 'black', 'non-sig' = 'white')) + 
      labs(
        title = titleText, 
        x = 'Gene names', 
        y = 'Shrunken log2(fold changes)', 
        tag = pageLayoutMtxList[[goCategory]]$panelTag[goTerm]
      ) + 
      theme_bw() + 
      theme(
        title = element_text(size = 7.5), 
        axis.title = element_text(size = 10), 
        axis.text.y = element_text(
          color = 'black', 
          size = ifelse(nrow(indivPlotData) >= 80, ifelse(nrow(indivPlotData ) >= 120, 3.5, 5), 4.5) + 2,
          hjust = 1, vjust = 0.5), 
        axis.text.x = element_text(size = 10, color = 'black'), 
        plot.tag = element_text(size = 13),
        plot.margin = margin(t = 0.05, l = 0.2, b = 0.05, r = 0.2, 'in')
      ) + 
      coord_flip()
    pList[[goTerm]] <- p
    if (length(unique(indivPlotData[, 'sig'])) == 2 & is.null(pListLegend)){
      pLegend <- p + 
        geom_bar(aes(x = reorder(hgnc_symbol, - log2FoldChange), y = log2FoldChange, fill = sig), 
                 stat = 'identity', 
                 color = 'black', size = 0.25,
                 show.legend = TRUE) + 
        theme(legend.position = 'bottom', legend.text = element_text(size = 10), legend.title = element_text(size = 10))
      pListLegend <- getPlotLegend(pLegend)
    }
  }
  layoutMtx <- pageLayoutMtxList[[goCategory]]$pageLayout
  cairo_pdf(file = getOutFilename(saveTo = '~/projects/lab-notebook/figs', proj = 'noct', name = paste('geneBarCharts.', goCategory, sep = ''), suffix = 'pdf'), width = 2  + 4 * ncol(layoutMtx), height = 11)
  grid.arrange(
               grobs = unname(pList), 
               layout_matrix = layoutMtx, 
               widths = rep(unit(4, 'in'), ncol(layoutMtx)),
               heights = rep(unit(9 / 6, 'in'), nrow(layoutMtx)),
               top = textGrob(label = paste('Supplementary Figure S', as.character(1 + i), sep = ''), x = 0, hjust = 0, gp = gpar(fontface = 'bold')), 
               bottom = pListLegend
  )
  dev.off()
}


# # back-up for functions developed for exploratory plotting
alignTitles <- function(gplot, title = 4, subtitle = 4, caption = 4){
  # # grab the saved ggplot2 object
  g <- ggplotGrob(gplot)
  # # find the object which provides the plot information for title, subtitle, and caption
  g$layout[which(g$layout$name == 'title'), ]$l <- title
  g$layout[which(g$layout$name == 'subtitle'), ]$l <- subtitle
  g$layout[which(g$layout$name == 'caption'), ]$l <- caption
  return(g)
}
addSignificanceSign <- function(gplot, red = TRUE, blue = TRUE, titleText, fz = 7.5){
  g <- ggplotGrob(gplot)
  layoutCount <- nrow(g$layout)
  # # mimic grob of title
  titleCount <- which(g$layout$name == 'title')
  titleGrob <- g$grobs[[titleCount]]
  sigGrob <- titleGrob
  # # mimic layout of title
  g$layout[titleCount, c('l', 'r')] <- g$layout[titleCount, 'l'] - 1
  titleLayout <- g$layout[titleCount, ]
  sigLayout <- titleLayout
  sigLayout[, c('t', 'b')] <- titleLayout[, 't'] + 1
  # # adjust position of signs
  lineCount <- stringr::str_count(titleText, '\n') + 1
  sigGrob$children[[1]]$vjust <- 1 - 1 / lineCount
  # # adjust text to add in
  bottomText <- gsub(x = titleText, pattern = '.*\n', replacement = '', perl = TRUE)
  if (red) {
    layoutCount <- layoutCount + 1
    # # add layout
    sigLayout[, 'name'] <- 'redAsterisk'
    g$layout[layoutCount, ] <- sigLayout
    # # add grob
    sigGrob$children[[1]]$label <- bquote(paste(phantom(.(bottomText)), bold('\u2217')))
    sigGrob$children[[1]]$gp$col <- 'red'
    g$grobs[[layoutCount]] <- sigGrob
    bottomText <- paste(bottomText, '\u2217', sep = '')
  }
  if (blue) {
    layoutCount <- layoutCount + 1
    # # add layout
    sigLayout[, 'name'] <- 'bluePlus'
    g$layout[layoutCount, ] <- sigLayout
    # # add grob
    sigGrob$children[[1]]$label <- bquote(paste(phantom(.(bottomText)), bold('+')))
    sigGrob$children[[1]]$gp$col <- 'blue'
    g$grobs[[layoutCount]] <- sigGrob
  }
  invisible(g)
}

