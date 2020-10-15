# # summarize iPAGE output, label GO terms of interest with categories, summrize as
# # supplementary table
# rm(list = ls())
load(paste(Sys.getenv('HOME'), '/data/noct/humanOE/revision/20200716-noct-goToGene.permTest.summarized.RData', sep = ''))
print(ls())

# # gene stats associated with GO terms
print(head(allGoEnsg))
pvmtx <- read.table(
  file = paste(Sys.getenv('HOME'), '/data/noct/humanOE/revision/z-iPAGE/20200605-noct-z_shrunken_input.tsv.sigBox.indp0_PAGE/pvmatrix.txt', sep = ''), sep = '\t', header = FALSE, skip = 1, stringsAsFactors = FALSE)
print(head(pvmtx))
getGOaccn <- function(x){
  return(strsplit(x = x, split = ' ', fixed = TRUE)[[1]][1])
}
pvmtxAccn <- sapply(X = pvmtx[, 'V1'], FUN = getGOaccn)
# print(head(pvmtxAccn))

# # need: median(t-stats), median(l2fc), median(abs(t-stats)), median(abs(l2fc))
geneStatsTbl <- allGoEnsg[, c('ENSG', 'log2FoldChange', 'lfcSE', 'GO')]
geneStatsTbl[, 't'] <- geneStatsTbl[, 'log2FoldChange'] / geneStatsTbl[, 'lfcSE']
# print(head(geneStatsTbl))
resTbl <- as.data.frame(matrix(data = NA, nrow = length(pvmtxAccn), ncol = 7, dimnames = NULL), stringsAsFactors = FALSE)
resTbl[, 1] <- unname(pvmtxAccn)
resTbl[, 2] <- pvmtx[, 'V1']
# print(str(resTbl))
for (i in 1:nrow(resTbl)){
  g <- pvmtxAccn[i]
  print(i)
  indivStatsTbl <- geneStatsTbl[geneStatsTbl[, 'GO'] == g, ]
  resTbl[i, 3] <- median(indivStatsTbl[, 't'])
  resTbl[i, 4] <- median(indivStatsTbl[, 'log2FoldChange'])
  resTbl[i, 5] <- median(abs(indivStatsTbl[, 't']))
  resTbl[i, 6] <- median(abs(indivStatsTbl[, 'log2FoldChange']))
}

goBone  <- c('GO:0045668', 'GO:0005518', 'GO:0002062')
goMito  <- c('GO:0005747', 'GO:0045333', 'GO:0007005', 'GO:0006631', 'GO:0070125', 'GO:0070469')
goNAD   <- c('GO:0008137', 'GO:0006120', 'GO:0051287')
goNeuro <- c('GO:0007218', 'GO:0005267', 'GO:0005244', 'GO:0051402', 'GO:0001580',
             'GO:0050909', 'GO:0071805', 'GO:0001843', 'GO:0051928', 'GO:0005262',
             'GO:0042053')
goBoneFat <- c('GO:0070328')
goCheckRevision <- c(goBone, goMito, goNAD, goNeuro, goBoneFat)
goCheckCategory <- c(rep('Bone', length(goBone)), rep('Mito', length(goMito)), rep('NAD', length(goNAD)), rep('Neuro', length(goNeuro)), rep('BoneFat', length(goBoneFat)))
names(goCheckCategory) <- goCheckRevision
# print(goCheckCategory)
resTbl[, 7] <- goCheckCategory[resTbl[, 1]]
resTbl[is.na(resTbl[, 7]), 7] <- 'Other'
# print(head(resTbl))
resTbl <- resTbl[order(resTbl[, 7], resTbl[, 3], decreasing = FALSE), ]
colnames(resTbl) <- c('GO Acccession', 'GO Description', 
                      'Median of t-statistic', 'Median of Log2 Fold Change', 
                      'Median of Abs Value of t-statistic', 'Median of Abs Value of Log2 Fold Change', 
                      'Category')
resTbl <- resTbl[, c(7, 1:6)]
write.table(resTbl, 
            file = getOutFilename(saveTo = paste(Sys.getenv('HOME'), '/projects/lab-notebook/tbls/', sep = ''), name = 'ipageTerms', 
                                  proj = 'noct', suffix = '.csv'), 
            sep = ',', row.names = FALSE, col.names = TRUE, quote = TRUE)
