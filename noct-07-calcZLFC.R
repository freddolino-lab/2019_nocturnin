#!/usr/bin/env Rscript
# # calculate the z-scores of log2 fold changes from DESeq2 to use as a more robust
# # iPAGE input

dataPath<- '~/data/noct/humanOE/revision/z-iPAGE/'
deseqResFilename <- paste(dataPath, '20191220-noct-patched.newDESeq2Shrinkage.csv', sep = '')
deseqResTbl <- read.table(deseqResFilename, header = TRUE, sep = ',', quote = '"', stringsAsFactors = FALSE)
# print(head(deseqResTbl))

modTbl <- deseqResTbl # modified table
modTbl[, 'stat_shrunken'] <- modTbl[, 'log2FoldChange'] / modTbl[, 'lfcSE']
# print(head(modTbl[, c(1:5, 8, 6:7)]))

source('~/projects/styles/R/general.R')
write.table(modTbl[, c(1:5, 8, 6:7)], getOutFilename(saveTo = dataPath, proj = 'noct', name = 'z_shrunken_full', suffix = '.csv'), sep = ',', row.names = FALSE, col.names = TRUE, quote = TRUE)
write.table(modTbl[, c('X', 'stat_shrunken')], getOutFilename(saveTo = dataPath, proj = 'noct', name = 'z_shrunken_input', suffix = '.tsv'), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
