# calculate the confidence intervals for log2FC for the targets of interest
# targets:
# FRRS1L
#   RORB
# PTPRZ1
# DLGAP1
# TMEFF2
#   NOCT
#    B2M
#   ACTB
#    TBP

source('~/projects/styles/R/general.R')

tbl <- read.table(file = '20191223-TPM.alignedBySTAR.rounded.DESeqStats.csv', sep = ',', header = TRUE, stringsAsFactors = FALSE)
# print(head(tbl))

genes <- c('FRRS1L', 'RORB', 'PTPRZ1', 'DLGAP1', 'TMEFF2', 'NOCT', 'B2M', 'ACTB', 'TBP')
oe <- c('pFC3F_GST', 'pFC3F_NOCT_del2-15')
fields <- c('hgnc_symbol', 'gene_id', 'padj', 'log2FoldChange', 'lfcSE', 'stat', 'pvalue')
# tblInterest <- tbl[tbl$hgnc_symbol %in% genes, ]
# print(tblInterest)
# print(tbl[match(table = tbl$hgnc_symbol, x = genes), ])
tblInterest <- tbl[match(table = tbl$hgnc_symbol, x = genes), fields]
tblOE <- tbl[tbl$gene_id %in% oe, fields]
tblInterest <- rbind(tblInterest, tblOE)
print(tblInterest)

rownames(tblInterest) <- NULL
# print(tblInterest)

z <- 1.96 # constant for using SE to calculate 95% confidence interval
tblInterest$CI_lo <- tblInterest$log2FoldChange - z * tblInterest$lfcSE
tblInterest$CI_hi <- tblInterest$log2FoldChange + z * tblInterest$lfcSE
print(tblInterest, digits = 4, sep = '\t')

write.table(tblInterest, 
			file = getOutFilename(saveTo = '~/data/noct/humanOE/06-tpm/', 
			                      name = 'target_CI', 
			                      suffix = 'csv', proj = 'noct'), 
			sep	= ',', quote = FALSE, 
			row.names = FALSE, col.names = TRUE)
# move on to Google Sheet for analysis combined with qPCR data and formatting
