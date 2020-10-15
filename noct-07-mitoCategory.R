# Project/code: Noct
# Start date: Fri Nov  8 14:48:26 2019
# Objective: format and plot TPM of human cell line overexpression data
# --------------
# Author: diaorch
# Modification date:  Fri Nov  8 14:48:26 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------
library('biomaRt')
ensembl92 = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', 
                       version = 92)
listAttributes(ensembl92) # 9 chromosome_name Chromosome/scaffold name feature_page

chrNameAttributes <- getBM(
  attributes = c('chromosome_name'), mart = ensembl92
)
print(dim(chrNameAttributes))
print(chrNameAttributes$chromosome_name)
# saw MT in chromosome names

mtGeneSet <- getBM(
  attributes = c('ensembl_gene_id', 'chromosome_name'),
  filters = 'chromosome_name' , values = 'MT', 
  mart = ensembl92)
print(head(mtGeneSet)) # mito-encoded genes
print(dim(mtGeneSet))
# https://www.ncbi.nlm.nih.gov/gene/?term=MT%5BCHR%5D+AND+human%5BORGN%5D
# if only count genes marked "[Homo sapiens (human)]", there are 37 too
# mtTranscriptSet <- getBM(
#   attributes = c('ensembl_transcript_id', 'chromosome_name'),
#   filters = 'chromosome_name' , values = 'MT', 
#   mart = ensembl92)
# print(head(mtTranscriptSet))
# rm(mtTranscriptSet)

# write.table(x = mtGeneSet, file = '06-tpm/20191108-mitoGeneTable.csv', 
#             sep = ',', quote = TRUE, col.names = TRUE, row.names = FALSE)
mtGeneSet <- read.table(file = '06-tpm/20191108-mitoGeneTable.csv', 
                        sep = ',', quote = '"', header = TRUE)
print(head(mtGeneSet))
print(dim(mtGeneSet))

tbl <- read.table(file = '06-tpm/20191223-TPM.alignedBySTAR.rounded.DESeqStats.csv', 
                  sep = ',', header = TRUE, stringsAsFactors = FALSE)
print(head(tbl))
print(tbl[startsWith(x = tbl$gene_id, prefix = 'p'), ])
print(tbl[tbl$hgnc_symbol == 'RORB', ])
print(nrow(tbl))
# print(sum(tbl$description == ''))
# print(head(tbl[tbl$description == '', ]))
# print(nrow(tbl[tbl$pvalue == '', ]))
# print(nrow(tbl[tbl$baseMean == '', ]))
print(tbl[!is.na(tbl$padj) & tbl$padj < 0.01, ])
print(nrow(tbl[!is.na(tbl$padj) & tbl$padj < 0.01, ]))
print(summary(tbl[!is.na(tbl$padj) & tbl$padj < 0.01, 'log2FoldChange']))
print(summary(tbl[!is.na(tbl$padj) & tbl$padj < 0.01, 'baseMean']))

deTbl <- tbl[!is.na(tbl$padj) & tbl$padj < 0.05, ]
print(head(deTbl))
deTbl[, 'TPM_mean'] <- (deTbl$Del_3F_11_S5_TPM + deTbl$Del_3F_12_S6_TPM + deTbl$Del_3F_4_S4_TPM + 
                          deTbl$GST_3F_2_S1_TPM + deTbl$GST_3F_8_S2_TPM + deTbl$GST_3F_9_S3_TPM) / 6
print(head(deTbl))
library('ggplot2')
library('gridExtra')
p1 <- ggplot(data = deTbl) + 
  geom_point(aes(x = TPM_mean, y = log2FoldChange), size = 0.5) + 
  geom_hline(yintercept = 0) + 
  labs(
    title = 'Expression changes vs expression level for significantly changed genes (adjusted p-value < 0.05):', 
    x = 'Mean TPMs of six samples') + 
  theme(
    title = element_text(size = 8.5)
  )

p1_transx <- ggplot(data = deTbl) + 
  geom_point(aes(x = TPM_mean, y = log2FoldChange), size = 0.5) + 
  geom_hline(yintercept = 0) + 
  scale_x_continuous(trans = 'log10') + 
  labs(title = 'With average TPM transformed by log10', 
       x = 'Mean TPMs of six samples') + 
  theme(
    title = element_text(size = 8.5)
  )

grid.arrange(p1, p1_transx, ncol = 1)

print(head(tbl))
goi <- c('NADK', 'NADK2', 'HDDC3') # NADK1, NADK2, MESH1
print(tbl[tbl$hgnc_symbol %in% goi, ])
goiTbl <- tbl[tbl$hgnc_symbol %in% goi, ]
goiTbl[, 'TPM_mean'] <- (goiTbl$Del_3F_11_S5_TPM + goiTbl$Del_3F_12_S6_TPM + goiTbl$Del_3F_4_S4_TPM + 
                          goiTbl$GST_3F_2_S1_TPM + goiTbl$GST_3F_8_S2_TPM + goiTbl$GST_3F_9_S3_TPM) / 6
print(goiTbl)

# email 20191101 "nuclear encoded"
# https://www.broadinstitute.org/files/shared/metabolism/mitocarta/human.mitocarta2.0.html
# mitoCarta <- read.table(file = '06-tpm/Human.MitoCarta2.0.col7.txt', 
#                         header = TRUE, sep = '\t') # all genes
# the following mitoCarta data is the sheet A of supplied .xlsx file
# "a collection of 1158 nuclear and mtDNA genes encoding proteins with strong support of mitochondrial localization"
# according to broad inistitute website above
mitoCarta <- read.table(file = '06-tpm/Human.MitoCarta2.0.SheetA.1158genes.csv', 
                        header = TRUE, quote = '"', sep = ',', stringsAsFactors = FALSE)
print(head(mitoCarta))
print(dim(mitoCarta))
print(tail(mitoCarta))
print(head(mitoCarta$EnsemblGeneID))

print(head(mtGeneSet))
print(mitoCarta[mitoCarta$hg19_Chromosome == 'chrM', ])
print(mtGeneSet$ensembl_gene_id[!(mtGeneSet$ensembl_gene_id %in% mitoCarta$EnsemblGeneID)])
print(dim(mtGeneSet))
print(dim(mitoCarta))
print(sum(mtGeneSet$ensembl_gene_id %in% mitoCarta$EnsemblGeneID))
print(sum(!(mtGeneSet$ensembl_gene_id %in% mitoCarta$EnsemblGeneID)))

allGeneCartaTbl <- read.table(file = '06-tpm/Human.corresponding.allGenes.csv', 
                              header = TRUE, quote = '"', sep = ',', stringsAsFactors = FALSE)
print(head(allGeneCartaTbl))
print(dim(allGeneCartaTbl))
print(head(allGeneCartaTbl$EnsemblGeneID))

nucMtTbl <- tbl
nucMtTbl[, 'category'] <- NA
for(i in 1:nrow(nucMtTbl)){
  # i <- 1
  tag = ''
  if (nucMtTbl$gene_id[i] %in% mtGeneSet$ensembl_gene_id){
    tag = 'Mitochondrially-encoded gene'
  } else if (!(nucMtTbl$gene_id[i] %in% mtGeneSet$ensembl_gene_id) & # not mt-encoded
             nucMtTbl$gene_id[i] %in% mitoCarta$EnsemblGeneID){ # high confidence mt-localization
    tag = 'Nuclear-encoded with high-confidence mitochondrially localized gene'
  } else if (!(nucMtTbl$gene_id[i] %in% mtGeneSet$ensembl_gene_id) & # not mt-encoded
             !(nucMtTbl$gene_id[i] %in% mitoCarta$EnsemblGeneID) & # not high confidence mt-localization
             nucMtTbl$gene_id[i] %in% allGeneCartaTbl$EnsemblGeneID){ # in the rest of Human Mito Carta data
    tag = 'Nuclear-encoded with low-confidence mitochondrially localized gene'
  } else{
    tag <- 'Unknown'
  }
  nucMtTbl$category[i] <- tag
}
print(table(nucMtTbl$category))
print(dim(mtGeneSet))
print(dim(mitoCarta))
print(dim(allGeneCartaTbl))
print(nucMtTbl[sample(x = 1:nrow(nucMtTbl), size = 5, replace = FALSE), ])

source('~/projects/styles/R/general.R')
write.table(x = nucMtTbl, file = getOutFilename(name = 'categorizedTable.addUnknown.TPM.DE', saveTo = '06-tpm', 
                                                proj = 'noct', suffix = '.csv'), 
            sep = ',', quote = TRUE, col.names = TRUE, row.names = FALSE)
