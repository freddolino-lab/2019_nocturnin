# Project/code: Noct, RNA-seq of NOCT overexpression
# Start date: Sun Oct 20 16:08:14 2019
# Objective: query and calculate total exon length for each ENSG gene from 
#            Ensembl to TPM calculation
# --------------
# Author: diaorch
# Modification date:  Sun Oct 20 16:08:14 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------

tpmPath <- '06-tpm/'

#### Get "gene length" information ####
# "gene length" is defined as the total length of exon(s) of a gene, which is 
# NOT: any one of the transcript length(s), or, start position - end postion 
# on chromosome (i.e. total lengths of gene)  
# alignment done on GRCh38
library('biomaRt')
listEnsembl() # current `Ensembl Genes` version is 98
ensembl98 = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
exonicInfo98 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_exon_id',
                                   'exon_chrom_start','exon_chrom_end'), 
                    # filters = 'ensembl_gene_id',
                    # values = c('ENSG00000000003', 'ENSG00000000005', 'ENSG00000000419'), 
                    mart = ensembl98, bmHeader = FALSE)
print(head(exonicInfo98))
# previous data analysis was done in version 92
ensembl92 = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', version = 92)
exonicInfo92 <- getBM(attributes = c('ensembl_gene_id', 'ensembl_exon_id',
                                     'exon_chrom_start','exon_chrom_end'), 
                      mart = ensembl92, bmHeader = FALSE)
print(head(exonicInfo92))
save(exonicInfo98, exonicInfo92,
     file = paste(tpmPath, '20191020-exonic_info.RData', sep = ''))

#### calculate exon-covered length (wrong method) ####
# BECAUSE exons can be overlapping, cannot add up all exon lengths as "total exon length"
# # exonicLength <- as.vector(x = exonicInfo$manual_exon_length <- exonicInfo$exon_chrom_end - exonicInfo$exon_chrom_start + 1, mode = 'integer')
# exonicLenSumPerGene <- aggregate(
#   x = exonicInfo$manual_exon_length, 
#   by = list(aggregate_gene_id = exonicInfo$ensembl_gene_id), FUN = sum)
# print(head(exonicLenSumPerGene)) 
# print(exonicInfo[exonicInfo$ensembl_gene_id == 'ENSG00000003989', ])
# print(exonicLenSumPerGene[exonicLenSumPerGene$aggregate_gene_id == 'ENSG00000003989', ])
# print(sum(exonicInfo$manual_exon_length[exonicInfo$ensembl_gene_id == 'ENSG00000003989']))

#### test data 1: manual example ####
# testData <- data.frame(
#   ensembl_gene_id = '0', 
#   ensembl_gene_id =  c('1', '2', '3', '4', '5', '6', '7', '8', '9'), 
#   exon_chrom_start = c(101, 120, 130, 191, 210, 220, 240, 240, 191), 
#   exon_chrom_end   = c(125, 150, 170, 240, 230, 230, 250, 245, 255)
# )
# testData$manual_exon_length <- testData$exon_chrom_end - testData$exon_chrom_start + 1
# print(testData)
# ense <- testData

#### test data 2: ENSG00000003989 total exon-covered length = 7905 ####
# ense <- exonicInfo[exonicInfo$ensembl_gene_id == 'ENSG00000003989', ]

#### calculate total exon length for all genes ####
exonicInfo <- exonicInfo92
exonicInfo$manual_exon_length <- exonicInfo$exon_chrom_end - exonicInfo$exon_chrom_start + 1
print(head(exonicInfo))

totalExonLenEnsg <- rep(NA, length(unique(exonicInfo$ensembl_gene_id)))
names(totalExonLenEnsg) <- unique(exonicInfo$ensembl_gene_id)
for (ensgAccn in unique(exonicInfo$ensembl_gene_id)){
  ense <- exonicInfo[exonicInfo$ensembl_gene_id == ensgAccn, ]
  ense <- ense[order(ense$exon_chrom_start), ]
  # print(head(ense))
  # print(ense)
  totalExonLen <- 0
  s <- -1 # start
  e <- -2 # end
  for (i in 1:nrow(ense)){
    # print(ense[i, ])
    if (ense$exon_chrom_start[i] > e) {
      totalExonLen <- totalExonLen + (e - s + 1)
      s <- ense$exon_chrom_start[i]
      e <- ense$exon_chrom_end[i]
      # print('move both')
    } else {
      if (ense$exon_chrom_start[i] <= e){
        # no move s
      }
      if (ense$exon_chrom_end[i] > e){
        e <- ense$exon_chrom_end[i]
        # print('move e')
      }
    }
    if (i == nrow(ense)){
      totalExonLen <- totalExonLen + (e - s + 1)
    }
  }
  # print(totalExonLen)
  totalExonLenEnsg[ensgAccn] <- totalExonLen
}
print(head(totalExonLenEnsg))
print(length(totalExonLenEnsg))
print(totalExonLenEnsg['ENSG00000003989'])
print(any(names(totalExonLenEnsg) == 'ENSG00000130489'))
# save(totalExonLenEnsg, 
#      file = paste(tpmPath, '20191020-exonic_total_length.Ensembl.RData', sep = ''))
write.table(x = totalExonLenEnsg, 
            file = paste(tpmPath, '20191020-exonic_total_length.Ensembl.csv', sep = ''), 
            sep = ',', col.names = FALSE, row.names = TRUE)

# add pFC3F-NOCT-del2-15 and pFC3F-GST
# exon length info pulled from Kelsey's info, corresponds with:
# ~/data/noct/humanOE/00-ref/customized/overexpression.gff3
totalExonLenPlasmid <- c((856 - 735 + 1) + (2064 - 990 + 1), 
                         (857 - 735 + 1) + (2661 - 990 + 1 ))
names(totalExonLenPlasmid) <- c('pFC3F_GST', 'pFC3F_NOCT_del2-15')
print(totalExonLenPlasmid)
totalExonLenAll <- c(totalExonLenPlasmid, totalExonLenEnsg)
print(head(totalExonLenAll))
save(totalExonLenEnsg, totalExonLenPlasmid, totalExonLenAll,
     file = paste(tpmPath, '20191020-exonic_total_length.RData', sep = ''))
write.table(x = totalExonLenAll, 
            file = paste(tpmPath, '20191020-exonic_total_length.customRef.csv', sep = ''), 
            sep = ',', col.names = FALSE, row.names = TRUE)
