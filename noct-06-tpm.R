# Project/code: Noct, RNA-seq of NOCT overexpression
# Start date: Sun Oct 20 15:44:36 2019
# Objective: calculate TPM for genes for each sample
# --------------
# Author: diaorch
# Modification date:  Sun Oct 20 15:44:36 2019
# Modification (if modified from another procjet/script):
#   0. original (not modified from other sources)
# --------------

# setwd('~/data/noct/humanOE/06-tpm/')

tpmPath <- '~/data/noct/humanOE/06-tpm/'
fragLenSuffix <- '.fistSegment.field9.csv'
filenameSet <- list.files(path = tpmPath, pattern = fragLenSuffix)
# print(filenameSet)

exonLengthFile <- '20191020-exonic_total_length.RData'

countPath <- '04-htseq/'
htseqSuffix <- '.stranded.htseq.count'
countFilenameSet <- list.files(path = countPath, pattern = htseqSuffix)
print(countFilenameSet)

deseqPath <- '05-deseq/' # path masked for code sharing
# deseqFile <- '20191001-volcano_plotting.RData'
deseqFile <- '20191223-noct-patched.deseq2_lfcShrink_ordered_table.csv'
deResTbl <- read.table(file = paste(deseqPath, deseqFile, sep = ''), 
                       header = TRUE, sep = ',', quote = '"', stringsAsFactors = FALSE)
print(head(deResTbl))

#### Calculate median fragment length ####

medianFragLen <- rep(x = NA, length(filenameSet))
names(medianFragLen) <- gsub(pattern = fragLenSuffix, replacement = '', x = filenameSet)
print(medianFragLen)

for (i in 1:length(filenameSet)){
  # i <- 1
  sampleName <- gsub(pattern = fragLenSuffix, replacement = '', x = filenameSet[i])
  
  fragLenFirstSegment <- read.table(
    file = paste(tpmPath, filenameSet[i], sep = ''), sep = ',', header = FALSE)
  # print(head(fragLenFirstSegment))
  
  fragLenSet <- abs(fragLenFirstSegment$V1)
  print(summary(fragLenSet))
  
  print(median(fragLenSet))
  medianFragLen[sampleName] <- median(fragLenSet)
}

save(medianFragLen, 
     file = paste(tpmPath, '20191020-median_fragment_lengths.RData', sep = ''))
# load(file = paste(tpmPath, '20191020-median_fragment_lengths.RData', sep = ''))
# print(medianFragLen)

#### Calculate TPM ####
# load median read length
load(paste(tpmPath, '20191020-median_fragment_lengths.RData', sep = ''))
# load gene exon-covered lengths
load(paste(tpmPath, exonLengthFile, sep = ''))
# load('20191020-exonic_info.RData')
# print(head(exonicInfo))
# print(any(exonicInfo92$ensembl_gene_id == 'ENSG00000130489'))

htseqAnnot <- c('__no_feature', '__ambiguous', '__too_low_aQual', 
                '__not_aligned', '__alignment_not_unique')

tpmSamples <- NULL

for (countFilename in countFilenameSet){
  # countFilename <- countFilenameSet[1]
  # print(countFilename)
  sampleName <- gsub(pattern = htseqSuffix, replacement = '', x = countFilename)
  
  fragLen <- medianFragLen[sampleName]
  print(fragLen)
  
  countTbl <- read.table(file = paste(countPath, countFilename, sep = ''), 
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE)
  colnames(countTbl) <- c('gene_id', 'htseq_count')
  # print(head(countTbl))
  # print(tail(countTbl))
  
  # print(dim(countTbl))
  countTbl$exon_length <- totalExonLenAll[countTbl$gene_id]
  # print(head(countTbl))
  # print(tail(countTbl))
  # print(nrow(countTbl))
  # print(sum(is.na(countTbl$exon_length)))
  # print(countTbl[is.na(countTbl$exon_length), ])
  countTblClean <- countTbl[!is.na(countTbl$exon_length), ]
  # print(countTbl[is.na(countTbl$exon_length), ])
  
  countTblClean$transcript_quant <- 
    countTblClean$htseq_count * fragLen / countTblClean$exon_length
  # print(head(countTblClean))
  
  totalTranscript <- sum(countTblClean$transcript_quant)
  sampleNameTpm <- paste(sampleName, '_TPM', sep = '')
  countTblClean[sampleNameTpm] <- 
    countTblClean$transcript_quant / totalTranscript * 10 ^ 6
  # print(head(countTblClean))
  
  # tpmSamples[[sampleName]] <- countTblClean$tpm
  if (is.null(tpmSamples)){
    tpmSamples <- data.frame(countTblClean[, c('gene_id', sampleNameTpm)])
  } else {
    tpmSamples <- merge(x = tpmSamples, 
                        y = countTblClean[, c('gene_id', sampleNameTpm)], 
                        by = 'gene_id', all = TRUE)
  }
  # print(head(tpmSamples))
}

print(dim(tpmSamples))
print(head(tpmSamples))

# formating
tpmSamplesRes <- format(tpmSamples, digit = 2, scientific = FALSE, trim = TRUE)
print(head(tpmSamplesRes))
tpmSamplesRes <- cbind(
  data.frame(gene_id = tpmSamplesRes$gene_id), 
  as.data.frame(sapply(X = tpmSamplesRes[2:ncol(tpmSamplesRes)], FUN = as.numeric)))
print(sapply(X = as.data.frame(tpmSamplesRes), FUN = mode))
print(head(tpmSamplesRes))
print(head(tpmSamplesRes$Del_3F_11_S5[1:6]))

# write.csv(x = tpmSamplesRes, 
#           file = paste(tpmPath, '20191020-TPM.alignedBySTAR.csv', sep = ''), 
#           # col.names = TRUE, 
#           row.names = FALSE, quote = FALSE)

#### annotate TPM table with gene info #####
# listMarts()
# ensembl <- useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
# listAttributes(mart = ensembl, page = 'feature_page', what = 'description')
# head(listAttributes(mart = ensembl, page = 'feature_page'), 20)
ensembl92 = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl', 
                       version = 92)
gnMapping <- getBM(
  attributes = c('ensembl_gene_id', 'hgnc_symbol', 'description'),
  filters = 'ensembl_gene_id' , values = tpmSamplesRes$gene_id, 
  mart = ensembl92)
print(head(gnMapping))

tpmSamplesRes <- merge(x = tpmSamplesRes, y = gnMapping, 
                       by.x = 'gene_id', by.y = 'ensembl_gene_id', 
                       all.x = TRUE, all.y = FALSE)
print(head(tpmSamplesRes))
print(tpmSamplesRes[1300:1305, ])
# ENSG00000230417 has two HGNC symbols

write.csv(x = merge(x = tpmSamples, y = gnMapping, 
                    by.x = 'gene_id', by.y = 'ensembl_gene_id', 
                    all.x = TRUE, all.y = FALSE), 
          # x = tpmSamplesRes,
          file = paste(tpmPath, '20191108-TPM.alignedBySTAR.csv', sep = ''), 
          # col.names = TRUE, 
          row.names = FALSE, quote = TRUE)
write.csv(x = tpmSamplesRes,
          file = paste(tpmPath, '20191108-TPM.alignedBySTAR.rounded.csv', sep = ''), 
          # col.names = TRUE, 
          row.names = FALSE, quote = TRUE)

#### DESeq info ####
# load(paste(deseqPath, deseqFile, sep = ''))
print(dim(deResTbl))
print(head(deResTbl))

tpmWirttenRes <- read.table(file = paste(tpmPath, '20191108-TPM.alignedBySTAR.rounded.csv', sep = ''), 
                            header = TRUE, sep = ',', stringsAsFactors = FALSE)
print(head(tpmWirttenRes))

tpmDeRes <- merge(x = tpmWirttenRes, y = deResTbl, 
                  by.x = 'gene_id', by.y = 'X', all.x = TRUE, all.y = FALSE)
print(head(tpmDeRes))
print(dim(tpmDeRes))
tpmDeRes <- tpmDeRes[order(tpmDeRes$gene_id), ]
print(head(tpmDeRes))
write.csv(x = tpmDeRes,
          file = paste(tpmPath, '20191223-TPM.alignedBySTAR.rounded.DESeqStats.csv', sep = ''), 
          # col.names = TRUE, 
          row.names = FALSE, quote = TRUE, na = '')
