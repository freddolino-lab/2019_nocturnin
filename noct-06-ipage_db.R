#!/usr/bin/env Rscript

# Project/code: Nocturnin Overexpression RNA-seq TPM analysis for DESeq2 results
# Start date: Tue Aug 21 15:29:34 2018
# --------------
# Author: diaorch
# Date: 20180821
# Modification (if modified from another procjet/script):
#   0. /home/diaorch/project/noct/humanOE/noct-05-deseq2.R
#   1. /home/diaorch/project/noct/humanOE/noct-05-ipage_db.R
# --------------

# iPAGE database 
library(biomaRt)
ensembl = useEnsembl(biomart = 'ensembl', dataset = 'hsapiens_gene_ensembl')
listAttributes(mart = ensembl, page = 'feature_page')
# 39                             go_id                          GO term accession feature_page
# 40                         name_1006                               GO term name feature_page
# 41                   definition_1006                         GO term definition feature_page
# 42                   go_linkage_type                      GO term evidence code feature_page
# 43                    namespace_1003                                  GO domain feature_page
# getBM(attributes = c('ensembl_gene_id', 'go_id'), 
#       filters = 'ensembl_gene_id', 
#       values = sigGeneTblLite$ensemblGeneId[1:5], 
#       mart = ensembl, bmHeader = FALSE)

# getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id',
#                      'go_id'), 
#       filters = 'ensembl_gene_id', 
#       values = sigGeneTblLite$ensemblGeneId[1:5], 
#       mart = ensembl, bmHeader = FALSE)
# 
# getBM(attributes = c('ensembl_gene_id', 'go_id', 'name_1006', 'namespace_1003'), 
#       filters = 'ensembl_gene_id', 
#       values = sigGeneTblLite$ensemblGeneId[1:5], 
#       mart = ensembl, bmHeader = FALSE)

all_gene_go <- getBM(attributes = c('ensembl_gene_id', 'go_id'), 
                     mart = ensembl, bmHeader = FALSE)
dim(all_gene_go)
head(all_gene_go, n = 30)
all_go <- getBM(attributes = c('go_id', 'name_1006', 'namespace_1003'), 
                mart = ensembl, bmHeader = FALSE)
dim(all_go)
head(all_go)
save(all_gene_go, all_go, file = '~/data/noct/humanOE/20180827-noct-iPAGE_database_biomart_ensembl_hsapiens.RData')

### process gene-to-go into iPAGE database's *_index.txt file
dim(all_gene_go)
all_gene_go_complete <- all_gene_go[all_gene_go$go_id != '', ]
dim(all_gene_go_complete)
head(all_gene_go_complete)
head(all_gene_go, n = 20)
# test <- all_gene_go_complete[1:50, ]
# test_go <- aggregate(x = test$go_id, by = list(gene = test$ensembl_gene_id), FUN = function(x) base::paste(x, collapse = '\t'))
# # 1  ENSG00000194717
# # 1  GO:1903231\tGO:0035195\tGO:0010629\tGO:0090051\tGO:1903588\tGO:1905652\tGO:0048147
# print(test)
# # 87  ENSG00000194717 GO:1903231
# # 88  ENSG00000194717 GO:0035195
# # 89  ENSG00000194717 GO:0010629
# # 90  ENSG00000194717 GO:0090051
# # 91  ENSG00000194717 GO:1903588
# # 92  ENSG00000194717 GO:1905652
# # 93  ENSG00000194717 GO:0048147
# test_entry <- apply(X = test_go, MARGIN = 1, FUN = paste, collapse = '\t')
# tmpFilename <- '~/tmp/20180827-noct-ipageDbTest_index.txt'
# con <- file(tmpFilename, 'w')
# writeLines(test_entry, con = con)
# close(con)
collapsed_go <- aggregate(x = all_gene_go_complete$go_id, by = list(gene = all_gene_go_complete$ensembl_gene_id), FUN = function(x) base::paste(x, collapse = '\t'))
dim(collapsed_go)
length(unique(all_gene_go_complete$ensembl_gene_id))
head(collapsed_go)
collapsed_entry <- apply(X = collapsed_go, MARGIN = 1, FUN = paste, collapse = '\t')
head(collapsed_entry)
writeFilename <- '~/data/noct/humanOE/20180827-noct-iPAGE_database_biomart_ensembl_hsapiens_index.txt'
con <- file(writeFilename, 'w')
writeLines(collapsed_entry, con = con)
close(con)

### process go information into iPAGE database's *_names.txt file
dim(all_go)
head(all_go)
all_go_marked <- all_go
all_go_marked[all_go_marked == ''] <- NA
complete_go <- all_go_marked[complete.cases(all_go_marked), ]
head(complete_go)
dim(complete_go)
print(all_go_marked[!complete.cases(all_go_marked), ])
abbr_match <- c('F', 'C', 'P')
names(abbr_match) <- c('molecular_function', 'cellular_component', 'biological_process')
# + F: molecular function
# + C: cellular component
# + P: biological process
domain_namespace_abbr <- abbr_match[complete_go$namespace_1003]
length(domain_namespace_abbr)
head(domain_namespace_abbr)
head(complete_go)
go_summary <- cbind(complete_go$go_id, complete_go$name_1006, domain_namespace_abbr)
colnames(go_summary) <- NULL
rownames(go_summary) <- NULL
head(go_summary)
write.table(x = go_summary, file = '~/data/noct/humanOE/20180827-noct-iPAGE_database_biomart_ensembl_hsapiens_names.txt', quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)

saveit <- function(..., file) {
  x <- list(...)
  # print(str(x))
  save(list=names(x), file=file, envir=list2env(x))
}
saveit(hsapiens_index = collapsed_entry, hsapiens_names = go_summary, file = '~/data/noct/humanOE/20180827-noct-iPAGE_database_biomart_ensembl_hsapiens.RData')
load(file = '~/data/noct/humanOE/20180827-noct-iPAGE_database_biomart_ensembl_hsapiens.RData')
dim(hsapiens_index)
head(hsapiens_index)
dim(hsapiens_names)
head(hsapiens_names)
