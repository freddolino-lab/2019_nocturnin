# 2019 Nocturnin RNA-seq analysis code

## Overview

This is the repository holding all formatted scripts used in bioinformatics analysis for:

**Differential processing and localization of human Nocturnin controls metabolism of mRNA and nicotinamide dinucleotide metabolites**
Elizabeth T. Abshire, Kelsey L. Hughes, Rucheng Diao, Sarah Pearce, Raymond C. Trievel, Joanna Rorbach, Peter L. Freddolino, Aaron C. Goldstrohm
*bioRxiv* 2020.01.12.903534; doi: https://doi.org/10.1101/2020.01.12.903534

Detailed methods and statistics are described in the Materials and Methods in the paper. Briefly, raw RNA-seq data was quality controlled and sequencing reads were trimming accordingly. The data was then aligned to a custom reference genome containing human reference genome, ERCC, and overexpression components. The alignment was quantified and analyzed for differential expression and functional enrichment.

In all codes, the user name of submitter/analyzer has been replaced by string `user`, and the user name of the principle investigator has been replaced by string `PI`.

## List of scripts

+ `DESeq2.1.61.1.patch.R`: patched function applied to DESeq2 (version 1.16.1) to include `lfcSE` in the shrinkage output;
+ `LICENSE`
+ `noct-00-envSetup.sh`: paths to scripts and tools used in analysis, and environmental variable set-ups;
+ `noct-00-linkFastq.sh`: script to create symbolic links to raw sequencing data at location of data analysis directory;
+ `noct-01-qc.sh`: script to perform quality control on raw sequencing reads;
+ `noct-02-cutAndTrim.sh`
+ `noct-02-qc.sh`
+ `noct-03-pipeline-runStarHtseq.sh`: pipelien script running STAR and HTSeq for alignment and quantification by calling `noct-03-runStar.sh` and `noct-04-quantHTSeq.sh`; 
+ `noct-03-buildStarRefGTF.py`: Python3 script to build annotation (GTF files) and sequences (FASTA) information for ERCC spike-ins and components of over-expression vectors;
+ `noct-03-runStar.sh`: script to build STAR reference and run STAR alignment;
+ `noct-04-quantHTSeq.sh`: script to run HTSeq quantification;
+ `noct-05-deseq2.htseq.R`: script to analyze HTSeq quantification results with DESeq2 for differential expression;
+ `noct-05-ipage_db.R`: script to query and build database files used in iPAGE analysis;
+ `noct-05-plotVolcanoGoGene.R`: scirpt to plot volcano plot (Figure 5c) and GO term specific gene expression level bar charts (Figure S2);
+ `noct-05-targetConfIntv.R`: script to calculate target confidence intervals for Figure 5d;
+ `noct-06-ma.R`: script to plot MA plot (Figure 5b);
+ `noct-06-mitoCategory.R`: script to add mitochondria-related categorization to TPM table;
+ `noct-06-notchedBox.R`: script to plot notched boxplot for expression change levels of each mitochondria-related ccategory (Figure S3);
+ `noct-06-totalExonLength.R`: script to query and calculate the total exon length for each gene for normalization in TPM calculation (`noct-06-tpm.R`);
+ `noct-06-tpm.R`: calculate the TPM of gene expression (Supporting Information Table S1).
+ `README.md`
