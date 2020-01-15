# 2019 Nocturnin RNA-seq analysis code

## Overview

This is the repository holding all formatted scripts used in bioinformatics analysis for:
Differential processing and localization of human Nocturnin controls metabolism of mRNA and nicotinamide dinucleotide metabolites
Elizabeth T. Abshire, Kelsey L. Hughes, Rucheng Diao, Sarah Pearce, Raymond C. Trievel, Joanna Rorbach, Peter L. Freddolino, Aaron C. Goldstrohm
bioRxiv 2020.01.12.903534; doi: https://doi.org/10.1101/2020.01.12.903534

Detail methods and statistics are described in the Materials and Methods in the paper. Briefly, raw RNA-seq data was quality controlled and sequencing reads were trimming accordingly. The data was then aligned to a custom reference genome containing human reference genome, ERCC, and overexpression components. The alignment was quantified and analyzed for differential expression and functional enrichment.

In all codes, the user name of submitter/analyzer has been replaced by string `user`, and the user name of the principle investigator has been replaced by string `PI`.

## List of scripts
+ `DESeq2.1.61.1.patch.R`
+ `LICENSE`
+ `noct-00-envSetup.sh`
+ `noct-00-linkFastq.sh`
+ `noct-00-pipeline-runStarHtseq.sh`
+ `noct-01-qc.sh`
+ `noct-02-cutAndTrim.sh`
+ `noct-02-qc.sh`
+ `noct-03-buildStarRefGTF.py`
+ `noct-03-runStar.sh`
+ `noct-04-quantHTSeq.sh`
+ `noct-05-deseq2.htseq.R`
+ `noct-05-ipage_db.R`
+ `noct-05-plotVolcanoGoGene.R`
+ `noct-05-targetConfIntv.R`
+ `noct-06-ma.R`
+ `noct-06-mitoCategory.R`
+ `noct-06-notchedBox.R`
+ `noct-06-totalExonLength.R`
+ `noct-06-tpm.R`
+ `README.md`
