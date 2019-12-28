#!/usr/env/bin bash 

# original user name of submitter has been replaced as "user"
# original user name of PI has been replaced as "PI"

wd="/home/user/data/noct/humanOE/"
projectDir='/home/user/projects/noct/humanOE/'

rawDataPath="/seqdata/Goldstrohm_Project/" # mock raw data path
rawSuffix="_001.fastq.gz"

# link fastq files
seqDir="$wd""00-seq/"
fastqSuffix=".fastq.gz"
r1="_R1"
r2="_R2"
fastqR1Suffix="$r1""$fastqSuffix"
fastqR2Suffix="$r2""$fastqSuffix"

# quality control on raw data
qcDir="$wd""01-qc/"

# fixed length cutting (Cutadapt)
cutDir="$wd""02-cutAndTrim/"
CUTADAPT='/home/PI/src/cutadapt-1.8.1/bin/cutadapt'
cutSuffix=".cut.fastq.gz"

# quality trimming (Trimmomatic)
TRIMMOMATICPATH='/home/PI/src/Trimmomatic-0.33/trimmomatic-0.33.jar' 
trimSuffix=".trim.fastq.gz"

# kallisto indexing
KALLISTOPATH='/home/PI/src/kallisto_linux-v0.43.0/kallisto'
refDir="$wd""00-ref/"

# kallisto quantification
# quantDir="$wd""03-quant/"
quantDir="$wd""03-04-kallisto/"

# STAR alignment, new transcript references
# 20180612
starPath='STAR'
starDir="$wd"'03-star/'
starGenomeDir="$refDir"'starRef/'
humanRefFasta='/srv/user/Homo_sapiens/Ensembl/GRCh38/Sequence/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
newRefName="$refDir"'customGenomePlasmid_061218'
starSortedBamSuffix='.Aligned.sortedByCoord.out.bam'

# HTSeq quantification
htPath='htseq-count'
htDir="$wd"'04-htseq/'
htSuffix='.stranded.htseq.count'

# calculate Transcript Per Million (TPM) values 20191020
tpmDir="$wd"'06-tpm/'
