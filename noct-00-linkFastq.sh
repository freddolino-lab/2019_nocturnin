#!/usr/bin/env bash

# rawDataPath="/seqdata/Goldstrohm_Project/"
# rawSuffix="_001.fastq.gz"

# WD="/home/user/data/noct/humanOE/"
# seqDir="$WD""00-seq/"
# fastqSuffix=".fastq.gz"

mkdir "$seqDir"

for filename in "$rawDataPath"*"$rawSuffix"; do
  # echo "$filename"
  basename=${filename#$rawDataPath}
  basename=${basename%$rawSuffix}
  # echo "$basename"
  ln -s "$filename" "$seqDir""$basename""$fastqSuffix"
done
