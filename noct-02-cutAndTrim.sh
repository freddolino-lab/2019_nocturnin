#!/usr/bin/env bash

# 20180430 Nocturnin Human cellline OE RNA-seq 
# this script is to perform read cutting (fixed lengths or adapter sequences) and trimming (quality), based on the inspection result from raw data quality control

source noct-00-envSetup.sh

NTHREAD=8
trimMinLen=26

mkdir "$cutDir"

for filename in "$seqDir"*"$fastqR1Suffix"; do
  # echo "$filename"
  filenameBase=${filename#$seqDir}
  filenameBase=${filenameBase%$fastqR1Suffix}
  echo "$filenameBase"
  printf 'processing file: %s\n' "$filename"

  # remove the first 10 bases of reads
  printf '+ cutting first 10 bases off Reads 1\n'
  "$CUTADAPT" -u 10 -o "$cutDir""$filenameBase""$r1""$cutSuffix" "$seqDir""$filenameBase""$r1""$fastqSuffix"
  printf '+ cutting first 10 bases off Reads 2\n'
  "$CUTADAPT" -u 10 -o "$cutDir""$filenameBase""$r2""$cutSuffix" "$seqDir""$filenameBase""$r2""$fastqSuffix"
  
  # quality trim
  printf '+ trimming by quality in paired-end mode\n'
  java -jar "$TRIMMOMATICPATH" PE -threads "$NTHREAD" -phred33 "$cutDir""$filenameBase""$r1""$cutSuffix" "$cutDir""$filenameBase""$r2""$cutSuffix" "$cutDir""$filenameBase""$r1"".paired""$trimSuffix" "$cutDir""$filenameBase""$r1"".unpaired""$trimSuffix" "$cutDir""$filenameBase""$r2"".paired""$trimSuffix" "$cutDir""$filenameBase""$r2"".unpaired""$trimSuffix" LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:"$trimMinLen"

done
