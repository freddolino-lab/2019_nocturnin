#!/usr/bin/env bash

# Run STAR alignment, reference version 201806, 
# i.e. genome + ERCC + NOCT and GST overexpressed transcripts with UTRs

source noct-00-envSetup.sh
runThread=8

# build reference 'genome' for alignment
python3 ~/projects/noct/humanOE/noct-03-buildStarRefGTF.py
echo 'combining custom genome FASTA'
zcat "$humanRefFasta" > "$newRefName"'.fasta'
cat "$newRefName"'.suppl.fasta' >> "$newRefName"'.fasta' 

# build STAR index
echo 'running STAR indexing'
mkdir "$starGenomeDir"
overhangLen='100'
"$starPath" --runThreadN "$runThread" --runMode genomeGenerate --genomeDir "$starGenomeDir" --genomeFastaFiles "$newRefName"'.fasta' --sjdbGTFfile "$newRefName".gtf --sjdbOverhang "$overhangLen" # --limitGenomeGenerateRAM 51000000000

# running STAR
echo 'running STAR alignment'
mkdir "$starDir"
for filename in "$seqDir"*"$fastqR1Suffix";
do
    # echo "$filename"
    filenameBase=${filename#$seqDir}
    filenameBase=${filenameBase%$fastqR1Suffix}
    echo "$filenameBase"
    # mkdir "$starDir""$filenameBase"
    printf 'processing file: %s\n' "$filename"
    echo ' = STAR alignment = '
    # if using sub-directory
    # remember to add '/' at the end when using filenameBase as output prefix
    # i.e. --outFileNamePrefix "$starDir""$filenameBase"'/'
    "$starPath" --runThreadN "$runThread" --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir "$starGenomeDir" --readFilesIn "$cutDir""$filenameBase""$r1".paired"$trimSuffix" "$cutDir""$filenameBase""$r2".paired"$trimSuffix" --outFileNamePrefix "$starDir""$filenameBase"'.' >"$starDir""$filenameBase".log 2>"$starDir""$filenameBase".err
    # index output sorted BAM file
    echo ' = samtools indexing = '
    samtools index -@ "$runThread" "$starDir""$filenameBase"'.Aligned.sortedByCoord.out.bam'
done

# Next step: run HTseq after alignment
# bash ~/projects/noct/humanOE/noct-04-quanHTSeq.sh
