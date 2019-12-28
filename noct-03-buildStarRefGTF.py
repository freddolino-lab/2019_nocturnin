#!/usr/env/bin python3

"""
build custom reference with plasmid and spike-in RNA sequences, for Nocturnin 
overexpression sequencing analysis. The custom reference includes:
New custom reference GTF: 
1) pFC3F GST plasmid exon-intron features, from customized file: overexpression.gtf
2) pFC3F NOCT del 2-15 plasmid exon-intron features, from customized file: overexpression.gtf
3) ERCC features 
4) Homo_sapiens.GRCh38.92 entries
5) overexpression.gff: manual inference of exon and intron status for 
   plasmid pFC3F-GST and pFC3F-NOCT-del-2-15
New custom reference FASTA SUPPLEMENTARY:
1) pFC3F GST plasmid sequence 
2) pFC3F NOCT del 2-15 plasmid sequence 
3) ERCC sequences
Use plasmid sequences from Kelsey on 20180612, see hard-coded inputs. 

USAGE:
source activate bio_env_py.3.4 # activate Anaconda environment with BioPython 
python3 noct-03-buildStarRefGTF.py

INPUT (HARD-CODED):
(each in specific path, see input in main)
pFC3F_GST_061218.gb
pFC3F_NOCT_del2-15_061218.gb
overexpression.gtf
ercc.fasta
Homo_sapiens.GRCh38.92.gtf.gz

OUTPUT (HARD-CODED): 
(in reference directory, see code)
customGenomePlasmid_061218.gtf
customGenomePlasmid_061218.suppl.fasta

NOTES:
GTF file comes from Ensemble. 
"""

from Bio import SeqIO
import fileinput
import gzip
import itertools

def detectFileFormat(filename):
    """
    return format name by examining the last suffix of file name
    """
    suffixFormatDict = {
                        'gb': 'GenBank', 
                        'fa': 'FASTA', 
                        'fasta': 'FASTA', 
                        'gff3': 'GFF3'}
    suffix = filename.split('.')[-1]
    return(suffixFormatDict[suffix])

def stripFilenameSuffix(filename, sep = '.'):
    """
    return the base name i.e. the file name stripped of suffix 
    (stripping only the last suffix after separator, default as '.')
    """
    return(filename.rsplit(sep, 1)[0])


def convertFastaSeqRecordToGTF(faRecord, sourceOfFile):
    """
    manually parse SeqRecord object (from FASTA format) to convert to GTF entries
    for each FASTA entry, 3 GTF lines are generated, each for gene, transcript, or exon
    level for the ERCC
    result is written as a list with only 1 string which contains all three lines, 
    with separating newlines
    """
    gtfGene = faRecord.id + '\t' + sourceOfFile + '\tgene\t1\t' + \
              str(len(record.seq) + 1) + '\t.\t+\t.\tgene_id \"' + faRecord.id + \
              '\"; gene_version \"0\"; gene_name \"' + faRecord.id + \
              '\"; gene_source \"' + sourceOfFile + \
              '\"; gene_biotype \"transcribed_spike_in\";'
    gtfTspt = faRecord.id + '\t' + sourceOfFile + '\ttranscript\t1\t' + \
              str(len(record.seq) + 1) + '\t.\t+\t.\tgene_id \"' + faRecord.id + \
              '\"; gene_version \"0\"; transcript_id \"' + faRecord.id + \
              '\"; transcript_version \"0' + \
              '\"; gene_name \"' + faRecord.id + '\"; gene_source \"' + sourceOfFile + \
              '\"; gene_biotype \"transcribed_spike_in'+ \
              '\"; transcript_name \"' + faRecord.id + \
              '\"; transcript_source \"' + sourceOfFile + \
              '\"; transcript_biotype \"spike_in_transcript' + \
              '\"; tag \"custom\"; transcript_support_level \"1\";'
    gtfExon = faRecord.id + '\t' + sourceOfFile + '\texon\t1\t' + \
              str(len(record.seq) + 1) + '\t.\t+\t.\tgene_id \"' + faRecord.id + \
              '\"; gene_version \"0\"; transcript_id \"' + faRecord.id + \
              '\"; transcript_version \"0\"; exon_number \"1' + \
              '\"; gene_name \"' + faRecord.id + '\"; gene_source \"' + sourceOfFile + \
              '\"; gene_biotype \"transcribed_spike_in'+ \
              '\"; transcript_name \"' + faRecord.id + \
              '\"; transcript_source \"' + sourceOfFile + \
              '\"; transcript_biotype \"spike_in_transcript' + \
              '\"; exon_id \"' + faRecord.id + '-01' + \
              '\"; exon_version \"0\"; tag \"custom\"; transcript_support_level \"1\";'
    gtfEntry = '\n'.join([gtfGene, gtfTspt, gtfExon])
    return([gtfEntry])

if __name__ == '__main__':
    refPath = '/home/user/data/noct/humanOE/00-ref/'
    gtfgzGRCh38Filename = ('/srv/user/Homo_sapiens/Ensembl/GRCh38/' 
                           'Annotation/Homo_sapiens.GRCh38.92.gtf.gz')
    labPathSet = [
                  'fromGoldstrohmLab/',  
                  'fromGoldstrohmLab/',  
                  'fromFreddolinoLab/',  
                  'customized/']  
    inputFilenameSet = [
                        'pFC3F_GST_061218.gb', 
                        'pFC3F_NOCT_del2-15_061218.gb',
                        'ercc.fasta',  
                        'overexpression.gtf']  
    inputFormatSet = [
                      'GenBank',
                      'GenBank',
                      'FASTA',  
                      'GTF']  
    biopyParamDict = {
                      'GenBank': 'genbank',
                      'FASTA': 'fasta',
                      'GFF3': 'gff3' }
    # outputFileFasta = open(refPath + outputFileBase + '.suppl.fasta', 'w')
    outputFastaRecords = None
    # outputFastaRecords =[] 
    
    # GST plasmid
    ## Sequence comes from `pFC3F_GST_061218.gb`
    ## Annotation comes from `overexpression.gtf` file
    ### Locate input files
    labPath = labPathSet[0]
    inputFilename = inputFilenameSet[0]
    inputFileFormat = biopyParamDict[inputFormatSet[0]]
    inputFileBase = stripFilenameSuffix(inputFilename)
    ### Read in genbank file as SeqIO record(s)
    record = SeqIO.read(
                        refPath + labPath + inputFilename, inputFileFormat)
    ### Define name of sequence, to make FASTA header readable  
    ### (parsing original file fields would not be helpful)
    record.id = 'pFC3F_GST'
    record.description = '(' + record.description + ')'
    ### Write FASTA file, for later merging into custom genome sequences
    SeqIO.write(record, refPath + inputFileBase + '.fasta', 'fasta')
    outputFastaRecords = record
     
    # NOCT del plasmid
    labPath = labPathSet[1]
    inputFilename = inputFilenameSet[1]
    inputFileFormat = biopyParamDict[inputFormatSet[1]]
    inputFileBase = stripFilenameSuffix(inputFilename)
    ## read in genbank file as SeqIO record(s)
    record = SeqIO.read(
                        refPath + labPath + inputFilename, inputFileFormat)
    ## define name of sequence (parsing original file fields not helpful)
    record.id = 'pFC3F_NOCT_del2-15'
    record.description = '(' + record.description + ')'
    ## write FASTA file, for later merging into custom genome sequences
    SeqIO.write(record, refPath + inputFileBase + '.fasta', 'fasta')
    # outputFastaRecords = outputFastaRecords + record 
    outputFastaRecords = [outputFastaRecords, record]
    
    # ERCC
    ## Input is FASTA
    ## Annotation needs to be generated
    ### Locate files
    labPath = labPathSet[2]
    inputFilename = inputFilenameSet[2]
    inputFileFormat = biopyParamDict[inputFormatSet[2]]
    ## Read in FASTA file as SeqIO record(s)
    records = SeqIO.parse(
                        refPath + labPath + inputFilename, inputFileFormat)
    gtfBodyERCC = list()
    ### Parse and reformat each fasta entries
    for record in records: 
        gtf = convertFastaSeqRecordToGTF(record, 'Freddolino_Lab')
        gtfBodyERCC.extend(gtf)
        # combine FASTA entries with plasmid FASTA
        outputFastaRecords.append(record)
    
    # Combine all GFF componets from all files 
    # into plain string, one for header segments, 
    # and one for all GFF main body lines
    labPath = labPathSet[3]
    inputFilename = inputFilenameSet[3]
    overexpGTF_file = open(refPath + labPath + inputFilename)
    overexpGTF = overexpGTF_file.readlines()
    overexpGTF_file.close()
    allBody = ''.join(overexpGTF) + '\n'.join(gtfBodyERCC)
    
    # reference genome GTF
    ## read in whole GTF
    outputFileBase = 'customGenomePlasmid_061218'
    lastSegment = '#!genebuild-last-updated 2018-01\n'
    outputFileGTF = open(refPath + outputFileBase + '.gtf', 'w')
    with gzip.open(gtfgzGRCh38Filename, 'rt') as refFileGTF:
        refGTF = refFileGTF.read()
    outputFileGTF.write(refGTF + allBody + '\n')
    outputFileGTF.close()
    
    # Supplementary sequences to reference genome fot customization
    ## Write a supplementary FASTA file conatining sequences for 
    ## GST plasmid, NOCT plasmid, and ERCC sequences
    handle = open(refPath + outputFileBase + '.suppl.fasta', 'w')
    writer = SeqIO.FastaIO.FastaWriter(handle, wrap = 50)
    writer.write_file(outputFastaRecords)
    handle.close()
    
