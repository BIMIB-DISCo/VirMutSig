#!/bin/bash
#  REQUIRED EXTERNAL TOOLS:
#    samtools 1.6
#    varscan 2.4.4

# SETTINGS
filename="samples.txt"

while read line; do

# performing variant calling (paired-end)
echo $line

genomeFa="/home/dramazzotti/VirMutSig/step1_data_preprocessing/reference/SARS-CoV-2-ANC.fasta"
resultsdir="/home/dramazzotti/VirMutSig/step1_data_preprocessing/AMPLICON/"
bamDir=${resultsdir}bamDir/${line}/
vcfDir=${resultsdir}vcfDir/${line}/
sampleName=$line

#===B1===

echo "samtools mpileup -- building mpileup"

samtools mpileup -f $genomeFa ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam --output ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup

#===E1===

#===B2===

echo "varscan -- calling SNPs"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
    echo "Error mkdir"
    exit 1
fi

varscan pileup2snp ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup --min-var-freq 0.01 --p-value 1 > ${vcfDir}${sampleName}_aln.sorted_no_duplicates.vcf

#===E2===

done < $filename
