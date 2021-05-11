#!/bin/bash
#  REQUIRED EXTERNAL TOOLS:
#    trimmomatic 0.39
#    bwa 0.7.17
#    samtools 1.6
#    picard 2.22.1

# SETTINGS
filename="samples.txt"

while read line; do

# performing variant calling (paired-end)
echo $line

basedir="/home/dramazzotti/VirMutSig/step1_data_preprocessing/RNAseq/"
genomeFa="/home/dramazzotti/VirMutSig/step1_data_preprocessing/reference/SARS-CoV-2-ANC.fasta"
resultsdir="/home/dramazzotti/VirMutSig/step1_data_preprocessing/RNAseq/"
fastqTrimDir=${resultsdir}trimmed/${line}/
bwaSamDir=${resultsdir}bwasam/${line}/
bamDir=${resultsdir}bamDir/${line}/
coverageDir=${resultsdir}coverage/${line}/
jobs=16
sampleName=$line

#===B1===

echo "trimmomatic -- trimming fastq"

mkdir -p $fastqTrimDir
if [ ! -d "$fastqTrimDir" ]; then
    echo "Error mkdir"
    exit 1
fi

trimmomatic SE -threads $jobs -phred33 -summary ${fastqTrimDir}${sampleName}_trim.summary -quiet ${basedir}${sampleName}.fastq.gz ${fastqTrimDir}${sampleName}.trim.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

#===E1===

#===B2===

echo "bwa mem -- mapping reads to a reference SARS-CoV-2"

mkdir -p $bwaSamDir
if [ ! -d "$bwaSamDir" ]; then
    echo "Error mkdir"
    exit 1
fi

bwa mem -t $jobs $genomeFa ${fastqTrimDir}${sampleName}.trim.fastq.gz > ${bwaSamDir}${sampleName}_aln.sam

#===E2===

#===B3===

echo "samtools -- building sorted bam"

mkdir -p $bamDir
if [ ! -d "$bamDir" ]; then
    echo "Error mkdir"
    exit 1
fi

samtools view -bT $genomeFa ${bwaSamDir}${sampleName}_aln.sam > ${bamDir}${sampleName}_aln.bam
samtools sort ${bamDir}${sampleName}_aln.bam > ${bamDir}${sampleName}_aln.sorted.bam

#===E3===

#===B4===

echo "picard -- removing duplicates"

picard MarkDuplicates I=${bamDir}${sampleName}_aln.sorted.bam O=${bamDir}${sampleName}_aln.sorted_no_duplicates.bam M=${bamDir}${sampleName}_aln.sorted_no_duplicates_metrics.txt REMOVE_DUPLICATES=true
samtools index -b ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam

#===E4===

#===B5===

mkdir -p $coverageDir
if [ ! -d "$coverageDir" ]; then
    echo "Error mkdir"
    exit 1
fi

echo "samtools depth -- extracting coverage information"

samtools depth -a ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam > ${coverageDir}${sampleName}_aln.sorted_no_duplicates.txt

#===E5===

done < $filename
