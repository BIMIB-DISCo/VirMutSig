#!/bin/bash
filename="samples.txt"
while read line; do
# creating fastq files (paired-end) using SRA Toolkit
echo $line
fastq-dump --split-files ${line}
gzip ${line}_1.fastq
gzip ${line}_2.fastq
done < $filename
