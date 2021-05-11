#!/bin/bash
filename="samples.txt"
while read line; do
# creating fastq files (single-end) using SRA Toolkit
echo $line
fastq-dump ${line}
gzip ${line}.fastq
done < $filename
