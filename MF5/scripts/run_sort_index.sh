#!/bin/bash
#Priscila Darakjian
#Script for Alignment with STAR 

for file in *.fastq.gz
do
  file_path="./"$file
  echo $file_path
  /home/exacloud/lustre1/PARC/darakjia/programs/samtools/bin/samtools sort ${file}_Aligned.out.bam -o ${file}_Aligned_sorted.bam
  /home/exacloud/lustre1/PARC/darakjia/programs/samtools/bin/samtools index ${file}_Aligned_sorted.bam
  file_path=""
 done
