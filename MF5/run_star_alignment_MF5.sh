#!/bin/bash
#Priscila Darakjian
#Script for Alignment with STAR 

for file in *.fastq.gz
do
  file_path="./"$file
  echo $file_path
  /home/exacloud/lustre1/PARC/darakjia/programs/STAR_2.3.0e/STAR --genomeDir /home/exacloud/lustre1/PARC/darakjia/genomes/mf5/ --readFilesIn $file_path --readFilesCommand zcat --outFileNamePrefix ${file}_ --runThreadN 12 --outFilterMismatchNmax 3 --outFilterMultimapNmax 1             
  /home/exacloud/lustre1/PARC/darakjia/programs/samtools/bin/samtools view -bS ${file}_Aligned.out.sam > ${file}_Aligned.out.bam
  /home/exacloud/lustre1/PARC/darakjia/programs/samtools/bin/samtools sort ${file}_Aligned.out.bam ${file}_Aligned.out
  /home/exacloud/lustre1/PARC/darakjia/programs/samtools/bin/samtools index ${file}_Aligned.out.bam
  file_path=""
 done
