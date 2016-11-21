#NOTE: I HARDCODED THE NAMES OF THE SAMPLES IN run_bedtools_coverage.sh
#      to make sure we know the order of the columns in the read coverage result file
#      since multiBamCov does not add a header to the file.

/lawrencedata/share/bedtools2/bin/bedtools multicov -S -split -bams /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10229_S9_L003_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10230_S11_L004_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10231_S23_L008_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10232_S20_L007_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10233_S3_L001_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10235_S15_L005_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10236_S2_L001_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10237_S5_L002_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10238_S21_L007_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/Alignment_MF5/RNA160728RH_Cy_10239_S17_L006_R1_001.fastq.gz_Aligned_sorted.bam -bed /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/MF5_reduced_non_overlapping_exons_sorted.bed > /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MF5/data/RNA160728RH_Cy_coverage_splitoption_stranded_oppositeDirection.txt