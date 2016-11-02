#NOTE: I HARDCODED THE NAMES OF THE SAMPLES IN run_bedtools_coverage.sh
#      to make sure we know the order of the columns in the read coverage result file
#      since multiBamCov does not add a header to the file.

/lawrencedata/share/bedtools2/bin/bedtools multicov -S -split -bams /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10208_S1_L001_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10209_S8_L003_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10210_S14_L005_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10210_S14_L005_R1_001.fastq.gz_Aligned.sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10211_S18_L006_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10212_S19_L007_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10213_S4_L002_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10214_S12_L004_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10215_S7_L003_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10216_S10_L004_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10217_S13_L005_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10218_S16_L006_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10220_S24_L008_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10221_S6_L002_R1_001.fastq.gz_Aligned_sorted.bam /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/Alignment_MacaM/RNA160728RH_Rh_10222_S22_L008_R1_001.fastq.gz_Aligned_sorted.bam -bed /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/MacaM_reduced_non_overlapping_exons_sorted.bed > /lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/MacaM/data/RNA160728RH_Rh_coverage_splitoption_stranded_oppositeDirection.txt