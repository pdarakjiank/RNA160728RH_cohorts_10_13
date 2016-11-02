library(GenomicFeatures)
library(Rsamtools)

# Read gtf files
setwd("/lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/")
gtf<-read.delim("/lawrencedata/sharedRNASeq/Rhesus_reference_genome/MacaM/MacaM/UNMC_rhesus_annotation_v7.6.8.gtf", sep="\t", header=FALSE, stringsAsFactors=FALSE)
gtf <- gtf[gtf$V3 == "exon",]

# Extract gene symbol id field (3) from column 9 of gtf
sym.gene.name.gtf<-sapply(strsplit(gtf$V9,";"),function(x) x[3])

# Create separate columns for ensembl ids and gene symbols
# I am not using ensembl ids in this instance but left this as it is for future use, if needed
gtf$gene.sym<-sym.gene.name.gtf

# Generate a GRanges object from gtf and nongtf data frames above for working with chr ranges
gtf.ranges<-GRanges(seqnames=Rle(gtf$V1),ranges=IRanges(start=gtf$V4,end=gtf$V5),strand=gtf$V7,gene.names=gtf$gene.sym)

# Find overlaps between exons 
gtf.range.ovls<-findOverlaps(gtf.ranges,drop.self=T,drop.redundant=T,type="any")

gtf.range.ovls.mat <- as.matrix(gtf.range.ovls)

gtf.range.ovls.dta <- data.frame(as.data.frame(gtf.ranges)[gtf.range.ovls.mat[,1],], as.data.frame(gtf.ranges)[gtf.range.ovls.mat[,2],])

# Verify if there are overlatps between exons of different genes
gtf.ovl.genes.inds <- which(gtf.range.ovls.dta$gene.names != gtf.range.ovls.dta$gene.names.1)

gtf.rm.inds <- as.integer(gtf.range.ovls.mat[gtf.ovl.genes.inds,])

# Prepare a list of the ranges in which exons are overlapping between genes
gtf.ranges.gene.ovl.df <- as.data.frame(gtf.ranges[unique(gtf.rm.inds)])

# Convert data frame to bed format so that we can run bedtools (intersectBed and subtractBed)
gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[c(1,2,3,6,4,5)]
# Remove chromosomes named NT_... and MT
#gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[- grep("MT", gtf.ranges.gene.ovl.df$seqnames),]
#gtf.ranges.gene.ovl.df<-gtf.ranges.gene.ovl.df[- grep("CHR_", gtf.ranges.gene.ovl.df$seqnames),]

# Remove " gene_symbol " string from gene.names column
gtf.ranges.gene.ovl.df$gene.names<-sub(" gene_symbol ","",gtf.ranges.gene.ovl.df$gene.names)


# This is the original list of exons
gtf.ranges.df.unique<-unique(as.data.frame(gtf.ranges))
gtf.ranges.df.unique<-gtf.ranges.df.unique[c(1,2,3,6,4,5)]
gtf.ranges.df.unique$gene.names<-sub(" gene_symbol ","",gtf.ranges.df.unique$gene.names)

# Generate bed files for exons overlapping between genes and for the original list of exons 
write.table(gtf.ranges.gene.ovl.df,"MacaM/data/overlapping_exons_btwn_genes.bed",quote=F,col.names=F,row.names=F,sep="\t")
write.table(gtf.ranges.df.unique,"MacaM/data/exons.bed",quote=F,col.names=F,row.names=F,sep="\t")

# This was run to remove from original list just the portions of exons that overlap between genes
#sort -k1,1 -k2,2n overlapping_exons_btwn_genes.bed > overlapping_exons_btwn_genes_sorted.bed
#sort -k1,1 -k2,2n exons.bed > exons_sorted.bed
#/lawrencedata/darakjia/bedtools-2.17.0/bin/intersectBed -a exons_sorted.bed -b overlapping_exons_btwn_genes_sorted.bed > intersect.bed
#/lawrencedata/darakjia/bedtools-2.17.0/bin/subtractBed -s -a exons_sorted.bed -b intersect.bed > non_overlapping_exons.bed

# Get back the file with no gene overlaps
no_gene_ovl_exons <- read.table("MacaM/data/non_overlapping_exons.bed")
# Convert it to genomic ranges
no.gene.ovl.exons.ranges<-GRanges(seqnames=Rle(no_gene_ovl_exons$V1),ranges=IRanges(start=no_gene_ovl_exons$V2,end=no_gene_ovl_exons$V3),strand=no_gene_ovl_exons$V7,gene.names=no_gene_ovl_exons$V5)

# Merge exon isoforms
no.gene.ovl.exons.ranges.red<-reduce(split(no.gene.ovl.exons.ranges,elementMetadata(no.gene.ovl.exons.ranges)$gene.names))

# Convert above objects into data frame
no.gene.ovl.exons.ranges.red.df<-as.data.frame(no.gene.ovl.exons.ranges.red)

# Re-arrange data frame in bed format and write object to file
no.gene.ovl.exons.ranges.red.df<-no.gene.ovl.exons.ranges.red.df[c(3,4,5,2,6,7)]
no.gene.ovl.exons.ranges.red.df$width<-0

write.table(no.gene.ovl.exons.ranges.red.df,"MacaM/data/MacaM_reduced_non_overlapping_exons.txt", sep="\t",quote=F,col.names=F,row.names=F)

#Run coverage with bedtools
#sort -k6 -k1,1 -k2,2n MacaM_reduced_non_overlapping_exons.txt > MacaM_reduced_non_overlapping_exons_sorted.bed

#nohup ./run_bedtools_coverage.sh > run_bedtools_coverage.log 2>&1 </dev/null &

