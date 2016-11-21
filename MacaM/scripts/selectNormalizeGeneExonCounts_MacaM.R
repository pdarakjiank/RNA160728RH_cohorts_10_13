library(matrixStats)      
library(edgeR)
library(WGCNA)
library(biomaRt)
library(plotrix)

library(foreach)
library(doMC)
registerDoMC()
library(proxy)
library(sgof)
library(multtest)
library(plyr)

getDoParWorkers()
options(cores=6)
getDoParWorkers()

# linux does not work with box so this code has a local data directory

setwd("/lawrencedata/ongoing_analyses/RNA160728RH/160908_D00735_0143_ACA7N4ANXX/RNA160728RH/")

source("MacaM/scripts/functionDefinitions.R")
try(dir.create("MacaM/analysis/resultsCoexpr_MacaM"), silent = F)
try( dir.create("MacaM/analysis/figuresCoexpr_MacaM"), silent = F)
try(dir.create("MacaM/analysis/resultsCoSplicEx_MacaM"), silent = F)
try( dir.create("MacaM/analysis/figuresCoSplicEx_MacaM"), silent = F)

# read raw data - can be improved by using read.table from original .txt file
exonReadsRaw=read.table("MacaM/data/RNA160728RH_Rh_coverage_splitoption_stranded_oppositeDirection.txt",header=F,sep="\t",stringsAsFactors = F)

# Change column 4 name to "gene_sym"
names(exonReadsRaw)[4]<-"gene_sym"

# Combine the chromosome number, start location, and gene symbol to create a unique id  
# column for each exon
exonReadsRaw$exon<-paste(exonReadsRaw$V1,exonReadsRaw$V2,exonReadsRaw$V3,exonReadsRaw$gene_sym,sep="_")

# Sample 10210 was listed twice (everything in both columns are
# exactly the same), so let's fix that by removing one of those columns.
exonReadsRaw<-exonReadsRaw[,-9]

# Create a data frame with gene symbol and exon read counts
exon_counts<-exonReadsRaw[,7:21]
exon_counts<-cbind(exonReadsRaw$gene_sym,exon_counts)

# Calculate the total counts for each gene for each sample
gene_counts<-ddply(exon_counts, 1, numcolwise(sum))

# Change the row names of the exon data frame to the exon unique ids (created above)
rownames(exon_counts)<-exon_counts$exon

# Remove the gene symbol and exon id columns from the exon data frame
exon_counts<-exon_counts[,2:15]

# Change the row names of the gene data frame to the gene symbols
names(gene_counts)[1]<-"gene_sym"
rownames(gene_counts)<-gene_counts$gene_sym

# Remove the gene symbol column from the gene data frame
gene_counts<-gene_counts[,2:15]

# Load sample names in the order they are in the coverage file
samples <- read.table("MacaM/scripts/bam_files.txt",header=F,stringsAsFactors = F)

# Again, as done above, need to remove sample 10210 which is listed twice
samples <- samples[,-3]

names(gene_counts)<-samples[1,]
names(exon_counts)<-samples[1,]

write.table(gene_counts,"MacaM/analysis/RNA160728RH_MacaM_gene_reads_not_normalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
write.table(exon_counts,"MacaM/analysis/RNA160728RH_MacaM_exon_reads_not_normalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
save.image("MacaM/data/Read_Counts_and_Normalization_MacaM.RData")

# IF YOU JUST NEED THE NORMALIZED COUNTS:
# calculate edgeR normalization factors and normalize the data - use all data not just selected
# these normalized data are not used for the DE analysis since edgeR's differential expression
# algorithm (used further down in this script) normalizes the data when calculating the DE
UQnormFactors_exons=calcNormFactors(exon_counts, method=c("upperquartile"))
UQnormFactors_genes=calcNormFactors(gene_counts, method=c("upperquartile"))

effectiveLibrarySizes_exons= UQnormFactors_exons*colSums(exon_counts)
effectiveLibrarySizes_genes= UQnormFactors_genes*colSums(gene_counts)

meanEffLibSize_exons=mean(effectiveLibrarySizes_exons)
meanEffLibSize_genes=mean(effectiveLibrarySizes_genes)

countNormFactor_exons= meanEffLibSize_exons/effectiveLibrarySizes_exons
countNormFactor_genes= meanEffLibSize_genes/effectiveLibrarySizes_genes

normalizedGeneCountsUQ_exons=0* exon_counts
normalizedGeneCountsUQ_genes=0* gene_counts

for (sample in 1:dim(normalizedGeneCountsUQ_exons)[2]){  
  normalizedGeneCountsUQ_exons[,sample]= exon_counts[, sample]* countNormFactor_exons[sample]	
}

for (sample in 1:dim(normalizedGeneCountsUQ_genes)[2]){  
  normalizedGeneCountsUQ_genes[,sample]= gene_counts[, sample]* countNormFactor_genes[sample]  
}

#normalizedGeneCountsUQ_exons =round(normalizedGeneCountsUQ_exons)
#normalizedGeneCountsUQ_genes =round(normalizedGeneCountsUQ_genes)

write.table(normalizedGeneCountsUQ_exons,"MacaM/analysis/RNA160728RH_MacaM_exon_reads_UQNormalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
write.table(normalizedGeneCountsUQ_genes,"MacaM/analysis/RNA160728RH_MacaM_gene_reads_UQNormalized.txt", sep="\t",quote=F,col.names=T,row.names=T)
save.image("MacaM/data/Read_Counts_and_Normalization_MacaM.RData")


########################################################################################
# IF YOU NEED TO CALCULATE DIFFERENTIAL EXPRESSION:
load("MacaM/data/Read_Counts_and_Normalization_MacaM.RData")

# read sample info
sampleKey = read.table("MacaM/data/rhesus_phenotype.txt",header=F,sep="\t")
sampleKey = as.data.frame(t(sampleKey[,2:15]))

samplesDrinkers<- subset(sampleKey,sampleKey$V2=="Drinker")
samplesControls<- subset(sampleKey,sampleKey$V2!="Drinker")

# this cleans up the subsets of the levels that remained from sampleKey to avoid further headaches
samplesDrinkers<- levels(droplevels(samplesDrinkers$V1))
samplesControls<- levels(droplevels(samplesControls$V1))


countsDrinkers=gene_counts[,substring(samples,16,20) %in% samplesDrinkers]
countsControls=gene_counts[,substring(samples,16,20) %in% samplesControls]

save(samplesControls,samplesDrinkers,countsControls,countsDrinkers,sampleKey, file="DE_MacaM.RData")

########################################################################################
load("MacaM/data/DE_MacaM.RData")
# divide the data in different groups and perform DE with edgeR

groupSelection=c(rep("countsDrinkers",dim(countsDrinkers)[2]),rep("countsControls",dim(countsControls)[2]))
groupSelection =factor(groupSelection)

d=DGEList(counts= cbind(countsDrinkers, countsControls), group= groupSelection)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d, dispersion="tagwise") 

d2 <- d[rowSums(1e+06 * d$counts/expandAsMatrix(d$samples$lib.size, dim(d)) > 1) >= 3, ]
pdf( "MDS_plot_Code_Noncode_normalized_gene_counts.pdf" , width = 7 , height = 7 )
plotMDS( d2 , main = "MDS Plot for Count Data - Drinkers x Controls", labels = d2$samples$group,cex=0.5 )
dev.off()

#from the old script (I did this to compare to the results from SGoF; no difference)
#####
de.calls <- decideTestsDGE(de.tgw, p=0.05, adjust.method="fdr") # returns 0 de calls; 104 if no adjustment
sum(de.calls)
resultsDEtotal=cbind(de.tgw$table, de.calls)

#####
#####
# use sgof package for multiple comparison correction
# results from sgof come out sorted but un-named (!!!!) so the original pvalues and geneNames need to be sorted

pValues_DE=de.tgw$table$PValue
geneNames<-rownames(de.tgw$table)
names(pValues_DE)=geneNames

sortIndexes=sort.int(pValues_DE, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNames[sortIndexes]

adjustedResults<-SGoF(u=pValues_DE)
summary(adjustedResults)

sortedAdjustedPvals_DE=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals_DE)=sortedGeneNames

fileConnSummary<-file("resultsCoexpr_MacaM/SummaryResultsCoexpr.txt",  open="wt")

writeLines(paste("Number of genes with FDR < 0.05: ", sum(sortedAdjustedPvals_DE < 0.05), sep=' '), fileConnSummary)
close(fileConnSummary)

# some sanity checks
plot(pValues_DE[sortedGeneNames], adjustedResults$data) # should be straight line

##############################################################################
# verify that smaller libraries are multiplied by bigger normalization factors
plot(colSums(gene_counts), countNormFactor_genes, xlab="Un-normalized library sizes", ylab="Normalization factors")

# plot original library sizes versus normalized
#range of data much smaller in normalized data
xylim=c(min(colSums(gene_counts), colSums(normalizedGeneCountsUQ_genes)), max(colSums(gene_counts), colSums(normalizedGeneCountsUQ_genes)))
plot(colSums(gene_counts), colSums(normalizedGeneCountsUQ_genes), xlim=xylim, ylim=xylim)

## standard deviation  of normalized libraries should me much smaller
sd(colSums(gene_counts))
sd(colSums(normalizedGeneCountsUQ_genes))


###################################################################################
# select genes with logCPM > 0 (equivalent to CPM>1) 

geneNamesHighCPM=geneNames[de.tgw$table[,"logCPM"]>0]
geneCountsHighCPM=normalizedGeneCountsUQ_genes[geneNamesHighCPM,]

meanCountsDrinkers=rowMeans(countsDrinkers[geneNamesHighCPM,])
meanCountsControls=rowMeans(countsControls[geneNamesHighCPM,])

sdCountsDrinkers=apply(countsDrinkers[geneNamesHighCPM,],1, sd) 
sdCountsControls=apply(countsControls[geneNamesHighCPM,],1, sd) 


##############################################################################
# find differentially variable genes

pvalVar=rep(1, length(geneNamesHighCPM))
names(pvalVar)=geneNamesHighCPM

for (gene in geneNamesHighCPM) 
{
  pvalVar[gene]=var.test(x=t(countsDrinkers[gene,]), y=t(countsControls[gene,]))$p.value
}

pvalVar[is.na(pvalVar)]=1

pValues=pvalVar
names(pValues)=geneNamesHighCPM

sortIndexes=sort.int(pValues, decreasing = F, index.return=T)$ix
sortedGeneNames=geneNamesHighCPM[sortIndexes]

adjustedResults<-SGoF(u=pValues)
summary(adjustedResults)

sortedAdjustedPvals_DV=adjustedResults$Adjusted.pvalues
names(sortedAdjustedPvals_DV)=sortedGeneNames

fileConnSummary<-file("resultsCoexpr_MacaM/SummaryResultsCoexpr.txt",  open="at")

writeLines(paste("Number of genes with DV < FDR=0.05: ", sum(sortedAdjustedPvals_DV<0.05), sep=' '), fileConnSummary)
close(fileConnSummary)

geneNamesDE=sortedGeneNames[sortedAdjustedPvals_DE < 0.05]
geneNamesDV=sortedGeneNames[sortedAdjustedPvals_DV < 0.05]

write.csv(geneNamesDE, file="resultsCoexpr_MacaM/geneNamesDE.csv")
write.csv(geneNamesDV, file="resultsCoexpr_MacaM/geneNamesDV.csv")

setdiff(geneNamesHighCPM, geneNames)

###########################################################################33
results_highCPMgenes=cbind(de.tgw$table[geneNamesHighCPM,], meanCountsControls[geneNamesHighCPM], meanCountsDrinkers[geneNamesHighCPM], pValues_DE[geneNamesHighCPM], sortedAdjustedPvals_DE[geneNamesHighCPM],sdCountsControls[geneNamesHighCPM],sdCountsDrinkers[geneNamesHighCPM],  pvalVar, sortedAdjustedPvals_DV[geneNamesHighCPM] )
results_highCPMgenes=round(results_highCPMgenes,3)
colnames(results_highCPMgenes)=c(colnames(de.tgw$table), c("mean counts controls", "mean counts drinkers", " p val DE", " adj p DE", "sd controls", "sd drinkers", "p val DV", "adj p val DV")) 
rownames(results_highCPMgenes)=geneNamesHighCPM
#this will be collected in Supplemental Table 1 
write.csv(results_highCPMgenes, file="resultsCoexpr_MacaM/resultsDEDV_highCPM_MacaM.csv")


#######################################################################################3
# possibly swithch to bicor correlation in the future
#adjCoexpr=adjacency(t(geneCountsHighCPM), corFnc = "bicor", type="unsigned", power=6)

adjCoexpr=adjacency(t(geneCountsHighCPM),  type="unsigned", power=6)
adjCoexpr[is.na(adjCoexpr)]=0

#sanity check, should be 0
sum(adjCoexpr<0, na.rm=T)
connCoexpr=rowSums(adjCoexpr, na.rm = T)-1
connCoexpr_WGCNA=softConnectivity(t(geneCountsHighCPM), type="unsigned", power=6)
plot(connCoexpr, connCoexpr_WGCNA)

names(connCoexpr)=geneNamesHighCPM

connCoexpr=connCoexpr[connCoexpr > 0]

connCoexpr[is.na(connCoexpr)]=0
sortedConn=sort(connCoexpr, decreasing = T)
sum(connCoexpr < 0, na.rm=T)

totalConn=sum(sortedConn)
cumulativeConnFraction=0*sortedConn
for (i in 1:length(sortedConn)){
  cumulativeConnFraction[i]=sum(sortedConn[1:i])/totalConn
}


lastGene=min(which(cumulativeConnFraction > .9))

highConnGenes=names(sortedConn)[1:lastGene]

plot(cumulativeConnFraction, xlab="", ylab="")
title(xlab="Number of genes \n ranked by connectivity", cex.lab=1.25, font.lab=2)
title(ylab="Fraction total connectivity", cex.lab=1.25, font.lab=2)
abline(v=lastGene)
abline(h=0.9)
text(x=8100, y=0.8, labels=paste(lastGene, " genes", sep=""))
text(x=3100, y=0.95, labels="90% connectivity captured")
title(main="Selecting the network size\n starting from genes with > 1 CPM", cex.lab=1.5, font.lab=2)
# export figure
#############################################################################################################
adjCoexprHighConn=adjCoexpr[highConnGenes,highConnGenes]

selectedGeneCounts=normalizedGeneCountsUQ[highConnGenes,]
###################################################################################
exonCounts=read.table("RNASeq019 HSCC Alex/SH/RNASeq019_SH_mm10_exon_reads_not_normalized.txt")


splitIDs=mapply(strsplit, rownames(exonCounts), MoreArgs=list(split="_", fixed = FALSE, perl = FALSE, useBytes = FALSE))
exonGeneName=unlist(lapply(splitIDs, "[[", 4))
exon_start=unlist(lapply(splitIDs, "[[", 2))

#sanity check

colnames(exonCounts)==colnames(geneReadsRaw)

# sanity check
sum(exonCounts[exonGeneName=="Drd2",])
sum(selectedGeneCounts["Drd2",])

normExonCounts=0* exonCounts
for (sample in 1:dim(exonCounts)[2]){
  normExonCounts[,sample]= exonCounts[, sample]* countNormFactor[sample]	
}
normExonCounts =round(normExonCounts)
rownames(normExonCounts)=rownames(exonCounts)

# select exons from genes with at least 1 CPM
exonCountsHighCounts=normExonCounts[which(exonGeneName %in% geneNamesHighCPM),]
exonGeneNamesHighCounts=exonGeneName[exonGeneName %in% geneNamesHighCPM]

# compute Canberra distances 
canberraListExons=foreach (geneName = geneNamesHighCPM, .inorder=T, .verbose = T) %dopar% {
  currExonCounts= exonCountsHighCounts[which(exonGeneNamesHighCounts==geneName),]	
  if (is.null(dim(currExonCounts))){
    exonDistMatrix=as.matrix(dist(as.matrix(currExonCounts), method="canberra"))
  } else {
    exonDistMatrix=as.matrix(dist(t(as.matrix(currExonCounts)), method="canberra"))
  }
  # colnames(exonDistMatrix)=exonColnames
  # rownames(exonDistMatrix)=exonColnames
  exonDistMatrix
  
}
names(canberraListExons)=geneNamesHighCPM
save(canberraListExons, file="resultsCoSplicEx_SH/canberraListExonsSH.RData")
 load("resultsCoSplicEx_SH/canberraListExonsSH.RData")
########################################################################################################

nGenes=length(canberraListExons)
gene_indexes=1:nGenes

# reformat the data so one can use WGCNA adjacency function to construct CoSplicEx adjacency matrix
lengthVector=length(as.vector(as.dist(canberraListExons[[1]])))
distData=matrix(data=0, nrow=lengthVector, ncol=length(canberraListExons))
colnames(distData)=names(canberraListExons)

for(gene in names(canberraListExons)) {
  distData[,gene]=as.vector(as.dist(canberraListExons[[gene]]))
} 

##########################################################################################################
# compute connectivity in large CoSplicEx network

#######################################################################################3
adjCoSplicEx_raw=adjacency(distData,  type="unsigned", power=6)

adjCoSplicEx_raw[is.na(adjCoSplicEx_raw)]=0

#sanity check, should be 0
sum(adjCoSplicEx_raw<0, na.rm=T)
connCoSplicEx=rowSums(adjCoSplicEx_raw, na.rm = T)-1

names(connCoSplicEx)=geneNamesHighCPM

connCoSplicEx=connCoSplicEx[connCoexpr > 0]

connCoSplicEx[is.na(connCoSplicEx)]=0
sortedConn=sort(connCoSplicEx, decreasing = T)
sum(connCoSplicEx < 0, na.rm=T)

totalConn=sum(sortedConn)
cumulativeConnFraction=0*sortedConn
for (i in 1:length(sortedConn)){
  cumulativeConnFraction[i]=sum(sortedConn[1:i])/totalConn
}


lastGene=min(which(cumulativeConnFraction > .9))

highConnGenes=names(sortedConn)[1:lastGene]

plot(cumulativeConnFraction, xlab="", ylab="")
title(xlab="Number of genes \n ranked by connectivity", cex.lab=1.25, font.lab=2)
title(ylab="Fraction total connectivity", cex.lab=1.25, font.lab=2)
abline(v=lastGene)
abline(h=0.9)
text(x=6100, y=0.8, labels=paste(lastGene, " genes", sep=""))
text(x=3100, y=0.95, labels="90% connectivity captured")
title(main="Selecting the network size\n starting from genes with > 1 CPM", cex.lab=1.5, font.lab=2)
# export figure
#############################################################################################################
adjCoSplicEx=adjCoSplicEx_raw[highConnGenes,highConnGenes]


canberraListSelected=canberraListExons[highConnGenes]

exonGeneNameSelected=highConnGenes
selectedExonCounts=normExonCounts[which(exonGeneName %in% exonGeneNameSelected),]

########################################################################################################

save(selectedGeneCounts, canberraListSelected,adjCoSplicEx,selectedExonCounts, exonGeneNameSelected, groupSelection, samplesHigh, samplesLow, file="selectedData_SH.RData")
#load("data/selectedData.RData")



############################################################################################
