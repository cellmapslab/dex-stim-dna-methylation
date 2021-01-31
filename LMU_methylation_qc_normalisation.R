library(minfi)
library(minfiData)
library(minfiDataEPIC)
library(RColorBrewer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(sva)
library(RPMM)
library(ggplot2)
library(missMethyl)
library(matrixStats)
library(wateRmelon)
library(reshape)
library(Hmisc)
library(xlsx)
library(impute)

setwd("/binder/jade/KSP-LMU/methylation/")

source("BMIQ_1.6_Teschendorff.R")

sessionInfo <- sessionInfo()
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#read in data
targets = read.table("idat_locations.txt", header=T)
RGset_LMU = read.metharray.exp(targets = targets, verbose=T)
save(RGset_LMU, file="RGset_LMU.rda")
targets_rerun = read.table("rerun/idat_locations_rerun.txt")
RGset_rerun = read.metharray.exp(targets = targets_rerun, verbose=T)
save(RGset, file="rerun/RGset_rerun.rda")
RGset_rerun = RGset_rerun[,c(1:8)] #remove tube samples, use ony repl
# RGset_rerun = RGset_rerun[,c(9:16)] for tube samples

#merge RGsets
RGset_LMU_all = combineArrays(RGset_LMU, RGset_rerun, verbose=T)
save(RGset_LMU_all, file="RGset_LMU_all.rda")
#dim(RGset_LMU_all): 1051539 648
#mv RGset_LMU_all.rda RGset_LMU_all_REPL.rda
#save(RGset_LMU_all, file="RData/RGset_LMU_all_TUBE.rda")

#get data
samplesheet1 = read.xlsx("/binder/common/methylation/raw_methylation/Projekt_M01039_LMU/M01039_methEPIC/M01039_01-07_samplesheet.xlsx",sheetIndex=1)
samplesheet2 = read.csv("/binder/common/methylation/raw_methylation/Projekt_M01039_LMU/M01039_WDH_EPIC_KSP_LMU_8samples/M01039_WDH_samplesheet.csv")
pheno = read.table("samplesheet_epic_kjp_final.txt", header=T)

header = samplesheet1[c(8),]
samplesheet1 = samplesheet1[-c(1:8),]
samplesheet2 = samplesheet2[c(11:18),] #samplesheet2 = samplesheet2[c(19:26),] for TUBE samples
colnames(samplesheet1) = header
colnames(samplesheet2) = header
samplesheet_all = rbind(samplesheet1, samplesheet2)
samplesheet_all = samplesheet_all[! samplesheet_all$Sentrix_ID == "203244490189",] #remove failed chip from first run
samplesheet_all$Sample_Name = gsub("_REP","",samplesheet_all$Sample_Name) #rename replicated samples for merge with phenotype file
# samplesheet_all$Sample_Name = gsub("_TUBE","",samplesheet_all$Sample_Name)
write.csv(samplesheet_all, file="samplesheet_all.csv",row.names=F)
# write.csv(samplesheet_all, file="samplesheet_all_TUBE.csv",row.names=F)
#mv samplesheet_all.csv samplesheet_all_REPL.csv

pheno_all = merge(pheno,samplesheet_all, by.x="Proben_name",by.y="Sample_Name")
pheno_all$MethID = paste(pheno_all$Sentrix_ID, pheno_all$Sentrix_Position, sep="_")
write.csv(pheno_all, file="pheno_and_samplesheet_final.csv")
# write.csv(pheno_all, file="pheno_and_samplesheet_final_TUBE.csv")
#dim(pheno_all): 640  18
#mv pheno_and_samplesheet_final.csv pheno_and_samplesheet_final_repl.csv

RGset_LMU_all_nodupl = RGset_LMU_all[,colnames(RGset_LMU_all) %in% pheno_all$MethID] #remove failed chip from RGset of first run
save(RGset_LMU_all_nodupl, file="RData/RGset_LMU_all_nodupl.rda")
#dim(RGset_LMU_all_nodupl): 1051539     640
#mv RGset_LMU_all_nodupl.rda RGset_LMU_all_nodupl_REPL.rda
#save(RGset_LMU_all_nodupl, file="RData/RGset_LMU_all_nodupl_TUBE.rda")

#save original data sets
Mset = preprocessRaw(RGset_LMU_all_nodupl) # cpG locus level, with 2 channels methylated/ unmethylated
save(Mset,file="RData/Mset_LMU_original.Rdata")
#save(Mset,file="RData/Mset_LMU_original_TUBE.Rdata")
#mv Mset_LMU_original.Rdata Mset_LMU_original_REPL.Rdata

RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)# CpG locus level, but not mapped to a genome, Beta and M values
save(RatioSet, file="RData/RatioSet_LMU_original.Rdata")
#save(RatioSet, file="RData/RatioSet_LMU_original_TUBE.Rdata")
#mv RatioSet_LMU_original.Rdata RatioSet_LMU_original_REPL.Rdata

RawBetas = getBeta(RatioSet) # Beta value matrix
save(RawBetas,file="RData/RawBetas_LMU_original.Rdata")
#mv RawBetas_LMU_original.Rdata RawBetas_LMU_original_REPL.Rdata

gRatioSet<-mapToGenome(RatioSet,mergeManifest=TRUE)
save(gRatioSet, file = "RData/gRatioSet_LMU_original.Rdata")
#mv gRatioSet_LMU_original.Rdata gRatioSet_LMU_original_REPL.Rdata

#check quality
#detection P-value
detP <- detectionP(RGset_LMU_all_nodupl)
head(detP)
save(detP, file = "RData/detP.Rdata")
#save(detP, file="RData/detP_TUBE.Rdata")
#mv detP.Rdata detP_REPL.Rdata

#plot mean detection P values
targets = read.csv("pheno_and_samplesheet_final.csv")
targets = targets[order(targets$MethID),]
detP = detP[,order(colnames(detP))]
#pdf("Reports/detectionP_TUBE.pdf", width = 8, height = 3)
pdf("Reports/detectionP.pdf", width = 8, height = 3)
pal <- brewer.pal(8,"Dark2")
barplot(colMeans(detP), col=pal[factor(targets$Proben_name)], las=2, cex.names=0.4,ylab="Mean detection p-values")
abline(h=0.01,col="red")
dev.off()
#mv detectionP.pdf detectionP_REPL.pdf

#filter samples with low det
keep <- colMeans(detP) < 0.05 #keep = 640
RGSet_qual <- RGset_LMU_all_nodupl[,keep]
save(RGSet_qual, file = "RData/RGSet_qual.Rdata")
#mv RGSet_qual.Rdata RGSet_qual_REPL.Rdata
#save(RGset_qual, file="RData/RGSet_qual_TUBE.Rdata")

#beta values for filtered RGset
RawBetas_qual <- RawBetas[,keep]
save(RawBetas_qual, file = "RData/RawBetas_qual.Rdata")
#mv RawBetas_qual.Rdata RawBetas_qual_REPL.Rdata

#filter targets
detP_qual <- detP[,keep]
save(detP_qual, file = "RData/detP_qual.Rdata")
#mv detectionP_qual.Rdata detectionP_qual_REPL.Rdata
#save(detP, file="RData/detP_qual_TUBE.Rdata")

#filter samplesheet
targets = read.csv("pheno_and_samplesheet_final.csv")
targets_qual = targets[targets$MethID %in% keep,]
write.csv(targets_qual, file="pheno_and_samplesheet_qual.csv")
#mv pheno_and_samplesheet_qual.csv pheno_and_samplesheet_qual_REPL.csv

#excluded after low dectectionP: "203249030186_R03C01", "203249030189_R01C01"

#quality control
qcReport(RGset_qual, sampGroups=targets_qual$Slides, pdf="Reports/qcReport.pdf")
#mv qcReport.pdf qcReport_REPL.pdf
#qcReport(RGset_qual, sampGroups=targets_qual$Slides, pdf="Reports/qcReport_TUBE.pdf")

# plot distribution artefacts
# pdf("Reports/beta_densities_TUBE.pdf")
pdf("Reports/beta_densities.pdf")
for (i in 1:ncol(RGSet_qual))
{titel<-paste(rownames(pData(RGSet_qual))[i])
densityPlot(as.matrix(RawBetas_qual[,i]),main=titel)
print(i)}
dev.off()
#mv beta_densities.pdf beta_densities_REPL.pdf

#predict Sex
load("RData/gRatioSet_original.Rdata")
predictedSex <- getSex(gRatioSet,cutoff=-2)
targets_qual = read.csv("pheno_and_samplesheet_qual_REPL.csv")
#targets_qual = read.csv("pheno_and_samplesheet_qual_TUBE.csv")
rownames(targets_qual) = targets_qual$MethID
predictedSex = merge(targets_qual, predictedSex, by="row.names")
predictedSex[predictedSex$Sex == "f",]$Sex = "F"
predictedSex[predictedSex$Sex == "m",]$Sex = "M"
save(predictedSex, file="RData/predictedSex_REPL.rda")
#save(predictedSex, file="RData/predictedSex_TUBE.rda")

#mismatches: 2 samples
predictedSex[! predictedSex$Sex == predictedSex$predictedSex,]
#202212330092_R06C01       MU1B3-1794801
#202232360147_R08C01       MU1B3-1759301

#exclude mismatches and save clean sets
out_index = c("202212330092_R06C01","202232360147_R08C01")
RGSet_clean <- RGSet_qual[ , !(colnames(RGSet_qual) %in% out_index)]
save(RGset_clean, file="RData/RGset_clean.rda")
#save(RGset_clean, file="RData/RGset_clean_TUBE.rda")
#mv RGset_clean.rda RGset_clean_REPL.rda

RawBetas_clean = RawBetas_qual[ ,!(colnames(RawBetas_qual) %in% out_index)]
save(RawBetas_clean, file = "RData/RawBetas_clean.Rdata")
#save(RawBetas_clean, file = "RData/RawBetas_clean_TUBE.Rdata")
#mv RawBetas_clean.Rdata RawBetas_clean_REPL.Rdata

detP_clean <- detP_qual[,!(colnames(detP_qual) %in% out_index)]
save(detP_clean, file = "RData/detP_clean.Rdata")
#save(detP_clean, file = "RData/detP_clean_TUBE.Rdata")
#mv detP_clean.Rdata detP_clean_TUBE.Rdata

## Get annotations probes:
annot = getAnnotation(RGSet_clean)
save(annot,file="RData/annot_clean.Rdata")
#save(annot,file="RData/annot_clean_TUBE.Rdata")
#mv annot_clean.Rdata annot_clean_REPL.Rdata

#normalization (standard procedure): QN + BMIQ / BMIQ
quantileN = preprocessQuantile(RGSet_clean)
save(quantileN, file="RData/quantileN.Rdata")
#save(quantileN, file="RData/quantileN_REPL.Rdata")
#mv quantileN.RData quantileN_TUBE.RData

quantileNBetas = getBeta(quantileN) # get betas
save(quantileNBetas,file="RData/quantileNBetas.Rdata")
#save(quantileNBetas,file="RData/quantileNBetas_TUBE.Rdata")
#mv quantileNBetas.Rdata quantileNBetas_REPL.Rdata

quantileNMs = getM(quantileN)
save(quantileNMs,file="RData/quantileNMs.Rdata")
#save(quantileNMs,file="RData/quantileNMs_TUBE.Rdata")
#mv quantileNMs.Rdata quantileNMs_REPL.Rdata


# BMIQ after quantile normalization:
probeType = as.data.frame(annot[rownames(quantileNBetas),c("Name","Type")])
probeType$probeType = ifelse(probeType$Type %in% "I",1,2)

BMIQ.quantileN = apply(quantileNBetas[,1:length(colnames(quantileNBetas))],2,function(a) BMIQ(a,probeType$probeType,plots=FALSE)$nbeta)

length(which(is.nan(BMIQ.quantileN))) # should be 0
save(BMIQ.quantileN, file="RData/BMIQ.quantileN.Rdata")
#save(BMIQ.quantileN, file="RData/BMIQ.quantileN_TUBE.Rdata")
#mv BMIQ.quantileN.Rdata BMIQ.quantileN_REPL.Rdata

#Check distributions before and after normalization
png(file="Reports/Beta_Distributions_QN_BMIQ_Norm.png")
#png(file="Reports/Beta_Distributions_QN_BMIQ_Norm_TUBE.png")
par(mfrow=c(1,2))
densityPlot(RawBetas_clean, sampGroups = targets_clean$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta")
densityPlot(BMIQ.quantileN, sampGroups = targets_clean$Slide, legend=FALSE, main = "quantileBMIQ adjusted Betas", xlab = "Beta")
dev.off()
#mv Beta_Distributions_QN_BMIQ_Norm.png Beta_Distributions_QN_BMIQ_Norm_REPL.png
#continue with REPL

#remove NAs: 6846
RawBetas_clean_noNA <- RawBetas_clean[rowSums(is.na(RawBetas_clean)) == 0,]

#BMIQ only
probeType = as.data.frame(annot[rownames(RawBetas_clean_noNA),c("Name","Type")])
probeType$probeType = ifelse(probeType$Type %in% "I",1,2)
BMIQ.only = apply(RawBetas_clean_noNA[,1:length(colnames(RawBetas_clean_noNA))],2,function(a) BMIQ(a,probeType$probeType,plots=FALSE)$nbeta)
length(which(is.nan(RawBetas_clean_noNA))) # should be 0
save(BMIQ.only, file="BMIQ_noNA.only.Rdata")

#Compare normalizations
png(file="Reports/Beta_Distributions_comparison_Norm_REPL.png")
par(mfrow=c(1,2))
densityPlot(BMIQ.only, sampGroups = targets_clean$Slide, legend=FALSE, main = "BMIQ adjusted Betas", xlab = "Beta")
densityPlot(BMIQ.QN, sampGroups = targets_clean$Slide, legend=FALSE, main = "quantileBMIQ adjusted Betas", xlab = "Beta")
dev.off()

#contiue with QN + BMIQ

#filter probes
load("RData/detP_clean_REPL.Rdata")
load("RData/BMIQ.quantileN_REPL.Rdata")
load("RData/quantileN_REPL.Rdata)
load("RData/nnot_clean_REPL.Rdata")
# ensure probes are in the same order in the BMIQ beta and detP objects
detP_clean_f <- detP_clean[match(rownames(BMIQ.quantileN),rownames(detP_clean)),]
# ensure probes are in the same order in the quantileN (GenomicRatio) and detP objects
detP_clean_ff <- detP_clean[match(featureNames(quantileN),rownames(detP_clean)),]
keep <- rowSums(detP_clean_f < 0.01) == ncol(BMIQ.quantileN)
#table(keep)
#FALSE   TRUE
#58387 807472
keep_ff <- rowSums(detP_clean_ff < 0.01) == ncol(quantileN)
#table(keep_ff)
#FALSE   TRUE
#58387 807472
BMIQ.quantileN_filtered <- BMIQ.quantileN[keep,] #dim(BMIQ.quantileN_filtered): 807472    636
quantileN_filtered <- quantileN[keep_ff,] #dim: 807472 636

# filter out probes on sex chromosomes
keep <- !(rownames(BMIQ.quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
#table(keep)
# FALSE   TRUE
# 17147 790325
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] #dim(BMIQ.quantileN_filtered): 790325    636
keep_ff <- !(featureNames(quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
#dim(keep_ff)
# FALSE   TRUE
# 17147 790325
quantileN_filtered <- quantileN_filtered[keep_ff,] #dim: 790325 636

# removal of probes where common SNPs may affect the CpG
quantileN_filtered <- dropLociWithSnps(quantileN_filtered) # this NEEEDs a GenomicRatioSet, beta matrix does not work
quantileN_filtered # dim: 766100 636
BMIQ.quantileN_filtered <-  BMIQ.quantileN_filtered[rownames(BMIQ.quantileN_filtered) %in% featureNames(quantileN_filtered),]
#dim(BMIQ.quantileN_filtered): 766100    636

# load Chen Probe annotation file and exclude SNP and Cross Hybridizing probes:
load("addFiles/ChenProbeIDs.rdata")
annot2$SNPs = annot2[,"EUR"]
index<-which(annot2$sex=="Exclude" | annot2$CRSH=="Exclude" | annot2$EUR=="Exclude")
length(index) #62520
exclude_Chen <-annot2[index,]

exclude_crosshyb <-read.table("addFiles/CpGs_crosshybridizing_EPIC.txt",sep="",header=F)
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_crosshyb$V1)
#table(keep)
#FALSE   TRUE
# 36818 729282
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
#dim(BMIQ.quantileN_filtered)  729282    636

keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_crosshyb$V1)
# table(keep_ff)
# FALSE   TRUE
# 36818 729282
quantileN_filtered <- quantileN_filtered[keep_ff,]
quantileN_filtered # dim: 729282 636

exclude_poly <-read.table("addFiles/CpGs_polymorphic_EPIC.txt",sep="",header=T)
index<-which(exclude_poly$EUR_AF>0.05) #n=10971
exclude_polym <- exclude_poly[index,]

keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_polym$IlmnID)
# table(keep)
# FALSE   TRUE
#  414 728868
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] #dim(BMIQ.quantileN_filtered): 728868    636
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_polym$IlmnID)
# table(keep_ff)
#  FALSE   TRUE
#   414 728868
quantileN_filtered <- quantileN_filtered[keep_ff,]

save(BMIQ.quantileN_filtered, file="RData/BMIQ.quantileN_filtered.Rdata") # QN+BMIQ normalized, probe filtered & sample-qced Betas.
save(quantileN_filtered, file="RData/quantileN_filtered.Rdata") # GenomicRatio of filtered QN normalized data
Betas_quantileN_filtered <- getBeta(quantileN_filtered) # beta values
save(Betas_quantileN_filtered, file="RData/Betas_quantileN_filtered.Rdata")
Ms_quantileN_filtered <- getM(quantileN_filtered)
save(Ms_quantileN_filtered, file="RData/Ms_quantileN_filtered.Rdata")

# combat for batch correction
# function for variance explained
rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

#get M-values
mval <- apply(BMIQ.quantileN_filtered, 2, function(x) log2((x)/(1-x)))

## Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars = as.matrix(rowVars(mval))
which(vars==0) # 0

## Replace all probes with no variance with NA and remove them from the normalized data set
vars[vars == 0] = NA
vars = na.omit(vars)
intersect = intersect(rownames(vars), rownames(mval))
print(length(intersect)) # 728868 probes without variance == 0

BMIQ.quantileN_filtered_batch = BMIQ.quantileN_filtered[intersect,]
mval = mval[intersect,]

## Ensure Objects are aligned
pd_clean = pData(quantileN_filtered)
table(ifelse(rownames(pd_clean) == colnames(mval),"Match","Off")) # All should match 636

#run pca
PCobj = prcomp(t(mval), retx = T, center = T, scale. = T)
save(PCobj, file="RData/PCobj.Rdata")
pdf("Reports/boxplot_PCA.pdf")
boxplot(PCobj$x,col="grey",frame=F)
dev.off()

#Scree plot to determine number of PCs to keep
pdf("Reports/screeplot_PCA.pdf")
plot(PCobj,type="line",cex.lab=1.5, cex.main=1.5)
dev.off()

# Extract the PCs from the PCobj object
PCs = PCobj$x
R = 10
propvar = summary(PCobj)$importance["Proportion of Variance", 1:R]
#    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
#0.05344 0.03656 0.03287 0.02492 0.01469 0.01447 0.00818 0.00759 0.00600 0.00509
cummvar = summary(PCobj)$importance["Cumulative Proportion", 1:R]
#    PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8     PC9    PC10
#0.05344 0.09000 0.12287 0.14779 0.16248 0.17694 0.18512 0.19271 0.19871 0.20380

R = 5 # choose eg 5 PCs - elbow
## Generate plots of the resulting PCs
# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs
pdf("Reports/PCA_variance_explained.pdf")
par(mfrow=c(1,2))
par(mar = c(5,5,4,2))
barplot(propvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
par(mar = c(5,5,4,2))
barplot(cummvar*100, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
dev.off()

PCs = PCobj$x
PCs =PCs[,1:R]
pd_clean = read.csv("pheno_and_samplesheet_clean_REPL.csv")
rownames(pd_clean) = pd_clean$MethID
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T)

pdf(file="Reports/PC_Variation_by_batch.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Platte), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sentrix_ID), xlab = "PC1", ylab = "PC2", main="Array")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$pos96well), xlab = "PC1", ylab = "PC2", main="Position")
dev.off()

#check for extreme outliers: none
o1 <- 3*sd(Prin.comp$PC1)
o2 <- 3*sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2)

#anova to check strong variations in PCs by Batch
anova(lm(Prin.comp$PC1~Prin.comp$Platte))
anova(lm(Prin.comp$PC2~Prin.comp$Platte))
anova(lm(Prin.comp$PC3~Prin.comp$Platte))
anova(lm(Prin.comp$PC4~Prin.comp$Platte))
anova(lm(Prin.comp$PC5~Prin.comp$Platte)) #2.209e-06 ***

anova(lm(Prin.comp$PC1~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC2~Prin.comp$Sentrix_ID)) #< 2.2e-16 ***
anova(lm(Prin.comp$PC3~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC4~Prin.comp$Sentrix_ID)) #1.299e-14 ***
anova(lm(Prin.comp$PC5~Prin.comp$Sentrix_ID))

anova(lm(Prin.comp$PC1~Prin.comp$pos96well))
anova(lm(Prin.comp$PC2~Prin.comp$pos96well))
anova(lm(Prin.comp$PC3~Prin.comp$pos96well))
anova(lm(Prin.comp$PC4~Prin.comp$pos96well)) #2.7e-13 ***
anova(lm(Prin.comp$PC5~Prin.comp$pos96well))

#first correct for plate
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
pd_clean = pd_clean[order(rownames(pd_clean)),]
table(ifelse(rownames(pd_clean) == rownames(PCobj$x),"Match","Off")) #match: 636

#median substitution for NA: median(pd_clean$risk_all, na.rm=T)
#pd_clean[is.na(pd_clean$risk_all),]$risk_all = 12
mod <- model.matrix(~Case.Control+risk_all, data=pd_clean)
M_combat_1plate = ComBat(mval,batch = pd_clean$Platte, mod = mod)
save(M_combat_1plate,file="RData/M_combat_1plate.Rdata")

# check if associations are still present
PCobj = prcomp(t(M_combat_1plate), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T)

pdf(file="Reports/PC_Variation_by_batch_afterCombat1.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Platte))
anova(lm(Prin.comp$PC2~Prin.comp$Platte))
anova(lm(Prin.comp$PC3~Prin.comp$Platte))
anova(lm(Prin.comp$PC4~Prin.comp$Platte))
anova(lm(Prin.comp$PC5~Prin.comp$Platte))

anova(lm(Prin.comp$PC1~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC2~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC3~Prin.comp$Sentrix_ID)) #0.003847
anova(lm(Prin.comp$PC4~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC5~Prin.comp$Sentrix_ID))

anova(lm(Prin.comp$PC1~Prin.comp$pos96well))
anova(lm(Prin.comp$PC2~Prin.comp$pos96well))
anova(lm(Prin.comp$PC3~Prin.comp$pos96well)) #< 2.2e-16 ***
anova(lm(Prin.comp$PC4~Prin.comp$pos96well))
anova(lm(Prin.comp$PC5~Prin.comp$pos96well))

#second round
mod <- model.matrix(~Case.Control+risk_all, data=pd_clean)
M_combat_2array = ComBat(M_combat_1plate,batch = pd_clean$Sentrix_ID, mod = mod)
save(M_combat_2array,file="RData/M_combat_2array.Rdata")

PCobj = prcomp(t(M_combat_2array), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
pd_clean = read.csv("pheno_and_samplesheet_clean_REPL.csv")
rownames(pd_clean) = pd_clean$MethID
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T)

# Check whether batches are still distinguished by first and second PC:
pdf(file="Reports/PC_Variation_by_batch_afterCombat2.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Platte))
anova(lm(Prin.comp$PC2~Prin.comp$Platte))
anova(lm(Prin.comp$PC3~Prin.comp$Platte))
anova(lm(Prin.comp$PC4~Prin.comp$Platte))
anova(lm(Prin.comp$PC5~Prin.comp$Platte))

anova(lm(Prin.comp$PC1~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC2~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC3~Prin.comp$Sentrix_ID)) #0.001
anova(lm(Prin.comp$PC4~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC5~Prin.comp$Sentrix_ID))

anova(lm(Prin.comp$PC1~Prin.comp$pos96well)) #< 2.2e-16 ***
anova(lm(Prin.comp$PC2~Prin.comp$pos96well)) #< 2.2e-16 ***
anova(lm(Prin.comp$PC3~Prin.comp$pos96well)) #< 2.2e-16 ***
anova(lm(Prin.comp$PC4~Prin.comp$pos96well))
anova(lm(Prin.comp$PC5~Prin.comp$pos96well)) #0.03601 *

pd_clean[is.na(pd_clean$risk_all),]$risk_all = 12
pd_clean = pd_clean[order(rownames(pd_clean)),]
mod <- model.matrix(~Case.Control+risk_all, data=pd_clean)
M_combat_3position = ComBat(M_combat_2array,batch = pd_clean$pos96well, mod = mod)
save(M_combat_3position,file="RData/M_combat_3position.Rdata")

PCobj = prcomp(t(M_combat_3position), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
pd_clean = read.csv("pheno_and_samplesheet_clean_REPL.csv")
rownames(pd_clean) = pd_clean$MethID
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T)

pdf(file="Reports/PC_Variation_by_batch_afterCombat3.pdf")
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

anova(lm(Prin.comp$PC1~Prin.comp$Platte))
anova(lm(Prin.comp$PC2~Prin.comp$Platte))
anova(lm(Prin.comp$PC3~Prin.comp$Platte))
anova(lm(Prin.comp$PC4~Prin.comp$Platte))
anova(lm(Prin.comp$PC5~Prin.comp$Platte))

anova(lm(Prin.comp$PC1~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC2~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC3~Prin.comp$Sentrix_ID)) #0.00117 **
anova(lm(Prin.comp$PC4~Prin.comp$Sentrix_ID))
anova(lm(Prin.comp$PC5~Prin.comp$Sentrix_ID))

anova(lm(Prin.comp$PC1~Prin.comp$pos96well))
anova(lm(Prin.comp$PC2~Prin.comp$pos96well))
anova(lm(Prin.comp$PC3~Prin.comp$pos96well))
anova(lm(Prin.comp$PC4~Prin.comp$pos96well))
anova(lm(Prin.comp$PC5~Prin.comp$pos96well))

#transform combat mvalues back
expit2 = function(x) 2^x/(1+2^x)
BMIQ.quantileN_combated = expit2(M_combat_3position)
dim(BMIQ.quantileN_combated) # 728868    636
save(BMIQ.quantileN_combated,file = "RData/BMIQ.quantileN_combated.Rdata")

png(file="Reports/BetaValue_Distributions_original_afterNorm_afterNormCombat3.png")
par(mfrow=c(3,1))
densityPlot(RawBetas_clean, sampGroups = pd_clean$Sentrix_ID, legend=FALSE, main = "Raw Betas", xlab = "Beta")
densityPlot(BMIQ.quantileN, sampGroups = pd_clean$Sentrix_ID, legend=FALSE, main = "quantileBMIQ adjusted Betas", xlab = "Beta")
densityPlot(BMIQ.quantileN_combated, sampGroups = pd_clean$Sentrix_ID, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off()

all.equal(colnames(BMIQ.quantileN_combated),rownames(pd_cleaned)) #TRUE
annotated_pd_clean =new("AnnotatedDataFrame", data= as.data.frame(pd_clean))
BMIQ.quantileN_combated_ExprSet = new("ExpressionSet", exprs= as.matrix(BMIQ.quantileN_combated), phenoData=annotated_pd_clean)
save(BMIQ.quantileN_combated_ExprSet,file="RData/BMIQ.quantileN_combated_ExprSet.Rdata")

genotypemethylationcoupling = cbind(pd_clean$Proben_name, pd_clean$MethID)
colnames(genotypemethylationcoupling) = c("genomeID","arrayID")


#create files for mixup mapper
write.table(BMIQ.quantileN_combated, file="files_for_mixup/Betas_combat_mixedupmapper.txt", sep="\t", quote=F, row.names=T, col.names=T)
# formating for mixupmapper
# sed -i "s/202193490115_R01C01/\t202193490115_R01C01/g"  Betas_combat_mixedupmapper.txt
# change plink to TRITYPER format
# java -jar /home/jade/src/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar -i ./data/LMU_genotypes -I PLINK_BED -o harmonized_geno -O TRITYPER -cf 0.95 -hf 0.0001 -mf 0.01 -mrf 0.5 -ip 0
# java -jar /home/jade/src/Mixup_Mapper/eqtl-mapping-pipeline.jar --mode normalize --in Betas_combat_mixedupmapper.txt --out finalbeta_norm --centerscale
# java -Xmx30g -Xms30g -jar /home/jade/src/Mixup_Mapper/eqtl-mapping-pipeline.jar --mode mixupmapper --in /binder/jade/KSP-LMU/genotype/GenotypeHarmonizer/harmonized_geno/ --out mixup_results --inexp files_for_mixup/finalbeta_norm/Betas_combat_mixedupmapper.ProbesCentered.txt.gz --inexpplatform EPIC --inexpannot files_for_mixup/annotation.txt --gte files_for_mixup/genotypemethylationcoupling.txt 2>&1 | tee mixup_results/mixupmapping.log

#cellcomposition
require(minfi)
require(FlowSorted.Blood.450k) #error with FlowSorted.Blood.EPIC
load("RData/RGset_clean_REPL.rda")
RGset_converted = convertArray(RGSet_clean, outType="IlluminaHumanMethylation450k")
cellcomps = estimateCellCounts(RGSet_converted, compositeCellType = "Blood", referencePlatform = "IlluminaHumanMethylation450k")
write.table(cellcomps, file="cellcomposition/LMU_houseman_cellcounts.csv",col.names=T,row.names=T,quote=F,sep=",")
