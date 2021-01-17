# 04_filtering

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

source(packages.fn)

load(detP_clean.fn)
load(quantileN.bmiq.fn)
load(quantileN.fn)
# load(rawBetas_clean.fn)
# load(pd_clean.fn)
load(annotated_data_clean.fn)

##--- Filter out probes that failed in one or more samples based on detection p < 0.01

#--- BMIQ quantileN
# ensure probes are in the same order in the BMIQ beta and detP objects
detP_clean_f <- detP_clean[match(rownames(BMIQ.quantileN), rownames(detP_clean)),]

# remove any probes that have failed in one or more samples
keep <- rowSums(detP_clean_f < 0.01) == ncol(BMIQ.quantileN) 
table(keep) # -> note output (excel)
# FALSE   TRUE
# 33269 832590

# exclude these probes from our normalized data
BMIQ.quantileN_filtered <- BMIQ.quantileN[keep,]
# matrix
dim(BMIQ.quantileN_filtered)

#--- quantileN
# ensure probes are in the same order in the quantileN (GenomicRatio) and detP objects
detP_clean_ff <- detP_clean[match(featureNames(quantileN), rownames(detP_clean)),]
keep_ff <- rowSums(detP_clean_ff < 0.01) == ncol(quantileN) 
table(keep_ff) 
# FALSE   TRUE
# 33269 832590

# exclude these probes from our normalized data
quantileN_filtered <- quantileN[keep_ff,]
# genomicRatioSet
quantileN_filtered

##--- Filter out probes on sex chromosomes

#--- BMIQ quantileN
keep <- !(rownames(BMIQ.quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
table(keep) # -> note output (excel)
# FALSE   TRUE
# 17665 814925
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,]
dim(BMIQ.quantileN_filtered)

#--- quantileN
keep_ff <- !(featureNames(quantileN_filtered) %in% annot$Name[annot$chr %in% c("chrX","chrY")])
# FALSE   TRUE
# 17665 814925
quantileN_filtered <- quantileN_filtered[keep_ff,]
quantileN_filtered

##### removal of probes where common SNPs may affect the CpG. 
# You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.
quantileN_filtered <- dropLociWithSnps(quantileN_filtered) # this NEEEDs a GenomicRatioSet, beta matrix does not work
quantileN_filtered # -> note again dims
# 
BMIQ.quantileN_filtered <-  BMIQ.quantileN_filtered[rownames(BMIQ.quantileN_filtered) %in% featureNames(quantileN_filtered),]
dim(BMIQ.quantileN_filtered)
# features = 788960    samples = 403

# load Chen Probe annotation file and exclude SNP and Cross Hybridizing probes:
load("ChenProbeIDs.rdata") #### Data from Chen et al (2013) PMID: 23314698
# loaded as annot2
annot2$SNPs = annot2[,"EUR"] #### Can change Global to "AFR", "AMR", "ASN", or "EUR" to match content specific allelic frequency; Use "Global" if population is very diverse
index<-which(annot2$sex=="Exclude" | annot2$CRSH=="Exclude" | annot2$EUR=="Exclude")
length(index) #62520
exclude_Chen <-annot2[index,]

#not in pipeline, but additionally use data from McCartney et al which are based on EPIC# PMID: 27330998
exclude_crosshyb <-read.table("addFiles/CpGs_crosshybridizing_EPIC.txt",sep="",header=F)
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep) # -> note output
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)
#
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep_ff) # 
quantileN_filtered <- quantileN_filtered[keep_ff,] 
quantileN_filtered # 

exclude_poly <-read.table("addFiles/CpGs_polymorphic_EPIC.txt",sep="",header=T)
index<-which(exclude_poly$EUR_AF>0.05) #n=10971
exclude_polym <- exclude_poly[index,]

keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_polym$IlmnID) 
table(keep) # -> note output
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_polym$IlmnID) 
table(keep_ff) 
quantileN_filtered <- quantileN_filtered[keep_ff,] 

### ### ###
save(BMIQ.quantileN_filtered, file="RData/BMIQ.quantileN_filtered.Rdata") # QN+BMIQ normalized, probe filtered & sample-qced Betas.
save(quantileN_filtered, file="RData/quantileN_filtered.Rdata") # GenomicRatio of filtered QN normalized data
Betas_quantileN_filtered <- getBeta(quantileN_filtered) # beta values
save(Betas_quantileN_filtered, file="RData/Betas_quantileN_filtered.Rdata")
Ms_quantileN_filtered <- getM(quantileN_filtered)
save(Ms_quantileN_filtered, file="RData/Ms_quantileN_filtered.Rdata")

## Check density plots after excluding the poorly-detected probes:
png(file="Reports/BetaValue_Distributions_Norm_quantile_Filter.png")
densityPlot(BMIQ.quantileN_filtered,sampGroups = pd_clean$Slide, legend=FALSE, main = "Post Filtering - Normalized Beta", xlab = "Beta")
dev.off() 
