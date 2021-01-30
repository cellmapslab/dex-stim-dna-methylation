# 04_filtering

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

source(packages.fn)

help.files.dir <- paste0(src.dir, "04_filtering/01_help_files/")

load(detP_clean.fn)
load(quantileN.bmiq.fn)
load(quantileN.fn)
load(annotated_data_clean.fn)

##--- 1. Filter out probes that failed in one or more samples based on detection p < 0.01

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

##--- 2. Filter out probes on sex chromosomes

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


##--- 3. Removal of probes where common SNPs may affect the CpG

# You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.
quantileN_filtered <- dropLociWithSnps(quantileN_filtered) # this NEEEDs a GenomicRatioSet, beta matrix does not work
quantileN_filtered # -> note again dims
# 
BMIQ.quantileN_filtered <-  BMIQ.quantileN_filtered[rownames(BMIQ.quantileN_filtered) %in% featureNames(quantileN_filtered),]
dim(BMIQ.quantileN_filtered)
# features = 788960    samples = 403


##--- 4. Load Chen Probe annotation file and exclude SNP and Cross Hybridizing probes:

#--- Cross Hybridizing probes

cg.croshyb.epic.fn <- paste0(help.files.dir, CpGs_crosshybridizing_EPIC.txt")
exclude_crosshyb   <-read.table(cg.croshyb.epic.fn, sep = "", header = F)

#--- BMIQ quantileN
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep) # -> note output
# FALSE   TRUE
# 37869 751091
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)

#--- quantileN
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_crosshyb$V1) 
table(keep_ff)

quantileN_filtered <- quantileN_filtered[keep_ff,] 
quantileN_filtered # dim: 751091 403

#--- Polymorphic probes

cg.polymorphic.epic.fn <- paste0(help.files.dir, CpGs_polymorphic_EPIC.txt")
exclude_poly           <- read.table(cg.polymorphic.epic.fn, sep = "", header = T)
index                  <- which(exclude_poly$EUR_AF > 0.05) #n = 10971
exclude_poly.ids       <- exclude_poly[index,]

#--- BMIQ quantileN
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_poly.ids$IlmnID) 
table(keep) # -> note output
#  FALSE   TRUE
# 447 750644
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)

#--- quantileN
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_poly.ids$IlmnID) 
table(keep_ff) 
quantileN_filtered <- quantileN_filtered[keep_ff,] 

#--- Chen Probe annotation

chen.probe.id.fn <- paste0(help.files.dir, "ChenProbeIDs.rdata")
load(chen.probe.id.fn) #### Data from Chen et al (2013) PMID: 23314698

# loaded as annot2
annot2$SNPs  <- annot2[, "EUR"] #### Can change Global to "AFR", "AMR", "ASN", or "EUR" to match content specific allelic frequency; Use "Global" if population is very diverse
index        <- which(annot2$sex == "Exclude" | annot2$CRSH == "Exclude" | annot2$EUR == "Exclude")
length(index) # 62520
exclude_Chen <- annot2[index, ]

#--- BMIQ quantileN
keep <- !(rownames(BMIQ.quantileN_filtered) %in% exclude_Chen$Name) 
table(keep) # -> note output
# FALSE   TRUE
# 10287 740357
BMIQ.quantileN_filtered <- BMIQ.quantileN_filtered[keep,] 
dim(BMIQ.quantileN_filtered)

#--- quantileN
keep_ff <- !(featureNames(quantileN_filtered) %in% exclude_Chen$Name) 
table(keep_ff) 
quantileN_filtered <- quantileN_filtered[keep_ff,] 

#not in pipeline, but additionally use data from McCartney et al which are based on EPIC# PMID: 27330998

##--- 5. Save data

#--- QN+BMIQ normalized, probe filtered & sample-qced Betas
save(BMIQ.quantileN_filtered, file = "/home/ahryhorzhevska/mpip/datasets/methyl/rData/BMIQ.quantileN_filtered.Rdata") 

#--- GenomicRatio of filtered QN normalized data
save(quantileN_filtered, file = "/home/ahryhorzhevska/mpip/datasets/methyl/rData/quantileN_filtered.Rdata") 

##--- 6. Check density plots after excluding the poorly-detected probes
png(file = paste0(report.dir, "BetaValue_Distributions_Norm_quantile_Filter.png"))
densityPlot(BMIQ.quantileN_filtered,sampGroups = pd_clean$Slide, legend = FALSE, main = "Post Filtering - Normalized Beta", xlab = "Beta")
dev.off() 


Betas_quantileN_filtered <- getBeta(quantileN_filtered) # beta values
save(Betas_quantileN_filtered, file = "/home/ahryhorzhevska/mpip/datasets/methyl/rData/Betas_quantileN_filtered.Rdata")

Ms_quantileN_filtered <- getM(quantileN_filtered)
save(Ms_quantileN_filtered, file = "/home/ahryhorzhevska/mpip/datasets/methyl/rData/Ms_quantileN_filtered.Rdata")


