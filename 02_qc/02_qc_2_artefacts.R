src.data.dir   <- '/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData/'
src.dir        <- '~/github/dex-stim-dna-methylation/'
report.dir     <- paste0(src.dir, "02_qc/reports/")

packages.fn    <- paste0(src.dir, "packages.R")
source(packages.fn)

rgSet_qc.fn    <- paste0(src.data.dir, "rgSet_qc.Rdata")
load(rgSet_qc.fn)

rgSet.fn       <- paste0(src.data.dir, "rgSet_dex.RData")
load(rgSet.fn)

rawBetas_qc.fn <- paste0(src.data.dir, "RawBetas_qc.Rdata")
load(rawBetas_qc.fn)  


##--- Distribution artefacts

#-- check individual densities
pdf(paste0(report.dir, "03_beta_densities_report.pdf"))

for (i in 1:ncol(rgSet_qc)){
	title <- paste(rownames(pData(rgSet_qc))[i])
	densityPlot(as.matrix(RawBetas_qc[,i]), main = title)
	print(i)}

dev.off()


##--- Sex mismatches

detP.fn	<- paste0(src.data.dir, "detP.Rdata")
load(detP.fn)

detP_qc.fn <- paste0(src.data.dir, "detP_ql.Rdata")
load(detP_qc.fn)

gRatio.fn <- paste0(src.data.dir, "gRatioSet_original.Rdata")
load(gRatio.fn)

pd.fn <-  paste0(src.data.dir, "pd_original.Rdata")
load(pd.fn)

#-- predict sex with methylation data
predictedSex <-getSex(gRatioSet, cutoff = -2)
sex          <- cbind(sampleNames(rgSet),
		as.character(pd.orig$Sample_Name),
		predictedSex$predictedSex,
		predictedSex$yMed,
		predictedSex$xMed,
		predictedSex$yMed-predictedSex$xMed)
sex          <- as.data.frame(sex)
names(sex)   <- c("ArrayID", "ID", "predictedSex", "yMed", "xMed", "yMed-xMed")
sex$yMed     <-as.numeric(as.character(sex$yMed))
sex$xMed     <-as.numeric(as.character(sex$xMed))
sex$ID       <- as.character(sex$ID)
save(sex, file = paste0(src.data.dir, "sex_predicted.Rdata"))


#-- get real sex

sample.sheet.fn <- paste0(src.dir, "00_samplesheets/samplesheet_epic_methylation_dex.csv")
targets         <- read.csv(sample.sheet.fn, sep = ';' )
gender          <- targets[, c("Sample_Name", "sex")]

sex.test.df              <- merge(sex, gender, by.x = "ID", by.y = "Sample_Name")
colnames(sex.test.df)[7] <- "realSex"
sex.test.df$realSex      <-  gsub("W", "F", sex.test.df$realSex)
(sex.match.ids           <- subset(sex$ID, !(sex$ID %in% gender$Sample_Name))) # no matched --> good

write.table(sex.test.df, paste0(src.data.dir, "sex_prediction_tbl.txt"), sep = "\t", quote = F, row.names = F)

table(sex.test.df$predictedSex, sex.test.df$realSex) # 135 FF, 268 MM, 1 realFpredM, total = 404

sex.test.df[sex.test.df$realSex == "F" & sex.test.df$predictedSex == "M", ]
# ID                      ArrayID     predictedSex    yMed      xMed         yMed-xMed       realSex
# 200705940062_R06C01 200705940062_R06C01   M       9.390169  11.0858   -1.69563538590411       F

# Checked with genotype data:
# MPIPSYKL_009498 MPIPSYKL_009498 0 0 2 -9 

##--- Sample exclusions:
  # 1. Distribution artefacts
  # 2. Poor quality (detection p-value)
  # 3. Sex mismatches

##-- Individuals to exclude (BEFORE NORMALIZATION) :
  # DONE: Already excluded samples based on detection p-value
  # TODO: Exclude the distribution artefacts and sex mismatches from the files


out.ids <- c("200705940062_R06C01") # sample to exclude

out.ids <- out.ids[!duplicated(out.ids)] 
i1      <- intersect(out.ids, colnames(rgSet_qc)) # see how many are still in RGSet_qual
i2      <- intersect(out.ids, colnames(rgSet)) # re-check that all identified IDs are in the full RGset (no misspelling)
setdiff(i2, i1) # ID(s) identified in detP step and again in later quality control

rgSet_clean    <- rgSet_qc[, !(colnames(rgSet_qc) %in% out.ids)]
rgSet_clean
save(rgSet_clean, file = paste0(src.data.dir, "rgSet_clean.Rdata"))

RawBetas_clean <-  RawBetas_qc[, !(colnames(RawBetas_qc) %in% out.ids)]
save(RawBetas_clean, file = paste0(src.data.dir, "RawBetas_clean.Rdata")) 

detP_clean     <- detP_ql[, !(colnames(detP_ql) %in% out.ids)] 
save(detP_clean, file = paste0(src.data.dir, "detP_clean.Rdata")) 

pd_clean       <- pData(rgSet_clean)
save(pd_clean, file = paste0(src.data.dir, "pd_clean.Rdata")) 

##--- Get annotations probes:
annot          <- getAnnotation(rgSet_clean)
save(annot, file = paste0(src.data.dir, "annotated_data_clean.Rdata"))
