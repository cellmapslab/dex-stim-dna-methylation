src.data.dir   <- '/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData/'
src.dir        <- '~/github/dex-stim-dna-methylation/'
report.dir     <- paste0(src.dir, "02_qc/reports/")

rgSet_qc.fn    <- paste0(src.data.dir, "rgSet_qc.Rdata")
load(rgSet_qc.fn)

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

gRatio.fn <- paste0(src.data.dir, "gRatioSet_original.Rdata")
load(gRatio.fn)
