src.data.dir  <- '/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData/'
report.dir    <- '/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/reports/"

rgSet.fn      <- paste0(src.data.dir, "rgSet_dex.RData")

source(packages.fn)

load(rgSet.fn)
rgSet <- rgSubSet

sample.sheet.fn <- paste0(src.dir, "samplesheet_epic_methylation_dex.csv")
targets         <- read.csv(sample.sheet.fn, sep = ',' )


##--- Calculate the detection p-values

detP <- detectionP(rgSet)
head(detP)
save(detP, file = "RData/detP.Rdata") 


##--- Examine mean detection p-values across all samples to identify any failed samples

pdf(paste0(report.dir, "detectionP.pdf"), width = 8, height = 3)
pal <- brewer.pal(8,"Dark2")
par(mfrow = c(1,2))

barplot(colMeans(detP), 
	col = pal[factor(targets$Sample_Group)], 
	las = 2, cex.names = 0.8, ylab = "Mean detection p-values")
abline(h = 0.05, col = "red")
legend("topleft", 
	legend=levels(factor(targets$Sample_Group)), fill=pal,
        bg="white")

barplot(colMeans(detP), col=pal[factor(targets$Sample_Group)], las=2, 
        cex.names=0.8, ylim=c(0,0.002), ylab="Mean detection p-values")
abline(h=0.05,col="red")
legend("topleft", legend=levels(factor(targets$Sample_Group)), fill=pal, 
       bg="white")
