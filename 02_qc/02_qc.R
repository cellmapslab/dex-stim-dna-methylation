src.data.dir <- '/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/rData/'
src.dir      <- '~/github/dex-stim-dna-methylation/'
report.dir   <- paste(getwd(), "reports", sep = "/")

sh.script.mkdir <- sprintf("cd %s ; 
                    
                    if [ ! -d reports ] 
                    then 
                      mkdir reports 
                    fi ; 
                    
                    cd reports ;
      
                    ", report.dir)

system(sh.script.mkdir)

packages.fn     <- paste0(src.dir, "packages.R")
source(packages.fn)

sample.sheet.fn <- paste0(src.dir, "samplesheet_epic_methylation_dex.csv")
targets         <- read.csv(sample.sheet.fn, sep = ',' )

rgSet.fn        <- paste0(src.data.dir, "rgSet_dex.RData")
load(rgSet.fn)

##--- Calculate the detection p-values

pvalues.fn <- paste0(src.data.dir, "detP.Rdata")
if(file.exists(pvalues.fn)) {
	load(pvalues.fn) } else {
	detP       <- detectionP(rgSet)
	head(detP)
	save(detP, file = pvalues.fn)
}

##--- Examine mean detection p-values across all samples to identify any failed samples

pdf(paste(report.dir, "detectionP.pdf", sep="/"), paper = "a4", width = 11, height = 8)
pal <- brewer.pal(12,"Set3")
par(mfrow = c(2,))

barplot(colMeans(detP), 
	col = pal[factor(targets$BeadChipPosition)], 
	las = 2, 
	xaxt = "n",
	border="transparent",
	xlab = "Samples", 
	ylab = "Mean detection p-values")
abline(h = 0.01, col = "red")
# legend("topleft", legend=levels(factor(targets$BeadChipPosition)), fill = pal, bg = "white")

barplot(colMeans(detP), 
	col = pal[factor(targets$BeadChipPosition)], 
	las = 2, 
	xaxt = "n", 
	border="transparent",
	ylim = c(0,0.002),
	xlab = "Samples" ,
	ylab = "Mean detection p-values")
abline(h = 0.001, col = "red")
# legend("topleft", legend = levels(factor(targets$BeadChipPosition)), fill = pal, bg = "white")

dev.off()


##---  Remove poor quality samples with a mean detection p-value _> 0.01_

# remove poor quality samples with a mean detection p-value >.05
keep     <- colMeans(detP) < 0.01 
rgSet_qc <- rgSet[, keep] # dim: 403 out of 404
save(rgSet_qc, file = paste0(src.data.dir, "rgSet_qc.Rdata")) 

# remove poor quality samples from raw betas
load(paste0(src.data.dir, "RawBetas.Rdata"))
rawBetas_ql <- RawBetas[, keep]
dim(RawBetas_ql)
save(rawBetas_ql, file = paste(src.data.dir, "rawBetas_ql.Rdata")) 

# remove poor quality samples from targets data = targets_qual
targets_qual <- targets[keep,] 
targets_qual[, 1:5] 

# remove poor quality samples from detection p-value table 
detP_ql <- detP[, keep] 
dim(detP_ql) 
save(detP_ql, file = paste0(src.data.dir, "detP_ql.Rdata"))

# -> note how many samples were excluded

# run minfi QC Report
qcReport(rgSet_ql, sampGroups = targets_ql$Slides, pdf = paste0(report.dir, "qcReport.pdf")) 
