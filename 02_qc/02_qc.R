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
detP       <- detectionP(rgSet)
head(detP)
save(detP, file = paste0(src.data.dir, "detP.Rdata")

load(pvalues.fn)

##--- Examine mean detection p-values across all samples to identify any failed samples

pdf(paste(report.dir, "detectionP.pdf", sep="/"), paper = "a4", width = 11, height = 8)
pal <- brewer.pal(12,"Set3")
par(mfrow = c(2,))

barplot(colMeans(detP), 
	col = pal[factor(targets$BeadChipPosition)], 
	las = 2, 
	cex.names = 0.2, 
	ylab = "Mean detection p-values")
abline(h = 0.01, col = "red")
# legend("topleft", legend=levels(factor(targets$BeadChipPosition)), fill = pal, bg = "white")

barplot(colMeans(detP), 
	col = pal[factor(targets$Sample_Group)], 
	las = 2, 
	cex.names = 0.4, 
	ylim = c(0,0.002), 
	ylab = "Mean detection p-values")
# abline(h = 0.001, col = "red")
# legend("topleft", legend = levels(factor(targets$BeadChipPosition)), fill = pal, bg = "white")

dev.off()
