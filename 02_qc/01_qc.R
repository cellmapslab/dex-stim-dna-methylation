args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

source(packages.fn)

sh.script.mkdir <- sprintf("cd %s ; 
                    
                    if [ ! -d reports ] 
                    then 
                      mkdir reports 
                    fi ; 
                    
                    cd reports ;
      
                    ", report.dir)

system(sh.script.mkdir)

targets         <- read.csv(sample.sheet.fn, sep = ';' )
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
par(mfrow = c(2, 1))

barplot(colMeans(detP), 
	col = pal[factor(targets$Sample_Group)], 
	las = 2, 
	xaxt = "n",
	border="transparent",
	xlab = "Samples", 
	ylab = "Mean detection p-values")
abline(h = 0.01, col = "red")
# legend("topleft", legend=levels(factor(targets$BeadChipPosition)), fill = pal, bg = "white")

barplot(colMeans(detP), 
	col = pal[factor(targets$Sample_Group)], 
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
# 200705940062_R06C01 to remove
rgSet_qc <- rgSet[, keep] # dim: 403 out of 404
save(rgSet_qc, file = paste0(src.data.dir, "rgSet_qc.Rdata")) 

# remove poor quality samples from raw betas
load(paste0(src.data.dir, "RawBetas_original.Rdata"))
RawBetas_qc <- RawBetas[, keep]
dim(RawBetas_qc)
save(RawBetas_qc, file = paste(src.data.dir, "RawBetas_qc.Rdata")) 

# remove poor quality samples from targets data = targets_qual
targets_qc <- targets[keep,] 
targets_qc[, 1:5] 

# remove poor quality samples from detection p-value table 
detP_ql <- detP[, keep] 
dim(detP_ql) # 866836 x 403 
save(detP_ql, file = paste0(src.data.dir, "detP_ql.Rdata"))

# -> note how many samples were excluded

# run minfi QC Report
names <- pData(rgSet_qc)$Sample_Name
groups <- pData(rgSet_qc)$Sample_Group
qcReport(rgSet_qc, sampNames = names, sampGroups = groups, 
		pdf = paste0(report.dir, "qcReport.pdf")) 
