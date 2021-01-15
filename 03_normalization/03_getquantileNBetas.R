args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# print(rgSet_clean.fn)
source(packages.fn)
source(bmiq.script.fn)

load(quantileN.fn)
load(rawBetas_clean.fn)

quantileNBetas <-  getBeta(quantileN) # get betas
save(quantileNBetas, file = quantileNBetas.fn)

# BMIQ after quantile normalization:
probeType           <- as.data.frame(annot[rownames(quantileNBetas), c("Name","Type")])
probeType$probeType <- ifelse(probeType$Type %in% "I", 1, 2)

BMIQ.quantileN      <- apply(quantileNBetas[, 1:length(colnames(quantileNBetas))], 2, 
				function(a) BMIQ(a,probeType$probeType,plots = FALSE)$nbeta)

length(which(is.nan(BMIQ.quantileN))) # should be 0
save(BMIQ.quantileN, file = paste0(src.data.dir, "BMIQ.quantileN.Rdata")

## Check distributions before and after chosen normalization
png(file = paste0(report.dir, "BetaValue_Distributions_Norm_Quantile.png"), width = 1400, height = 700, pointsize = 12)
par(mfrow = c(1,2))
densityPlot(RawBetas_clean, sampGroups = pd_clean$Slide, legend=FALSE, main = "Raw Betas", xlab = "Beta")
densityPlot(quantileNBetas, sampGroups = pd_clean$Slide, legend=FALSE, main = "Quantile Adjusted Betas", xlab = "Beta")
densityPlot(BMIQ.quantileN, sampGroups = pd_clean$Slide, legend=FALSE, main = "Quantile-BMIQ Adjusted Betas", xlab = "Beta", ylim=c(0,5))
dev.off()
