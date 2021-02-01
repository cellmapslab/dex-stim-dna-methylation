# 09 Cell Composition

# 1. Set up env and load data

library(minfi)

library(FlowSorted.Blood.EPIC)
library(FlowSorted.Blood.450k) 

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

x      <- load(beta.combat.expr.set.fn)
rg.set <- get(x)
rm(x)

# 2. Estimate cell type proportion

rg.set.converted <- convertArray(rg.set, outType="IlluminaHumanMethylation450k")
cell.counts      <- estimateCellCounts(rg.set.converted,
                                       compositeCellType = "Blood",
                                       referencePlatform = "IlluminaHumanMethylation450k")

# Error with FlowSorted.Blood.EPIC
# cell.counts      <- estimateCellCounts(rg.set, 
#                                       compositeCellType = "Blood", 
#                                       referencePlatform = "IlluminaHumanMethylationEPIC")

write.csv(cell.counts, file = paste0(report.dir, "DEX-stim-array-human-cellcounts.csv"), col.names = T, row.names = T, quote = F, sep = ";")
