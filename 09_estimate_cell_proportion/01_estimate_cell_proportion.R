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

# x                    <- load(beta.combat.expr.set.fn)
# beta.combat.expr.set <- get(x)
rgset.fn    <- paste0(src.final.data.dir, "dex_methyl_rgset_final.rds")
rgset       <- readRDS(rgset.fn)

# 2. Estimate cell type proportion

# Error with FlowSorted.Blood.EPIC
# Could not find reference data package for compositeCellType 'Blood' and referencePlatform 'EPIC' (inferred package name is 'FlowSorted.Blood.EPIC')

rgset.converted <- convertArray(rgset, outType="IlluminaHumanMethylation450k")

pdf(file = paste0(report.dir, "Cell_type_estimation.pdf"))
cell.counts            <- estimateCellCounts(rgset.converted,
                                             compositeCellType = "Blood",
                                             referencePlatform = "IlluminaHumanMethylation450k", 
                                             meanPlot = T)

dev.off()

write.table(cell.counts, file = paste0(report.dir, "DEX-stim-array-human-cellcounts.csv"), col.names = T, row.names = T, quote = F, sep = ";")
