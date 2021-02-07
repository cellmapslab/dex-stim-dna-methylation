# 09 Cell Composition

# 1. Set up env and load data

library(minfi)

library(FlowSorted.Blood.EPIC)
library(FlowSorted.Blood.450k) 

library(ggplot2)
library(RColorBrewer)

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# x                    <- load(beta.combat.expr.set.fn)
# beta.combat.expr.set <- get(x)
rgset.fn    <- paste0(src.final.data.dir, "dex_methyl_rgset_final.rds")
rgset       <- readRDS(rgset.fn)

pheno.fn     <- paste0(src.final.data.dir, "dex_methyl_phenotype.rds")
pheno        <- readRDS(pheno.fn)

# 2. Estimate cell type proportion

# Error with FlowSorted.Blood.EPIC
# Could not find reference data package for compositeCellType 'Blood' and referencePlatform 'EPIC' (inferred package name is 'FlowSorted.Blood.EPIC')

rgset.converted <- convertArray(rgset, outType="IlluminaHumanMethylation450k")

pdf(file = paste0(report.dir, "dex_stim_cell_type_estimation_plot.pdf"))
cell.counts            <- estimateCellCounts(rgset.converted,
                                             compositeCellType = "Blood",
                                             referencePlatform = "IlluminaHumanMethylation450k", 
                                             meanPlot = T)

dev.off()

write.table(cell.counts, file = paste0(report.dir, "dex_stim_array_human_cellcounts.csv"), col.names = T, row.names = T, quote = F, sep = ";")

# 3. Plot cell type proportions by group (dex, veh)

pal <- brewer.pal(8, "Set2")

pdf(file = paste0(report.dir, "dex_stim_cell_type_estimation_boxplot_by_group.pdf"))
cell.counts.dex <- cell.counts[pheno$Sample_Group == "dex",]
cell.counts.veh <- cell.counts[pheno$Sample_Group == "veh",]
boxplot(cell.counts.dex, at = 0:5 * 3 + 1, 
        xlim = c(0, 18), ylim = range(cell.counts.dex, cell.counts.veh), 
        xaxt = "n", main = "", ylab = "Cell type proportion", col = pal[1])
boxplot(cell.counts.veh, at = 0:5 * 3 + 2, xaxt = "n", col = pal[2], add = TRUE)
axis(1, at = 0:5 * 3 + 1.5, labels = colnames(cell.counts.dex), tick = TRUE)
legend("topleft", legend = c("dex", "veh"), fill = pal)
dev.off()