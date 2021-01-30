# 05_batch_effect

################################################################################
# control for batch effects with combat          
################################################################################
##### This is added to the pipeline
# We use combat to remove batch effects
# see e.g. Wen Bin Goh et al. 2017 for importance of batch effects removement 

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])
# or
input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

source(packages.fn)
source("functions.R")

load(pd_clean.fn)
load(bmiq.quantileN.fn)


##--- Combat to remove batch effects

mval <- apply(BMIQ.quantileN, 2, function(x) log2((x)/(1-x))) # M values

#-- Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars <- as.matrix(rowVars(mval))
which(vars == 0) # 0

#-- Replace all probes with no variance with NA and remove them from the normalized data set
vars[vars == 0] <- NA # 0
vars            <- na.omit(vars)
intersect       <- intersect(rownames(vars), rownames(mval))
print(length(intersect)) # probes without variance == 0

BMIQ.quantileN_batch <- BMIQ.quantileN[intersect, ]
mval                 <- mval[intersect,]

#-- Ensure Objects are aligned
table(ifelse(rownames(pd_clean) == colnames(mval),"Match","Off")) # All should match

## Check variation in array data associated with batch (ie. Slide/plate/box)
## Run a principle component analysis to determine if there are any remaining batch effects following data normalization.

PCobj <- prcomp(t(mval), retx = T, center = T, scale. = T)

# Scree plot for number of PCs

pdf(paste0(report.dir, "screeplot_PCA.pdf"))
fviz_eig(PCobj, ncp = 20)
dev.off()


# Extract the proportion of variability and cumulative proportion of varibility explained by the top R PCs

R <- 4
propvar <- round(summary(PCobj)$importance["Proportion of Variance", 1:R] * 100, 2 ) # -> write in xlsx
cummvar <- round(summary(PCobj)$importance["Cumulative Proportion", 1:R] * 100, 2) # -> write in xlsx
t(propvar); t(cummvar)

R <- 4

##---  Plot of PCA individal map by Batch, group (dex, veh) and sex 

PCs <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PC_Variation_by_batch_before_correction.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA before correction by Plate")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA before correction by Slide")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA before correction by Array")
dev.off()

##--- Find extreme outliers

o1 <- 1.5 * sd(Prin.comp$PC1)
o2 <- 1.5 * sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2) # 0

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Sample_Plate) 
})

anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
t(data.frame(anova.plate.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.           5.684620e-01
# PC2.Pr..F.           1.111118e-10
# PC3.Pr..F.           3.841619e-10
# PC4.Pr..F.           9.813305e-01

#-- for Slide

models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(Prin.comp$Slide))
})

anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
t(data.frame(anova.slide.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.  9.256443e-01
# PC2.Pr..F.  8.592451e-28
# PC3.Pr..F.  2.135583e-51
# PC4.Pr..F.  4.661745e-06

#-- for Array

models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})

anova.array.tbl <- sapply(models.array, anova, simplify = F)
t(data.frame(anova.array.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.    6.916667e-01
# PC2.Pr..F.    1.389933e-10
# PC3.Pr..F.    8.727043e-07
# PC4.Pr..F.    2.173907e-03

##--- Batch correction

#-- 1. Combat Correction for plate

# Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
model.mtrx <- model.matrix(~1, data = pd_clean)

# Run ComBat to remove most significant batch effects itertively
M_combat_1plate <- ComBat(mval, batch = pd_clean$Sample_Plate, mod = model.mtrx)
# save(M_combat_1plate, file = paste0(src.data.dir, "M_combat_1plate.Rdata"))

# Check to see if batch effect was successfully removed
PCobj     <- prcomp(t(M_combat_1plate), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PCA_ComBat_Correction_on_plate.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA ComBat correction on plate")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA ComBat correction on plate")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA ComBat correction on plate")
dev.off()

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Sample_Plate) 
})
anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
t(data.frame(anova.plate.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.            0.365975675
# PC2.Pr..F.            0.000145862
# PC3.Pr..F.            0.784264780
# PC4.Pr..F.            0.210116932

#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
t(data.frame(anova.slide.tbl)[1, c(5, 10, 15, 20)])
 
# PC1.Pr..F.                             9.691042e-01
# PC2.Pr..F.                             2.145165e-04
# PC3.Pr..F.                             4.931034e-05
# PC4.Pr..F.                             2.159736e-10

#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
t(data.frame(anova.array.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.    7.887227e-01
# PC2.Pr..F.    2.925067e-21
# PC3.Pr..F.    7.016281e-03
# PC4.Pr..F.    8.761940e-05

#-- 2. Combat Correction for Slide

mod             <- model.matrix(~1, data=pd_clean)
M_combat_2slide <- ComBat(M_combat_1plate,batch = pd_clean$Slide, mod = mod)
# save(M_combat_2slide, file = m.combat.2slide.fn)       

# Check to see if batch effect was succesfully removed
PCobj     <- prcomp(t(M_combat_2slide), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PCA_ComBat_Correction_on_plate_slide.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA after ComBat correction on both Plate and Slide")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA after ComBat correction on both Plate and Slide")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA after ComBat correction on both Plate and Slide")
dev.off()

o1 <- 1.5 * sd(Prin.comp$PC1)
o2 <- 1.5 * sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2) # 0

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Sample_Plate) 
})
anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
t(data.frame(anova.plate.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.  0.6299620
# PC2.Pr..F.  0.3854852
# PC3.Pr..F.  0.8519780
# PC4.Pr..F.  0.7124605

#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
t(data.frame(anova.slide.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.  0.99964921
# PC2.Pr..F.  0.99992277
# PC3.Pr..F.  0.07359928
# PC4.Pr..F.  0.14081736

#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
t(data.frame(anova.array.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.    6.654944e-01
# PC2.Pr..F.    1.180811e-32
# PC3.Pr..F.    2.159084e-02
# PC4.Pr..F.    5.682310e-02

#-- 3. Combat Correction for Array

mod             <- model.matrix(~1, data = pd_clean)
M_combat_3array <- ComBat(M_combat_2slide, batch = pd_clean$Array, mod = mod)
# save(M_combat_3array, file = m.combat.3array.fn)       

# Check to see if batch effect was succesfully removed
PCobj     <- prcomp(t(M_combat_3array), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PCA_ComBat_Correction_on_plate_slide_array.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA after final Combat correction")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA after final Combat correction")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA after final Combat correction")
dev.off()

o1 <- sd(Prin.comp$PC1)
o2 <- sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2) # 0

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Sample_Plate) 
})
anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
t(data.frame(anova.plate.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.              0.6193071
# PC2.Pr..F.              0.2321777
# PC3.Pr..F.              0.8041418
# PC4.Pr..F.              0.4536821

#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
t(data.frame(anova.slide.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.                               0.99957009
# PC2.Pr..F.                               0.85802028
# PC3.Pr..F.                               0.09660401
# PC4.Pr..F.                               0.16911922

#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
t(data.frame(anova.array.tbl)[1, c(5, 10, 15, 20)])

# PC1.Pr..F.       0.9614062
# PC2.Pr..F.       0.1781217
# PC3.Pr..F.       0.9595276
# PC4.Pr..F.       0.8281536

##--- Convert the batch-adjusted M-values back into betas:

expit2                      <- function(x) 2 ^ x / (1 + 2 ^ x)
betas.bmiq.quantilen.combat <- expit2(M_combat_3array)
dim(Betas_combated) # 865859    403

# Plot final densities
pdf(file = paste0(report.dir, "BetaValue_Distributions_after_ComBat.pdf"))
densityPlot(betas.bmiq.quantilen.combat, sampGroups = pd_clean$Slide, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Slide", xlab = "Beta")
densityPlot(betas.bmiq.quantilen.combat, sampGroups = pd_clean$Slide, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Slide", xlab = "Beta")
densityPlot(betas.bmiq.quantilen.combat, sampGroups = pd_clean$Array, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Array", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated), rownames(pd_clean)) #TRUE

annotated_pd_clean      <- new("AnnotatedDataFrame", data = as.data.frame(pd_clean)) #extend to AnnotatedDataFrame (required for ESet)
bmiq.quantilen.expr.set <- new("ExpressionSet", exprs = as.matrix(betas.bmiq.quantilen.combat), phenoData = annotated_pd_clean)
saveRDS(bmiq.quantilen.expr.set, file = paste0(src.data.dir, "bmiq_quantilen_expr_set.rds"))

# Save normalized and batch-adjusted beta values
saveRDS(betas.bmiq.quantilen.combat, file = paste0(src.data.dir, "betas_bmiq_quantilen_combat.rds"))
