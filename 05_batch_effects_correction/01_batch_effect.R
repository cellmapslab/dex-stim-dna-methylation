# 05 Batch effects correction

# 1. Load data

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

# OR
input.parameters.fn <- "input_parameters.csv"
  
input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

source(packages.fn)
source("05_batch_effects_correction/00_functions.R")

load(pd_clean.fn)
load(bmiq.quantileN.filtered.fn)

# 2. PCA 

mval <- apply(BMIQ.quantileN_filtered, 2, function(x) log2((x)/(1-x))) # M values

# Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars <- as.matrix(rowVars(mval))
which(vars == 0) # 0

# Replace all probes with no variance with NA and remove them from the normalized data set
vars[vars == 0] <- NA # 0
vars            <- na.omit(vars)
intersect       <- intersect(rownames(vars), rownames(mval))
print(length(intersect)) # probes without variance == 0

BMIQ.quantileN_filtered_batch <- BMIQ.quantileN_filtered[intersect, ]
mval                          <- mval[intersect,]

# Ensure Objects are aligned
table(ifelse(rownames(pd_clean) == colnames(mval),"Match","Off")) # All should match

# Check variation in array data associated with batch (ie. Slide/plate/box)
# Run PCA to determine if there are any remaining batch effects following data normalization

pc.obj <- prcomp(t(mval), retx = T, center = T, scale. = T)
save(pc.obj, file = pcobj.fn)

x      <- load(pcobj.fn)
pc.obj <- get(x)
rm(x)

# Scree plot to determine number of PCs to keep
pdf(paste0(report.dir, "screeplot_PCA.pdf"))
# plot(PCobj, type="line",cex.lab = 1.5, cex.main = 1.5) 
fviz_eig(pc.obj, ncp = 20)
dev.off()

# Extract the proportion of variability and cumulative proportion of 
# varibility explained by the top R PCs.
R <- 15
propvar <- round(summary(pc.obj)$importance["Proportion of Variance", 1:R] * 100, 2 ) # -> write in xlsx
cummvar <- round(summary(pc.obj)$importance["Cumulative Proportion", 1:R] * 100, 2) # -> write in xlsx
t(propvar); t(cumvar)

# Choose n first PCs that explains the most variation
R <- 6 # choose eg 6 PCs
prin.comp <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

# Find extreme outliers

o1 <- 1.5 * sd(prin.comp$PC1)
o2 <- 1.5 * sd(prin.comp$PC2)
which(abs(prin.comp$PC1) > o1 && abs(prin.comp$PC2) > o2) # 0

# PCA and ANOVA result before correction 

pdf.fn <- paste0(report.dir, "00_PCA-map_ANOVA-res_before_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

# 3. ComBat correction for Plate

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.plate <- ComBat(mval, batch = pd_clean$Sample_Plate, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.plate), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "01_PCA-map_ANOVA-res_PLATE_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

# 4. ComBat correction for Array

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.array <- ComBat(m.combat.plate, batch = pd_clean$Array, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.array), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "02_PCA-map_ANOVA-res_ARRAY_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

# 5. ComBat correction for Slide

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.slide <- ComBat(m.combat.array, batch = pd_clean$Slide, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.slide), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "03_PCA-map_ANOVA-res_SLIDE_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

# 6. Convert the batch-adjusted M-values back into betas:

m.combat.final <- m.combat.slide 

expit2         <- function(x) 2 ^ x / (1 + 2 ^ x)
Betas_combated <- expit2(m.combat.final)
dim(Betas_combated) # for now final betas!
# 740357    403

# 7. Plot final densities

pdf(file = paste0(report.dir, "04_BetaComBated_Distributions_Plot.pdf"))
densityPlot(Betas_combated, 
            sampGroups = pd_clean$Plate, 
            legend = FALSE, 
            main = "PostQC - Normalized and Batch Corrected Beta by Plate", xlab = "Beta")
densityPlot(Betas_combated, 
            sampGroups = pd_clean$Array, 
            legend = FALSE, 
            main = "PostQC - Normalized and Batch Corrected Beta by Array", xlab = "Beta")
densityPlot(Betas_combated, 
            sampGroups = pd_clean$Slide, 
            legend = FALSE, 
            main = "PostQC - Normalized and Batch Corrected Beta by Slide", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated), rownames(pd_clean)) #TRUE

# 8. Extend to AnnotatedDataFrame (required for ESet)

annotated_pd_clean     <- new("AnnotatedDataFrame", data= as.data.frame(pd_clean)) 
Betas_combated_ExprSet <- new("ExpressionSet", exprs = as.matrix(Betas_combated), phenoData = annotated_pd_clean)

# 9. Save normalized and batch-adjusted beta values and expr. set

save(Betas_combated, file = beta.combat.fn)
save(Betas_combated_ExprSet, file = beta.combat.expr.set.fn)

#--- END
