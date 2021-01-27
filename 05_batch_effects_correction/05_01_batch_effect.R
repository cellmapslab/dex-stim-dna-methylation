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
source("05_batch_effects_correction/05_00_functions.R")

report.dir <- paste0(src.dir, "05_batch_effects_correction/reports/")

load(pd_clean.fn)
load(bmiq.quantileN.filtered.fn)

##--- PCA 

mval <- apply(BMIQ.quantileN_filtered, 2, function(x) log2((x)/(1-x))) # M values

#-- Calculate the variance of each probe and remove any with a variance of 0 prior to Combat.
vars <- as.matrix(rowVars(mval))
which(vars == 0) # 0

#-- Replace all probes with no variance with NA and remove them from the normalized data set
vars[vars == 0] <- NA # 0
vars            <- na.omit(vars)
intersect       <- intersect(rownames(vars), rownames(mval))
print(length(intersect)) # probes without variance == 0

BMIQ.quantileN_filtered_batch <- BMIQ.quantileN_filtered[intersect, ]
mval                          <- mval[intersect,]

#-- Ensure Objects are aligned
table(ifelse(rownames(pd_clean) == colnames(mval),"Match","Off")) # All should match

## Check variation in array data associated with batch (ie. Slide/plate/box)
## Run a principle component analysis to determine if there are any remaining batch effects following data normalization.

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

##--- Find extreme outliers

o1 <- 1.5 * sd(prin.comp$PC1)
o2 <- 1.5 * sd(prin.comp$PC2)
which(abs(prin.comp$PC1) > o1 && abs(prin.comp$PC2) > o2) # 0

##--- 

pdf.fn <- paste0(report.dir, "00_PCA-map_ANOVA-res_before_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

##--- ComBat correction for Plate

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.plate <- ComBat(mval, batch = pd_clean$Sample_Plate, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.plate), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "01_PCA-map_ANOVA-res_PLATE_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

##--- ComBat correction for Array

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.array <- ComBat(m.combat.plate, batch = pd_clean$Array, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.array), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "02_PCA-map_ANOVA-res_ARRAY_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

##--- ComBat correction for Slide

model.mtrx     <- model.matrix(~1, data = pd_clean)
m.combat.slide <- ComBat(m.combat.array, batch = pd_clean$Slide, mod = model.mtrx)
pc.obj         <- prcomp(t(m.combat.slide), retx = T, center = T, scale. = T)
prin.comp      <- merge(pc.obj$x[, 1:R], pd_clean, by = "row.names", all = T) 

pdf.fn <- paste0(report.dir, "03_PCA-map_ANOVA-res_SLIDE_correction.pdf")
GetPCAnovaReport(pc.obj, prin.comp, R, pdf.fn) 

##--- Convert the batch-adjusted M-values back into betas:

m.combat.final <- m.combat.slide 

expit2         <- function(x) 2 ^ x / (1 + 2 ^ x)
Betas_combated <- expit2(m.combat.final)
dim(Betas_combated) # for now final betas!
# 740357    403

#plot final densities
pdf(file = paste0(report.dir, "04_BetaComBated_Distributions_Plot.pdf"))
densityPlot(Betas_combated, sampGroups = pd_clean$Plate, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Plate", xlab = "Beta")
densityPlot(Betas_combated, sampGroups = pd_clean$Array, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Array", xlab = "Beta")
densityPlot(Betas_combated, sampGroups = pd_clean$Slide, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Slide", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated), rownames(pd_clean)) #TRUE

annotated_pd_clean     <- new("AnnotatedDataFrame", data= as.data.frame(pd_clean)) #extend to AnnotatedDataFrame (required for ESet)
Betas_combated_ExprSet <- new("ExpressionSet", exprs = as.matrix(Betas_combated), phenoData = annotated_pd_clean)
save(Betas_combated_ExprSet, file = beta.combat.expr.set.fn)

## Save normalized and batch-adjusted beta values
save(Betas_combated, file = beta.combat.fn)


#--- END















#-- 1. Combat Correction for plate

# Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
model.mtrx <- model.matrix(~1, data = pd_clean)

# Run ComBat to remove most significant batch effects itertively
M_combat_1plate <- ComBat(mval, batch = pd_clean$Sample_Plate, mod = model.mtrx)
save(M_combat_1plate, file = paste0(src.data.dir, "M_combat_1plate.Rdata"))

# Check to see if batch effect was succesfully removed
PCobj     <- prcomp(t(M_combat_1plate), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PC_Variation_by_batch_after_combated_plate.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA after Combat correction on plate")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA after Combat correction on plate")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA after Combat correction on plate")
dev.off()

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Sample_Plate) 
})
anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
anova.plate.tbl # PC1, PC5

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   295600   73900   1.092 0.3601
# Residuals              398 26935483   67677
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4   682842  170710  5.7962 0.0001528 ***
#   Residuals              398 11722011   29452
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   27432  6857.9  0.4006 0.8082
# Residuals              398 6813947 17120.5
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value  Pr(>F)
# Prin.comp$Sample_Plate   4  136412   34103  3.2126 0.01296 *
#   Residuals              398 4224910   10615
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   39502  9875.5  1.1084 0.3521
# Residuals              398 3545899  8909.3
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  202669   50667  6.2946 6.414e-05 ***
#   Residuals              398 3203607    8049
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
anova.slide.tbl # PC1, PC2

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  2302133   46043  0.6501  0.968
# Residuals                                352 24928950   70821
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50 2679984   53600  1.9401 0.0003237
# Residuals                                352 9724868   27627
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50 1588963   31779  2.1297 4.319e-05
# Residuals                                352 5252416   14922
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50 1660119   33202  4.3267 < 2.2e-16
# Residuals                                352 2701203    7674
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50 1084764 21695.3  3.0539 9.204e-10
# Residuals                                352 2500637  7104.1
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50 1089150 21783.0  3.3091 4.146e-11
# Residuals                                352 2317126  6582.7
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
anova.array.tbl # PC1, PC5

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   258415   36916  0.5406 0.8036
# Residuals       395 26972668   68285
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7 3081329  440190  18.649 < 2.2e-16 ***
#   Residuals       395 9323523   23604
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value   Pr(>F)
# Prin.comp$Array   7  345772   49396  3.0038 0.004384 **
#   Residuals       395 6495607   16445
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7  394804   56401  5.6166 3.432e-06 ***
#   Residuals       395 3966517   10042
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   53741  7677.3  0.8587 0.5394
# Residuals       395 3531661  8940.9
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value  Pr(>F)
# Prin.comp$Array   7  124876 17839.5  2.1474 0.03802 *
#   Residuals       395 3281400  8307.3
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-- 2. Combat Correction for Slide

mod             <- model.matrix(~1, data=pd_clean)
M_combat_2slide <- ComBat(M_combat_1plate,batch = pd_clean$Slide, mod = mod)
# save(M_combat_2slide, file = m.combat.2slide.fn)       

# Check to see if batch effect was succesfully removed
PCobj     <- prcomp(t(M_combat_2slide), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PC_Variation_by_batch_after_combated_plate_slide.pdf"))
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Sample_Plate), batch.legend.title = "Plate", title = "PCA after Combat correction on both Plate and Slide")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Slide), batch.legend.title = "Slide", legend.pos = 'right', title = "PCA after Combat correction on both Plate and Slide")
PlotPCADensity(PCobj, Prin.comp, batch = as.character(Prin.comp$Array), batch.legend.title = "Array", title = "PCA after Combat correction on both Plate and Slide")
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
anova.plate.tbl # PC1 - PC6
# 
# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   190371   47593  0.6426 0.6324
# Residuals              398 29478452   74066
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   126085   31521  1.0928 0.3597
# Residuals              398 11480362   28845
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   28067  7016.8  0.5114 0.7274
# Residuals              398 5460408 13719.6
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   37460  9365.1  1.1303 0.3417
# Residuals              398 3297469  8285.1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   46635 11658.8  1.8712 0.1147
# Residuals              398 2479740  6230.5
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   18310  4577.5   0.805 0.5225
# Residuals              398 2263092  5686.2

#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
anova.slide.tbl #P1-P6

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  1748986   34980   0.441 0.9997
# Residuals                                352 27919838   79318
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50   638073   12762  0.4095 0.9999
# Residuals                                352 10968373   31160
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  811630   16233  1.2217 0.1557
# Residuals                                352 4676846   13286
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  322855  6457.1  0.7546 0.8879
# Residuals                                352 3012074  8557.0
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  301277  6025.5  0.9532  0.567
# Residuals                                352 2225098  6321.3
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  322710  6454.2  1.1599 0.2239
# Residuals                                352 1958693  5564.5

#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
anova.array.tbl # PC1, PC6

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   356099   50871  0.6855 0.6844
# Residuals       395 29312724   74209
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7 3995441  570777  29.622 < 2.2e-16 ***
#   Residuals       395 7611006   19268
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value   Pr(>F)
# Prin.comp$Array   7  293114   41873  3.1836 0.002733 **
#   Residuals       395 5195361   13153
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7  677952   96850  14.398 < 2.2e-16 ***
#   Residuals       395 2656977    6727
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7  835697  119385  27.892 < 2.2e-16 ***
#   Residuals       395 1690678    4280
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   43910  6272.9  1.1074 0.3575
# Residuals       395 2237492  5664.5

#-- 3. Combat Correction for Array

mod             <- model.matrix(~1, data = pd_clean)
M_combat_3array <- ComBat(M_combat_2slide, batch = pd_clean$Array, mod = mod)
# save(M_combat_3array, file = m.combat.3array.fn)       

# Check to see if batch effect was succesfully removed
PCobj     <- prcomp(t(M_combat_3array), retx = T, center = T, scale. = T)
PCs       <- PCobj$x[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PC_Variation_by_batch_after_combated_plate_slide_array.pdf"))
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
anova.plate.tbl # PC1-PC6\{PC4}

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   200023   50006  0.6575 0.6219
# Residuals              398 30267925   76050
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4  101034   25258  1.3858 0.2381
# Residuals              398 7253980   18226
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   53234   13308   1.082  0.365
# Residuals              398 4895449   12300
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  124327 31081.6  5.2473 0.0003968 ***
#   Residuals              398 2357496  5923.4
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   21089  5272.3  0.9154 0.4548
# Residuals              398 2292232  5759.4
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4    6526  1631.5  0.3124 0.8697
# Residuals              398 2078659  5222.8

#-- for Slide
models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})
anova.slide.tbl <- sapply(models.slide, anova, simplify = F)
anova.slide.tbl #PC1-PC3, PC5

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  1820943   36419  0.4475 0.9996
# Residuals                                352 28647004   81384
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  716019   14320  0.7593 0.8827
# Residuals                                352 6638995   18861
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  704504   14090  1.1686 0.2132
# Residuals                                352 4244179   12057
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  576673 11533.5   2.131 4.263e-05
# Residuals                                352 1905149  5412.4
# 
# as.factor(as.character(Prin.comp$Slide)) ***
#   Residuals
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  346921  6938.4   1.242 0.1371
# Residuals                                352 1966400  5586.4
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value  Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  422920  8458.4  1.7911 0.00145 **
#   Residuals                                352 1662265  4722.3
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-- for Array
models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})
anova.array.tbl <- sapply(models.array, anova, simplify = F)
anova.array.tbl # PC1-PC5

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   150781   21540  0.2806 0.9614
# Residuals       395 30317166   76752
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7  185015   26431  1.4561 0.1815
# Residuals       395 7169999   18152
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   67762  9680.3  0.7834 0.6015
# Residuals       395 4880921 12356.8
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   63374  9053.5  1.4787 0.1732
# Residuals       395 2418448  6122.7
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   43828  6261.1  1.0897 0.3689
# Residuals       395 2269494  5745.6
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value  Pr(>F)
# Prin.comp$Array   7   79225 11317.9  2.2286 0.03124 *
#   Residuals       395 2005960  5078.4
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#-> if nothing significant ok, leave it like that.


##--- Convert the batch-adjusted M-values back into betas:

expit2         <- function(x) 2 ^ x / (1 + 2 ^ x)
Betas_combated <- expit2(M_combat_3array)
dim(Betas_combated) # for now final betas!
# 740357    403

#plot final densities
pdf(file = paste0(report.dir, "BetaValue_Distributions_afterNormCombat.pdf"))
densityPlot(Betas_combated, sampGroups = pd_clean$Slide, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Slide", xlab = "Beta")
densityPlot(Betas_combated, sampGroups = pd_clean$Array, legend = FALSE, main = "PostQC - Normalized and Batch Corrected Beta by Array", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated), rownames(pd_clean)) #TRUE

annotated_pd_clean     <- new("AnnotatedDataFrame", data= as.data.frame(pd_clean)) #extend to AnnotatedDataFrame (required for ESet)
Betas_combated_ExprSet <- new("ExpressionSet", exprs = as.matrix(Betas_combated), phenoData = annotated_pd_clean)
save(Betas_combated_ExprSet, file = beta.combat.expr.set.fn)

## Save normalized and batch-adjusted beta values
save(Betas_combated, file = beta.combat.fn)
