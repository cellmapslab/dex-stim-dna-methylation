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

report.dir <- paste0(src.dir, "05_batch_effect/reports/")

load(pd_clean.fn)
load(bmiq.quantileN.filtered.fn)
# load(quantileN.bmiq.fn)
# load(quantileN.fn)
# load(annotated_data_clean.fn)

##--- Combat to remove batch effects

#-- Function that will calculate the variance of each row
rowVars <- function(x, na.rm = FALSE, dims = 1, unbiased = TRUE, SumSquares = FALSE, twopass = FALSE) {
  if (SumSquares) return(rowSums(x^2, na.rm, dims))
  N <- rowSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N - 1 else N
  if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else
    sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
  (rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
}

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

PCobj <- prcomp(t(mval), retx = T, center = T, scale. = T)
save(PCobj, file = pcobj.fn)

load(pd_clean.fn)
load(pcobj.fn)
# pdf(paste0(report.dir, "boxplot_PCA.pdf"))
# boxplot(PCobj$x, col = "grey",frame=F)
# dev.off()

# Can use Scree plot to determine number of PCs to keep
pdf(paste0(report.dir, "screeplot_PCA.pdf"))
# plot(PCobj, type="line",cex.lab = 1.5, cex.main = 1.5) 
fviz_eig(PCobj, ncp = 20)
dev.off()

# Extract the PCs from the PCobj object
pcs.mtrx <- PCobj$x

# Extract the proportion of variability and cumulative proportion of 
# varibility explained by the top R PCs.
R <- 15
propvar <- round(summary(PCobj)$importance["Proportion of Variance", 1:R] * 100, 2 ) # -> write in xlsx
cummvar <- round(summary(PCobj)$importance["Cumulative Proportion", 1:R] * 100, 2) # -> write in xlsx
t(propvar); t(cumvar)

R <- 6 # choose eg 6 PCs
## Generate plots of the resulting PCs
# Plot of the proportion of variability explained by the top R PCs
# Plot of the cummulative proportion of variability explained by the top R PCs
# pdf(paste0(report.dir, "PCA_variance_explained.pdf"))
# par(mfrow = c(1,2))	
# par(mar = c(5, 5, 4, 2))
# barplot(propvar, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Variation Explained (%)", cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
# par(mar = c(5, 5, 4, 2))
# barplot(cummvar, xlab = paste("Top", R, "PCs", sep = " "), ylab = "Cumulative Variation Explained (%)",cex.axis = 1.5, cex.lab = 1.8, cex.names = 1.5)
# dev.off()

##---  Plot of pca individal map by Batch, group (dex, veh) and sex 

PCs <- PCobj$x
PCs <- PCs[, 1:R]
Prin.comp <- merge(PCs, pd_clean, by = "row.names", all = T) 

pdf(paste0(report.dir, "PC_Variation_by_batch.pdf"))

#-- by Plate
pca.plate <- fviz_pca_ind(PCobj,
             col.ind = as.factor(Prin.comp$Sample_Plate),
             geom = "point",
             repel = T)
ggpar(pca.plate,
      title = "PCA: DEX-Methylation data",
      subtitle = "by Plate",
      legend.title = "Plate", legend = "bottom")

#-- by Slide
pca.plate <- fviz_pca_ind(PCobj,
                          col.ind = as.factor(as.character(Prin.comp$Slide)),
                          geom = "point",
                          repel = T)
ggpar(pca.plate,
      title = "PCA: DEX-Methylation data",
      subtitle = "by Slide",
      legend = "none")

#-- by Array
pca.plate <- fviz_pca_ind(PCobj,
                          col.ind = as.factor(Prin.comp$Array),
                          geom = "point",
                          repel = T)
ggpar(pca.plate,
      title = "PCA: DEX-Methylation data",
      subtitle = "by Array",
      legend.title = "Array", legend = "bottom")

#-- by Group (dex, veh)
pca.plate <- fviz_pca_ind(PCobj,
                          col.ind = as.factor(Prin.comp$Sample_Group),
                          geom = "point",
                          repel = T)
ggpar(pca.plate,
      title = "PCA: DEX-Methylation data",
      subtitle = "by Group",
      legend.title = "Group", legend = "bottom")

#-- by sex (dex, veh)
pca.plate <- fviz_pca_ind(PCobj,
                          col.ind = as.factor(Prin.comp$sex),
                          geom = "point",
                          repel = T)
ggpar(pca.plate,
      title = "PCA: DEX-Methylation data",
      subtitle = "by Gender",
      legend.title = "Gender", legend = "bottom")

dev.off()


##--- Find extreme outliers

o1 <- 3 * sd(Prin.comp$PC1)
o2 <- 3 * sd(Prin.comp$PC2)
which(abs(Prin.comp$PC1) > o1 && abs(Prin.comp$PC2) > o2) # 0

##--- ANOVA for detection of variation between PCs and batch

#-- for Plate 

models.plate <- apply(PCs, 2, function(pc){
 lm(pc ~ Prin.comp$Sample_Plate) 
})

anova.plate.tbl <- sapply(models.plate, anova, simplify = F)
anova.plate.tbl # PC1, PC4
# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4   203341   50835  0.7748 0.5421
# Residuals              398 26113161   65611
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  1844092  461023  12.954 6.327e-10 ***
#   Residuals              398 14163944   35588
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  770089  192522  11.444 8.398e-09 ***
#   Residuals              398 6695598   16823
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Sample_Plate   4    7546  1886.6  0.1095 0.9792
# Residuals              398 6856884 17228.4
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  452102  113026  13.021 5.647e-10 ***
#   Residuals              398 3454663    8680
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Sample_Plate   4  197778   49445  5.8353 0.0001428 ***
#   Residuals              398 3372405    8473
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


#-- for Slide

models.slide <- apply(PCs, 2, function(pc){
  lm(pc ~ as.factor(as.character(Prin.comp$Slide)))
})

anova.slide.tbl <- sapply(models.slide, anova, simplify = F)

# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# as.factor(as.character(Prin.comp$Slide))  50  2366031   47321  0.6955 0.9414
# Residuals                                352 23950471   68041


#-- for Array

models.array <- apply(PCs, 2, function(pc){
  lm(pc ~ Prin.comp$Array) 
})

anova.array.tbl <- sapply(models.array, anova, simplify = F)
anova.array.tbl # PC1, PC5, PC6
# 
# $PC1
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7   286711   40959  0.6215 0.7382
# Residuals       395 26029791   65898
# 
# $PC2
# Analysis of Variance Table
# 
# Response: pc
# Df   Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7  2393835  341976   9.922 2.004e-11 ***
#   Residuals       395 13614200   34466
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC3
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value    Pr(>F)
# Prin.comp$Array   7  631113   90159  5.2107 1.068e-05 ***
#   Residuals       395 6834573   17303
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC4
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value   Pr(>F)
# Prin.comp$Array   7  391861   55980  3.4163 0.001473 **
#   Residuals       395 6472569   16386
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# $PC5
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7  102178 14596.8  1.5155 0.1603
# Residuals       395 3804588  9631.9
# 
# $PC6
# Analysis of Variance Table
# 
# Response: pc
# Df  Sum Sq Mean Sq F value Pr(>F)
# Prin.comp$Array   7  113967 16281.0  1.8607 0.0747 .
# Residuals       395 3456217  8749.9


##--- Batch correction
#first correct for plate
## Model matrix for batch-corrections (May need to adjust model matrix to 'protect' coefficients (study specific)):
model.mtrx <- model.matrix(~1, data = pd_clean)

## Run ComBat to remove most significant batch effects itertively
M_combat_1plate <- ComBat(mval, batch = pd_clean$Sample_Plate, mod = model.mtrx)
save(M_combat_1plate, file = paste0(src.data.dir, "M_combat_1plate.Rdata"))

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_1plate), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(paste0(report.dir, "PC_Variation_by_batch_afterCombat1.pdf"))
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array))
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 

#second correction for slide
mod <- model.matrix(~1, data=pd_clean)
M_combat_2slide = ComBat(M_combat_1plate,batch = pd_clean$Slide, mod = mod)
save(M_combat_2slide, file = m.combat.2slide.fn)       

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_2slide), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(paste0(report.dir, "PC_Variation_by_batch_afterCombat2.pdf"))
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:
# for plate
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide))
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array)) 
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 

#third correction for array
mod <- model.matrix(~1, data=pd_clean)
M_combat_3array = ComBat(M_combat_2slide,batch = pd_clean$Array, mod = mod)
save(M_combat_3array, file = m.combat.3array.fn)       

## Check to see if batch effect was succesfully removed
PCobj = prcomp(t(M_combat_3array), retx = T, center = T, scale. = T)
PCs = PCobj$x
PCs =PCs[,1:R]
Prin.comp<-merge(PCs,pd_clean, by = "row.names",all=T) 

#### Check whether batches are still distinguished by first and second PC:
pdf(paste0(report.dir, "PC_Variation_by_batch_afterCombat3.pdf"))
par(mfrow=c(3,1))
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Sample_Plate), xlab = "PC1", ylab = "PC2", main="Plate")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Slide), xlab = "PC1", ylab = "PC2", main="Slide")
plot(Prin.comp$PC1,Prin.comp$PC2,pch=16, col=as.factor(Prin.comp$Array), xlab = "PC1", ylab = "PC2", main="Array")
dev.off()

## Can further test via ANOVA, to look for strong variations in PCs by Batch:

# for plate
anova(lm(Prin.comp$PC1~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC2~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC3~Prin.comp$Sample_Plate)) 
anova(lm(Prin.comp$PC4~Prin.comp$Sample_Plate)) 

#for slide
anova(lm(Prin.comp$PC1~Prin.comp$Slide))
anova(lm(Prin.comp$PC2~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC3~Prin.comp$Slide)) 
anova(lm(Prin.comp$PC4~Prin.comp$Slide)) 

#for array
anova(lm(Prin.comp$PC1~Prin.comp$Array)) 
anova(lm(Prin.comp$PC2~Prin.comp$Array)) 
anova(lm(Prin.comp$PC3~Prin.comp$Array)) 
anova(lm(Prin.comp$PC4~Prin.comp$Array)) 


#-> if nothing significant ok, leave it like that.


#################################################################################
## Convert the batch-adjusted M-values back into betas:

expit2 = function(x) 2^x/(1+2^x)
Betas_combated = expit2(M_combat_3array)
dim(Betas_combated) # for now final betas!
# 716331    264

## Save normalized and batch-adjusted beta values
save(Betas_combated, file = beta.combat.fn)

#plot final densities
png(paste0(report.dir, "BetaValue_Distributions_afterNormCombat.png",width=700,height=700,pointsize=12))
densityPlot(Betas_combated, sampGroups = pd_clean$Slide, legend=FALSE, main = "PostQC - Normalized and Batch Corrected Beta", xlab = "Beta")
dev.off() 

all.equal(colnames(Betas_combated),rownames(pd_clean)) #TRUE
annotated_pd_clean =new("AnnotatedDataFrame", data= as.data.frame(pd_clean)) #extend to AnnotatedDataFrame (required for ESet)

Betas_combated_ExprSet = new("ExpressionSet", exprs = as.matrix(Betas_combated), phenoData=annotated_pd_clean)
save(Betas_combated_ExprSet, file = beta.combat.expr.set.fn)
