# 7. Check if there are any mixups

# 1. MixupMapper Genotype data

mixup.res.fn   <- "/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/out_mixupmapper/MixupMapper/BestMatchPerGenotype.txt"
mixup.res.fn   <- "/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/out_mixupmapper_removed_mixups/MixupMapper/BestMatchPerGenotype.txt"
mixup.genotype <- read.table(mixup.res.fn, stringsAsFactors = FALSE)

colnames(mixup.genotype) <- mixup.genotype[1, ]
mixup.genotype           <- mixup.genotype[-1,]

table(mixup.genotype$Mixup)
# false  true
# 199     2

mixup.genotype[mixup.genotype$Mixup == "true", ]          

# Genotype        OriginalLinkedTrait      OriginalLinkedTraitScore BestMatchingTrait    BestMatchingTraitScore Mixup
# MPIPSYKL_007875  200712160042_R01C01      -4.9666422637773895     200720060022_R04C01  -12.023177359249031  true
# MPIPSYKL_007893  200720060022_R04C01      -5.231460612404898      200712160042_R01C01   -6.322773218809481  true

# Sample_ID	Sample_Plate	Sample_Well	project	person	Basename	Array	Slide	Sample_Name	Source_Dir	status	sex	Sample_Group	age	bmi
# R08267	EPIC_Juni2016_06	D2	DEX	MPIPSYKL_007893	200720060022_R04C01		1	W	veh	66	212.671.029.149.316
# R08343	EPIC_Juni2016_07	A4	DEX	MPIPSYKL_007875	200712160042_R01C01		1	W	veh	49	25

idividuals.to.exclude <- mixup.genotype[mixup.genotype$Mixup == "true", "Genotype"] 

# 2. Exclude samples

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)


for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# MixupsMapper was run only for VEH. However, we should also exclude the same individuals from DEX

x     <- load(pd_clean.fn)
pheno <- get(x)

samples.id.to.exclude <- pheno[pheno$person %in% idividuals.to.exclude, "Sample_Name"]

# Load data from which samples should be excluded
x                     <- load(beta.combat.fn)
beta.combat.mtrx      <- get(x)
x                     <- load(bmiq.quantileN.fn)
bmiq.quantileN        <- get(x)
x                     <- load(bmiq.quantileN.filtered.fn)
bmiq.quantileN.filter <- get(x)
x                     <- load(beta.combat.expr.set.fn)
beta.combat.expr.set  <- get(x) 
x                     <- load(rgSet_clean.fn)
rgset.clean           <- get(x)
x                     <- load(detP_clean.fn)
detP                  <- get(x)

rm(x)
# Filter out mix-ups

pheno.mixups                 <- pheno[!(pheno$Sample_Name %in% samples.id.to.exclude), ] # 399 x 16
detP.mixups                  <- detP[, !(colnames(detP) %in% samples.id.to.exclude)] # 866836 x 399
bmiq.quantileN.mixups        <- bmiq.quantileN[, !(colnames(bmiq.quantileN) %in% samples.id.to.exclude)] #  865859 x 399
bmiq.quantileN.filter.mixups <- bmiq.quantileN.filter[, !(colnames(bmiq.quantileN.filter) %in% samples.id.to.exclude)] # 740357 x 399
beta.combat.mtrx.mixups      <- beta.combat.mtrx[, !(colnames(beta.combat.mtrx) %in% samples.id.to.exclude)] # 740357 x 399
beta.combat.expr.set.mixups  <- beta.combat.expr.set[ , !(colnames(beta.combat.expr.set) %in% samples.id.to.exclude)] # 740357 x 399
rgset.clean.mixups           <- rgset.clean[ , !(colnames(rgset.clean) %in% samples.id.to.exclude)] # 1052641 x 399

# Save cleaned data

src.final.data.dir<- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/10_final_qc_data/"

saveRDS(pheno.mixups, paste0(src.final.data.dir, "dex_methyl_phenotype.rds")) # final phenotype data after excluding poor qc and mix-ups samples
saveRDS(detP.mixups, paste0(src.final.data.dir, "dex_methyl_detP.rds")) # final p-values matrix after excluding poor qc and mix-ups samples but not probes
saveRDS(bmiq.quantileN.mixups, paste0(src.final.data.dir, "dex_methyl_bmiq_quantileN.rds")) #final beta matrix after quantile + BMIQ normalization, samples filtering and mix-ups removal
saveRDS(bmiq.quantileN.filter.mixups, paste0(src.final.data.dir, "dex_methyl_bmiq_quantileN_filtered.rds")) # final beta matrix after normalization, probes and samples filtering and mix-ups removal
saveRDS(beta.combat.mtrx.mixups, paste0(src.final.data.dir, "dex_methyl_beta_combat_mtrx.rds")) # final beta matrix after normalization, probes and samples filtering, batch correction and mix-ups removal
saveRDS(beta.combat.expr.set.mixups, paste0(src.final.data.dir, "dex_methyl_beta_combat_exprset.rds")) # final beta expression set after normalization, probes and samples filtering, batch correction and mix-ups removal
saveRDS(rgset.clean.mixups, paste0(src.final.data.dir, "dex_methyl_rgset_final.rds"))  # final RGChannel Set after removing poor qc samples and mix-ups, number of probes is the same as in initial 

