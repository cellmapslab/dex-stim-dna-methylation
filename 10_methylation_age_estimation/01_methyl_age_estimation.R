# 10 Methylation age estimation

# 1. Set up env and load data

library(minfi)
library(ENmix)
library(SummarizedExperiment)
library(ggplot2)
library(tidyverse)

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)


for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

beta.mtrx.fn <- paste0(src.final.data.dir, "dex_methyl_beta_combat_mtrx.rds")
beta.mtrx    <- readRDS(beta.mtrx.fn)

pheno.fn     <- paste0(src.final.data.dir, "dex_methyl_phenotype.rds")
pheno        <- readRDS(pheno.fn)

# 2. Estimate methylation age

meth.age      <- methyAge(beta.mtrx, 
                          type = "all", 
                          fastImputation = F,
                          normalize = F)

# 3. Merge estimated age with real and write out results

real.age.tbl   <- data.frame(SampleID = pheno$Sample_Name, realAge = pheno$age, Individual = pheno$person )
merged.age.tbl <- merge(meth.age, real.age.tbl)

write.table(merged.age.tbl, file = paste0(report.dir, "DEX-stim-array-human-meth-age.csv"), col.names = T, row.names = T, quote = F, sep = ";")

# 4. Relate DNAm age to chronological age

MedianAbsDevHelp <- function(x, y) median(abs(x-y), na.rm=TRUE)
MedianAbsDev     <- function(dnam.age.method, real.age) signif(MedianAbsDevHelp(dnam.age.method, real.age), 2)


pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_Relation.pdf"))

err <- MedianAbsDev (merged.age.tbl$mAge_Hovath, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = mAge_Hovath, y = realAge,  color = pheno$Sample_Group)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hovath) vs Chronological age, err = ", err),
       x = "DNAm Age Hovath", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl$mAge_Hannum, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = mAge_Hannum, y = realAge, color = pheno$Sample_Group)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hannum) vs Chronological age, err = ", err),
       x = "DNAm Age Hannum", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl$PhenoAge, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = PhenoAge, y = realAge, color = pheno$Sample_Group)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Pheno) vs Chronological age, err = ", err),
       x = "DNAm Pheno Age ", y = "Chronological Age")

dev.off()

# 5. DNAm age dependence on DEX

group.samples.tbl <- data.frame(SampleID = pheno$Sample_Name, Group = pheno$Sample_Group)
merged.age.tbl    <- merge(merged.age.tbl, group.samples.tbl)

real.age.tbl  <- merged.age.tbl %>% select(SampleID, Age = realAge, Group) %>% mutate(Type = "ChronAge") 
pheno.age.tbl <- merged.age.tbl %>% select(SampleID, Age = PhenoAge, Group) %>% mutate(Type = "PhenoAge")

pheno.real.age.tbl <- rbind(real.age.tbl, pheno.age.tbl)

pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_DEX_Relation.pdf"))
ggplot(data = pheno.real.age.tbl, aes(x = Type, y = Age, fill = Group)) +
  geom_bar(stat = "identity", fill = "steelblue")+
  geom_text(aes(label = Group), vjust = 1.6, color = "white", size = 3.5)+
  theme_minimal()
dev.off()


pal <- brewer.pal(8, "Set2")
pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_DEX_Relation.pdf"))
age.tbl.dex <- merged.age.tbl[pheno$Sample_Group == "dex", -1]
age.tbl.veh <- merged.age.tbl[pheno$Sample_Group == "veh", -1]
boxplot(age.tbl.dex, xlim = c(0, 12), ylim = range(age.tbl.dex, age.tbl.veh), 
        at = 0:3 * 3 + 0.5, xaxt = "n", main = "", ylab = "Age", col = pal[5])
boxplot(age.tbl.veh, at = 0:3 * 3 + 1.5, xaxt = "n", col = pal[6], add = TRUE)
axis(1, at = 0:3 * 3 + 1.5, labels = colnames(age.tbl.dex), tick = TRUE)
legend("topleft", legend = c("dex", "veh"), fill = pal[5:6])
dev.off()

# 6.  Estimate methylation age only for veh

betas.colnames <- colnames(beta.mtrx)
samples.veh.id <- pheno$Sample_Name[pheno$Sample_Group == "veh"]

betas.mtrx.veh <- beta.mtrx[ , samples.veh.id ]

meth.age.veh   <- methyAge(betas.mtrx.veh, 
                          type = "all", 
                          fastImputation = F,
                          normalize = F)

merged.age.tbl.veh <- merge(meth.age.veh, real.age.tbl)

write.table(merged.age.tbl.veh, file = "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/12_DNAm_age/DEX-stim-array-human-meth-age-veh.csv", col.names = T, row.names = F, quote = F, sep = ";")
write.csv2(merged.age.tbl.veh, file = paste0(report.dir, "/DEX-stim-array-human-meth-age-veh.csv"), col.names = T, row.names = F, quote = F, sep = ";")

# Relate DNAm age to chronological age

MedianAbsDevHelp <- function(x, y) median(abs(x-y), na.rm=TRUE)
MedianAbsDev     <- function(dnam.age.method, real.age) signif(MedianAbsDevHelp(dnam.age.method, real.age), 2)


pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_Relation_VEH.pdf"))

err <- MedianAbsDev (merged.age.tbl.veh$mAge_Hovath, merged.age.tbl.veh$realAge)
ggplot(merged.age.tbl.veh, aes(x = mAge_Hovath, y = realAge)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hovath) vs Chronological age, err = ", err),
       x = "DNAm Age Hovath", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl.veh$mAge_Hannum, merged.age.tbl.veh$realAge)
ggplot(merged.age.tbl.veh, aes(x = mAge_Hannum, y = realAge)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hannum) vs Chronological age, err = ", err),
       x = "DNAm Age Hannum", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl.veh$PhenoAge, merged.age.tbl.veh$realAge)
ggplot(merged.age.tbl.veh, aes(x = PhenoAge, y = realAge)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Pheno) vs Chronological age, err = ", err),
       x = "DNAm Pheno Age ", y = "Chronological Age")

dev.off()
