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

x         <- load(beta.combat.fn)
beta.mtrx <- get(x)

x     <- load(pd_clean.fn)
pheno <- get(x)

rm(x)

# 2. Estimate methylation age

meth.age      <- methyAge(beta.mtrx, 
                          type = "all", 
                           fastImputation = F,
                          normalize = F)

# 3. Merge estimated age with real and write out results

real.age.tbl   <- data.frame(SampleID = pheno$Sample_Name, realAge = pheno$age)
merged.age.tbl <- merge(meth.age, real.age.tbl)

write.table(merged.age.tbl, file = paste0(report.dir, "DEX-stim-array-human-meth-age.csv"), col.names = T, row.names = T, quote = F, sep = ";")

# 4. Relate DNAm age to chronological age

MedianAbsDevHelp <- function(x, y) median(abs(x-y), na.rm=TRUE)
MedianAbsDev     <- function(dnam.age.method, real.age) signif(MedianAbsDevHelp(dnam.age.method, real.age), 2)


pdf(file = paste0(report.dir, "DNAm_Age_and_Chronological_Age_Relation.pdf"))

err <- MedianAbsDev (merged.age.tbl$mAge_Hovath, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = mAge_Hovath, y = realAge)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hovath) vs Chronological age, err = ", err),
       x = "DNAm Age Hovath", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl$mAge_Hannum, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = mAge_Hannum, y = realAge)) + 
  geom_point() +
  geom_smooth() +
  labs(title = paste("DNAm age (Hannum) vs Chronological age, err = ", err),
       x = "DNAm Age Hannum", y = "Chronological Age")

err <- MedianAbsDev (merged.age.tbl$PhenoAge, merged.age.tbl$realAge)
ggplot(merged.age.tbl, aes(x = PhenoAge, y = realAge)) + 
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