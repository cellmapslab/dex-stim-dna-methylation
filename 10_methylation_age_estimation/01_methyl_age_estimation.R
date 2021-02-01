# 10 Methylation age estimation

# 1. Set up env and load data

library(minfi)
library(ENmix)
library(SummarizedExperiment)

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)


for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

x         <- load(beta.combat.fn)
beta.mtrx <- get(x)
rm(x)

# OR
# beta.mtrx <- readR(paste0(src.data.dir, "betas_mtrx_after_gap_extreme_outliers_na.rds"))

meth.age      <- methyAge(beta.mtrx, 
                          type = "all", 
                          # fastImputaion = F,
                          normalize = F)

write.csv(meth.age, file = paste0(report.dir, "DEX-stim-array-human-meth-age.csv"), col.names = T, row.names = T, quote = F, sep = ";")
