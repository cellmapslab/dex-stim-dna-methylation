library("minfi")

src.dir <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/"

sample.sheet.fn <- paste0(src.dir, "samplesheet_epic_methylation_dex.csv")
rgSet.fn        <- paste0(src.dir, "RGSet_original.Rdata")

targets         <- read.csv(sample.sheet.fn, sep = ',' )

load(rgSet.fn)

rgSet           <- RGSet
rm(RGSet)

rgSubSet        <- rgSet[, targets$SampleName]

save(rgSubSet, file = paste0(src.dir, "rgSet_dex.RData"))

