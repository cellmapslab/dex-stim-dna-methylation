# 7. Check if there are any mixups

# 1. MixupMapper Genotype data

mixup.res.fn   <- "/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/out_mixupmapper/MixupMapper/BestMatchPerGenotype.txt"
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

# 2. Exclude samples

input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)


for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# Load data from which samples should be excluded
x                     <- load(beta.combat.fn)
beta.combat.mtrx      <- get(x)

x                     <- load(bmiq.quantileN.filtered.fn)
bmiq.quantileN.filter <- get(x)

x                     <- load(beta.combat.expr.set.fn)
beta.combat.expr.set  <- get(x) 

x     <- load(pd_clean.fn)
pheno <- get(x)

200720060022_R03C01