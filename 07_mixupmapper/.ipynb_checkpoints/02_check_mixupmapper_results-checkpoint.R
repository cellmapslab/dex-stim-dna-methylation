# 7. Check if there are any mixups

# 1. MixupMapper Genotype data

mixup.genotype <- read.table("/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/out_mixupmapper/MixupMapper/BestMatchPerGenotype.txt")

colnames(mixup.genotype) <- mixup.genotype[1, ]
mixup.genotype           <- mixup.genotype[-1,]

table(mixup.genotype$Mixup)
# false  true
# 199     2

mixup.genotype[mixup.genotype$Mixup == "true", ]          

# Genotype        OriginalLinkedTrait      OriginalLinkedTraitScore BestMatchingTrait    BestMatchingTraitScore Mixup
# MPIPSYKL_007875  200712160042_R01C01      -4.9666422637773895     200720060022_R04C01  -12.023177359249031  true
# MPIPSYKL_007893  200720060022_R04C01      -5.231460612404898      200712160042_R01C01   -6.322773218809481  true

#Sample_ID	Sample_Plate	Sample_Well	project	person	Basename	Array	Slide	Sample_Name	Source_Dir	status	sex	Sample_Group	age	bmi
R08267	EPIC_Juni2016_06	D2	DEX	MPIPSYKL_007893	200720060022_R04C01		1	W	veh	66	212.671.029.149.316
R08343	EPIC_Juni2016_07	A4	DEX	MPIPSYKL_007875	200712160042_R01C01		1	W	veh	49	25

