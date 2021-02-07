# 7. Prepare methylation data fro MixupMapper

# 1. Load data
input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

dnam.mixupmapper.dir <- "/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/"

#2. Genotype - phenotype coupling

# pd_clean.fn <- paste0(src.final.data.dir, "dex_methyl_phenotype.rds")
# pd_clean    <- readRDS(pd_clean.fn)
load(pd_clean.fn) # pd_clean
# should be only 2 columns: Sample_Name & ArrayID
genomeID           <- data.frame(pd_clean$person, pd_clean$Sample_Name)[pd_clean$Sample_Group == "veh", ]
colnames(genomeID) <- c("IndividualID", "ArrayID")

genotypemethylationcoupling <- genomeID
write.table(genotypemethylationcoupling, 
            file = paste0(dnam.mixupmapper.dir, "genotypemethylationcoupling_removed_mixups.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

# 3. Batch-adjusted beta values to txt format 

# Only veh samples, because genotype are only veh

# beta.mtrx.fn   <- paste0(src.final.data.dir, "dex_methyl_beta_combat_mtrx.rds")
# Betas_combated <- readRDS(beta.mtrx.fn)
load(beta.combat.fn) # Beta_combated
load(pd_clean.fn)

betas.colnames <- colnames(Betas_combated)
samples.veh.id <- pd_clean$Sample_Name[pd_clean$Sample_Group == "veh"]

betas.combated.veh <- Betas_combated[ , samples.veh.id ]

# Format for methylation data for MixupMapper: IDs in columns, probes in rows, first column-name is an empty tab!
first.col.name <- colnames(betas.combated.veh)[1]
colnames(betas.combated.veh)[1] <- paste0('\t', first.col.name)

write.table(betas.combated.veh, 
            file = paste0(dnam.mixupmapper.dir, "betas_combat_veh_mixupmapper_removed_mixups.txt"), 
            sep = "\t", quote = F, row.names = T, col.names = T)

# ALternative way to adjusting first column name in the shell (time consuming)
sh.script <- sprintf("BETAS_FILENAME='betas_combat_veh_mixupmapper.txt'; \
                      sce' ; \
                      FIRST_SAMPLE=$(cat $BETAS_FILENAME | head -n1 | awk '{print $1;}') ;\
                      sed 's/$FIRST_SAMPLE/\t$FIRST_SAMPLE/' $BETAS_FILENAME > $BETAS_FILENAME_final.txt")
system(sh.script)
