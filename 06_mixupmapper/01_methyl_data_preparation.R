
args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])
# or
input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# source(packages.fn)
# source("functions.R")

report.dir           <- paste0(src.dir, "06_mixupmapper/01_reports/")
dnam.mixupmapper.dir <- "/home/ahryhorzhevska/mpip/datasets/methylation/mixupmapper/"

##--- Genotype - phenotype coupling

load(pd_clean.fn) # pd_clean
# should be only 2 columns: Sample_Name & ArrayID
genomeID           <- data.frame(pd_clean$person, pd_clean$Sample_Name)[pd_clean$Sample_Group == "veh", ]
colnames(genomeID) <- c("IndividualID", "ArrayID")

genotypemethylationcoupling <- genomeID
write.table(genotypemethylationcoupling, 
            file = paste0(dnam.mixupmapper.dir, "genotypemethylationcoupling.txt"), 
            row.names = F, col.names = F, quote = F, sep = "\t")

##--- batch-adjusted beta values to txt format 

# Only veh samples, because genotype are only veh

load(beta.combat.fn) # Beta_combated
load(pd_clean.fn)

betas.colnames <- colnames(Betas_combated)
samples.veh.id <- pd_clean$Sample_Name[pd_clean$Sample_Group == "veh"]

betas.combated.veh <- Betas_combated[ , samples.veh.id ]

write.table(betas.combated.veh, 
            file = paste0(dnam.mixupmapper.dir, "betas_combat_veh_mixupmapper.txt"), 
            sep = "\t", quote = F, row.names = T, col.names = T)

sh.script <- sprintf("BETAS_FILENAME='betas_combat_veh_mixupmapper.txt'; \
                      sce' ; \
                      FIRST_SAMPLE=$(cat $BETAS_FILENAME | head -n1 | awk '{print $1;}') ;\
                      sed 's/$FIRST_SAMPLE/\t$FIRST_SAMPLE/' $BETAS_FILENAME > $BETAS_FILENAME_final.txt")
system(sh.script)
