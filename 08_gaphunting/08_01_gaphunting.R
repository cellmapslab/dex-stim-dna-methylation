# 07 Batch effects correction

# 1. Load data

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

# OR
input.parameters.fn <- "input_parameters.csv"

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

library(minfi)
# source("05_batch_effects_correction/05_00_functions.R")

report.dir <- paste0(src.dir, "08_gaphunting/01_reports/")

x <- load(beta.combat.fn)
betas.combated <- get(x)
rm(x) 

# 2. Set up input parameters

betas.mtrx <- betas.combated
n.samples  <- ncol(betas.mtrx)
threshold  <- 0.1 #c(0.025, 0.05, 0.1)
cutoff     <- ifelse(5 > (0.025 * n.samples), 5/n.samples ,0.0025)

# Check if there is any missing value
# Gaphunter function cannot handle missing value
table(is.na(betas.mtrx)) # none
                     
# 3. Runnning gaphunter

# Keep outlier-driven gap signals in the results
gap.res.with.outliers     <- gaphunter(betas.mtrx, keepOutliers = TRUE, threshold = threshold, outCutoff = cutoff) # note result in the excel

signal.probes.all  <- gap.res.with.outliers$proberesults # signals by groups (Groups, Group_1, .... Group_k) 
signal.samples.all <- gap.res.with.outliers$sampleresults # signals by individuals

# Combine signal probes and samples into one table
# Merge by probe name
signals.all        <- merge(signal.probes.all, signal.samples.all, by.x = "row.names", by.y = "row.names")

# Gaphunter result without outliers
gap.res.without.outliers  <- gaphunter(betas.mtrx, keepOutliers = FALSE, threshold = threshold, outCutoff = cutoff) 
signal.probes.without.outliers <- gap.res.without.outliers$proberesults

# Get outlier
outliers <- setdiff(rownames(signal.probes.all), rownames(signal.probes.without.outliers))

all.outlier.signals           <- signals.all[which(signals.all[, 1] %in% outliers),] 
no.groups                     <- max(all.outlier.signals$Groups)
start.sample.pos              <- ncol(all.outlier.signals) - n.samples + 1
all.outlier.signals           <- all.outlier.signals[, c(1:(no.groups + 2), start.sample.pos:ncol(all.outlier.signals))]
colnames(all.outlier.signals) <- c("Probes", colnames(all.outlier.signals)[2:(no.groups + 2)], colnames(betas.mtrx))

# Replace outliers in beta values matrix with NA
betas.mtrx.outlier.na <- NULL
for(p in 1:nrow(all.outlier.signals)){
  group.number   <- as.numeric(strsplit(names(which.max(all.outlier.signals[p, 3:(no.groups + 2)])),"")[[1]][6])
  df             <- all.outlier.signals[p, c(rep(TRUE, (no.groups + 2)), all.outlier.signals[p, -c(1:(no.groups + 2))] != group.number)]
  outlier.sample <- colnames(df)[-c(1:(no.groups + 2))]
  betas.mtrx.tmp <- betas.mtrx[which(rownames(betas.mtrx) %in% df$Probes),]
  betas.mtrx.tmp[names(betas.mtrx.tmp) %in% outlier.sample] <- NA
  
  betas.mtrx.outlier.na <- rbind(betas.mtrx.outlier.na, c(Probes = df$Probes, betas.mtrx.tmp))
}

rownames(betas.mtrx.outlier.na) <- betas.mtrx.outlier.na[, 1]
betas.mtrx.tmp                  <- data.frame(betas.mtrx.outlier.na[, -1], check.names = FALSE)
betas.mtrx.gapped.outliers.na   <- as.data.frame(betas.mtrx)

# Final betas matrix (saved as betas_after_gap in .Rda format) with all outliers detected as NAs 
betas.mtrx.gapped.outliers.na [match(rownames(betas.mtrx.tmp) ,rownames(betas.mtrx.gapped.outliers.na )), ] <- betas.mtrx.tmp
dim(betas.mtrx.gapped.outliers.na)
saveRDS(betas.mtrx.gapped.outliers.na, file = paste0(src.data.dir, "betas_mtrx_gapped_outliers_na.rds"))

# Beta matrix from original betas combated mtrx with additional NAs for extreme outliers.
betas.mtrx[is.na(betas.mtrx.tmp)] <- NA; 
betas.mtrx.gapped.extreme.outliers.na <- betas.mtrx
saveRDS(betas.mtrx.gapped.extreme.outliers.na, file = paste0(src.data.dir, "betas_mtrx_gapped_extreme_outliers_na.rds")) ## Use this for further analysis

dim(betas.mtrx.gapped.extreme.outliers.na) 

# Create summary files for probes 
probes.na.before.gap  <- data.frame(Probes = rownames(betas.mtrx), 
                                    Nr_NA_before_gap = apply(betas.mtrx, 1, function(x) sum(is.na(x))))
probes.na.after.gap <- data.frame(Probes = rownames(betas.mtrx.gapped.outliers.na), 
                                  Nr_NA_after_gap = apply(betas.mtrx.gapped.outliers.na, 1, function(x) sum(is.na(x))))

robes.na.before.after <- merge(probes.na.before.gap, probes.na.after.gap, by.x = "Probes", by.y = "Probes", all.x = TRUE, all.y = TRUE)
probes.na.before.after$outlier_gaphunter <- probes.na.before.after$Nr_NA_after_gap - probes.na.before.after$Nr_NA_before_gap
probes.na.before.after$percent_outlier_gaphunter <- (probes.na.before.after$outlier_gaphunter/(ncol(betas.mtrx.gapped.outliers.na) - 
                                                                                                 probes.na.before.after$Nr_NA_before_gap)) * 100
write.csv(probes.na.before.after , paste0(report.dir, "01_Gaphunter_CpGs_Summary_threshold_", threshold, ".csv"), quote = F, row.names = F, sep = "\t")

# Create summary files for samples 
samples.na.before.gap <- data.frame(Samples = colnames(betas.mtrx), 
                                    Nr_NA_before_gap = apply(betas.mtrx, 2, function(x) sum(is.na(x))))

samples.na.after.gap <- data.frame(Samples = colnames(betas.mtrx.gapped.outliers.na),
                                   Nr_NA_after_gap = apply(betas.mtrx.gapped.outliers.na, 2, function(x) sum(is.na(x))))
samples.na.before.after <- merge(samples.na.before.gap, samples.na.after.gap, by.x = "Samples", by.y = "Samples", all.x = TRUE, all.y = TRUE)
samples.na.before.after$outlier_gaphunter <- samples.na.before.after$Nr_NA_after_gap - samples.na.before.after$Nr_NA_before_gap
samples.na.before.after$percent_outlier_gaphunter <- (samples.na.before.after$outlier_gaphunter / (nrow(betas.mtrx.gapped.outliers.na) -
                                                                                                   samples.na.before.after$Nr_NA_before_gap)) * 100
write.csv(samples.na.before.after, paste0(report.dir, "02_Gaphunter_Samples_Summary_threshold_", threshold, ".csv"))  

