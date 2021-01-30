args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# print(rgSet_clean.fn)

source(packages.fn)
source(bmiq.script.fn)

load(quantileN.filtered.fn)

Betas_quantileN_filtered <- getBeta(quantileN_filtered) # beta values
save(Betas_quantileN_filtered, file = betas.quantileN.filtered.fn)

