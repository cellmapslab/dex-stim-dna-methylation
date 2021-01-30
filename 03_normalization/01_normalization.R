# 03_normalisation

args                <- commandArgs(T)
input.parameters.fn <- as.character(args[1])

input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

# print(rgSet_clean.fn)
source(packages.fn)

load(rgSet_clean.fn)
# load(detP_clean.fn)
# load(rawBetas_clean.fn)
# load(pd_clean.fn)
# load(annotated_data_clean.fn)

quantileN <- preprocessQuantile(rgSet_clean)
save(quantileN, file = quantileN.fn)