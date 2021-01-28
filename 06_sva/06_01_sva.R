# 06 Surrogate Variable Analysis (SVA)

# Reference: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf

library(sva)

input.parameters.fn <- "input_parameters.csv"
input.parameters    <- read.csv(input.parameters.fn, sep = ";", header = F )
input.parameters    <- as.data.frame(input.parameters)

for (row in 1:nrow(input.parameters))
  assign(input.parameters[row, 1], input.parameters[row, 2])

report.dir <- paste0(src.dir, "06_sva/01_reports/")

# 1. Setting up the data from an ExpressionSet
x     <- load(pd_clean.fn)
pheno <- get(x)

x        <- load(beta.combat.expr.set.fn) #Betas_combated_ExprSet
expr.set <- get(x)

rm(x)

# 2. Create the full and null model matrices
mod  <- model.matrix(~Sample_Group, data = pheno)
mod0 <- model.matrix(~ 1, data = pheno))

# 3. Applying the 'sva' to estimate surrogate variables
n.sv  <- num.sv(expr.set, mod, method="leek")
svobj <- sva(expr.set, mod, mod0, n.sv = n.sv)


