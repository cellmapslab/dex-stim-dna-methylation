if(!require(BiocManager)) {
install.packages("BiocManager"); require(BiocManager)}

if(!require(minfi)) {
  BiocManager::install("minfi"); require(minfi)}

if(!require(future)) {
  install.packages("future"); require(future)}

if(!require(RColorBrewer)) {
  install.packages("RColorBrewer"); require(RColorBrewer)}

