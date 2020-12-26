
src.dir      <- "/binder/mgp/workspace/2020_DexStim_Array_Human/methylation/"
src.data.dir <- "/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/" 

packages.fn <- paste0(src.dir, "packages.R")
rgSet.fn    <- paste0(src.data.dir, "rgSet_dex.RData")

source(packages.fn)

load(rgSet.fn)
rgSet <- rgSubSet

sh.script.mkdir <- sprintf("cd %s ; 
                    
                    if [ ! -d rData ] 
                    then 
                      mkdir rData 
                    fi ; 
                    
                    cd RData ;
      
                    ", src.data.dir)

system(sh.script.mkdir)


## Run in parallel

plan(multiprocess, workers = availableCores())

GetPhenoData <- function(rg.set){
  pd.orig <- pData(rg.set) # phenotype data
  save(pd.orig, file = "pd_original.Rdata")
}

GetRatioSet <- function(rg.set){
  Mset = preprocessRaw(rg.set)  # cpG locus level, with 2 channels methylated/ unmethylated
  save(Mset, file = "Mset_original.Rdata")

  RatioSet = ratioConvert(Mset, what = "both", keepCN = TRUE)# CpG locus level, but not mapped to a genome, Beta and M values
  save(RatioSet, file = "RatioSet_original.Rdata")

  return (RatioSet)

}


job.rg.set.1 %<-% GetPhenoData(rgSet)
job.rg.set.2 %<-% GetRatioSet(rgSet)

rg.set.list <- lapply(ls(pattern = "job.rg.set"), get)

ratio.set   <- rg.set.list[[2]] 


GetBetas <- function(ratio.set){
  RawBetas = getBeta(ratio.set) # Beta value matrix
  save(RawBetas,file = "RawBetas_original.Rdata")
}

GetgRatioSet <- function(ratio.set){
  gRatioSet =mapToGenome(ratio.set, mergeManifest=TRUE)
  save(gRatioSet, file = "gRatioSet_original.Rdata")
}

job.ratio.set.1 %<-% GetBetas(ratio.set)
job.ratio.set.2 %<-% GetgRatioSet(ratio.set)

ratio.set.list <- lapply(ls(pattern = "job.ratio.set"), get) 





