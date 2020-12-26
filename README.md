# __DEX-stimulated DNAm arrays processing__

**Input data:**

- iData data : `/binder/mgp/workspace/2020_DexStim_Array_Human/methylation/`
- RData: `/binder/mgp/datasets/2020_DexStim_Array_Human/methylation/`

    * `rgSet_dex`: raw data from the dex IDAT files; organized at the probe (not CpG locus) level and has two channels (Red and Green).
    
    The following data are obtained from __rgSet_dex__:
    
    * `Mset_original`: data organized by the CpG locus level, but not mapped to a genome and has two channels, _Meth (methylated)_ and _Unmeth (unmethylated)_;
    * `RatioSet_original`: data organized by the CpG locus level, but not mapped to a genome, and has at least one of two channels, _Beta_ and/or _M (logratio of Beta)_;
    * `RawBetas`: _Beta_ values matrix, obtained from _RatioSet_;
    * `gRatioSet_original`: the same as _RatioSet_ but mapped to a genome;
    * `pd_original`: phenotype data

