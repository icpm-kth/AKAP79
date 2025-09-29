#!/usr/bin/env Rscript

## turn conservatuion law analysis on:
# ./makeModel.R --cla
# ./makeModel.R -c
## off by default

a <- commandArgs(trailingOnly=TRUE)
if (any(grepl("--?c(la)?",a))){
	conservationLawAnalysis <- TRUE
} else {
	conservationLawAnalysis <- FALSE
}
files <- dir(".",pattern="tsv$")
sb <- SBtabVFGEN::sbtab_from_tsv(files)
## make the model variable
m <- SBtabVFGEN::sbtab_to_vfgen(sb,cla=conservationLawAnalysis)
saveRDS(m,"mcmc/AKAP79.RDS")
message("model variable saved to [mcmc/AKAP79.RDS]")
ex <- SBtabVFGEN::sbtab.data(sb,m$conservationLaws)
saveRDS(ex,"mcmc/experiments.RDS")
message("model experiments saved to [mcmc/experiments.RDS]")

C <- uqsa::generateCode(m)
cat(C,sep="\n",file="C/AKAP79_gvf.c")
message("C code written to [C/AKAP79_gvf.c]")

R <- uqsa::generateRCode(m)
cat(R,sep="\n",file="R/AKAP79.R")
message("R code written to [R/AKAP79.R]")
