library(gsl)
library(GenSA)
library(parallel)
library(dplyr)
library(tidyr)
library(tibble)
library(readr)
library(magrittr)

par.names <- c("a","b","g1","g2","n","dt")
source("functions.R")

args = commandArgs(trailingOnly=TRUE)
chunk.file=args[1]
current.gene=args[2]
max.time=as.numeric(args[3])
no.cores=as.numeric(args[4])
data.file=args[5]
out.file=args[6]

print("reading chunks")

chunks=read_csv(chunk.file)

chunks %>%
    filter(gene==current.gene) ->
    chunks

fit.data <- read_csv(data.file)

fit.data %>%
    filter(batch==1) %>%
    filter(gene == current.gene) ->
    fit.data

fit.data %>%
    left_join(chunks,by="gene") %>%
    group_by(gene) %>%
    mutate(a=2^min(mny),b=2^max(mny),g1=1e-1,g2=1e-2,n=1,dt=0) %>%
    mutate(expression=mny,sd.expression=sdy) ->
    fit.input.data

print("fitting")

out.pars <- do.untreated.fitting(data.matrix=fit.input.data,max.time=max.time,par.names=par.names)

print("writing results")
write_csv(x=out.pars,path=out.file,col_names=F)
