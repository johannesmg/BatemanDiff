library(dplyr)
library(tidyr)
library(tibble)
library(readr)

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
source("functions.R")


args = commandArgs(trailingOnly=TRUE)
no.chunks=as.numeric(args[1])
data.file=args[2]
out.file=args[3]
par.bound.file=args[4]


par.bounds <- load.bounds.from.file(par.bound.file)

all.genes <- unique(as.character(read.csv(data.file)[,"gene"]))
no.genes=length(all.genes)
chunk.delims <- unlist(lapply(1:no.chunks,rep,ceiling(no.genes/no.chunks)))[1:no.genes]

lo.pars <- matrix(rep(par.bounds[["lo.pars"]],no.genes),ncol=length(par.names),byrow=TRUE)
colnames(lo.pars) <- paste(par.names,"lo",sep="_")

hi.pars <- matrix(rep(par.bounds[["hi.pars"]],no.genes),ncol=length(par.names),byrow=TRUE)
colnames(hi.pars) <- paste(par.names,"hi",sep="_")

gene.frame <- data_frame(gene=all.genes,seed=ceiling(1e4*runif(no.genes)),chunk=chunk.delims)

chunk.matrix <- bind_cols(gene.frame,as_data_frame(lo.pars),as_data_frame(hi.pars))

write.csv(chunk.matrix,file=out.file,col.names=TRUE,row.names=FALSE)
