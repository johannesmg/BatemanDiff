library(parallel)
library(gsl)
library(lhs)                                       
library(dplyr)
library(tidyr)
library(readr)

source("functions.R")

args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
untreated.fit.file=args[2]
parameter.bounds.file=args[3]
out.file=args[4]
no.cores=as.numeric(args[5])
no.reps=as.numeric(args[6])

seed.no=ceiling(runif(1)*1e5)

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
no.pars=length(par.names)
parm.bounds <- load.bounds.from.file(parameter.bounds.file)
names(parm.bounds) <- c("lo","hi","start")


bounds.frame=tibble(parameter=par.names,lo=parm.bounds[["lo"]],hi=parm.bounds[["hi"]])
untreated.fit=read_csv(untreated.fit.file,col_names=F)
colnames(untreated.fit) <- c("gene",par.names[7:12],"cost")
current.gene=untreated.fit$gene[1]


untreated.fit %>%
    mutate(a=ifelse(a==0,1,a)) %>%
    mutate(b=ifelse(b==0,1,b)) %>%
    select(-cost,-gene) %>%
    gather("parameter","value",a:dt) %>%
    mutate(lo.observed=value/5,hi.observed=value*5) %>%
    select(-value) ->
    observed.frame
    
observed.frame %>%
    filter(parameter %in% c("g1","g2")) %>%
    mutate(parameter=paste0(parameter,"p")) %>%
    bind_rows(observed.frame) ->
    observed.frame

bounds.frame %>%
    left_join(observed.frame) %>%
    mutate(lo.observed=ifelse(is.na(lo.observed),lo,lo.observed)) %>%
    mutate(hi.observed=ifelse(is.na(hi.observed),hi,hi.observed)) %>%
    mutate(lo.observed=pmax(lo.observed,lo)) %>%
    mutate(hi.observed=pmin(hi.observed,hi)) ->
    bounds.frame

parm.bounds$start <- NA

all.data=read_csv(data.file,guess_max=3e4)
all.data=all.data %>%
  ungroup %>%
  filter(gene == current.gene)

all.data <- all.data[,c("t","mny","sdy","conc","batch","gene")]
colnames(all.data) <- c("t","expression","sd.expression","conc","batch","gene")

all.data %>%
select(-gene)->
all.data

all.data=as.matrix(all.data)

in.pars = provide.starting.parameters(current.gene,seed.no,no.reps=no.reps,no.pars=no.pars,par.names=par.names,bounds.frame,random.gen="randomLHS")

out.pars <- mclapply(1:nrow(in.pars),function(x){parm.bounds[["start"]] <- as.numeric(in.pars[x,2:13]);do.treated.fitting(data.matrix=all.data,parameter.vector=parm.bounds)},mc.cores=no.cores)

out.pars <- do.call("rbind",out.pars)
colnames(out.pars) <- c(par.names,"cost")

out.pars <- bind_cols(tibble(gene=rep(current.gene,nrow(out.pars)),seed=rep(seed.no,nrow(out.pars)),no.reps=rep(no.reps,nrow(out.pars))),as_tibble(out.pars))

out.pars=out.pars[which.min(out.pars$cost)[1],]

write_csv(out.pars,path=out.file,col_names=FALSE)

