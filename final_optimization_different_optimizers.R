library(parallel)
library(gsl)
library(lhs)                                       
library(dplyr)
library(tidyr)
library(readr)
library(minpack.lm)
library(dfoptim)
library(GenSA)

source("functions.R")
source("gradient_functions.R")

args = commandArgs(trailingOnly=TRUE)
data.file=args[1]
in.file=args[2]
parameter.bounds.file=args[3]
out.file=args[4]
no.cores=as.numeric(args[5])
no.reps.2=as.numeric(args[6])
optim.method=args[7]
fnscale=as.numeric(as.character(args[8]))
factr=as.numeric(as.character(args[9]))
maxit=as.numeric(as.character(args[10]))
maxtime=as.numeric(as.character(args[11]))

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
no.pars=length(par.names)

parm.bounds <- load.bounds.from.file(parameter.bounds.file)
names(parm.bounds) <- c("lo","hi","start")

bounds.frame=tibble(parameter=par.names,lo=parm.bounds[["lo"]],hi=parm.bounds[["hi"]])


in.pars.close=read_csv(in.file,col_names=TRUE)
current.gene=in.pars.close$gene[1]
seed.no=in.pars.close$seed[1]

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

in.pars.close=in.pars.close[which.min(in.pars.close$cost)[1],]

in.pars.close %>%
    select(-seed,-no.reps,-cost) %>%
    gather("parameter","value",g1p:dt) %>%
    mutate(lo.observed=value-value*1e-3,hi.observed=value+value*1e-3) ->
    observed.frame

bounds.frame %>%
    select(parameter,lo,hi) %>%
    left_join(observed.frame) %>%
    mutate(lo.observed=ifelse(is.na(lo.observed),lo,lo.observed)) %>%
    mutate(hi.observed=ifelse(is.na(hi.observed),hi,hi.observed)) %>%
    mutate(lo.observed=pmax(lo.observed,lo)) %>%
    mutate(hi.observed=pmin(hi.observed,hi)) ->
    bounds.frame.new


in.pars=provide.starting.parameters(current.gene,seed.no,no.reps=no.reps.2,no.pars=no.pars,par.names=par.names,bounds.frame.new,random.gen="randomLHS")

if (optim.method=="lbfgsb"){
    optim.arguments=list(trace=0,fnscale=fnscale,factr=factr,maxit=maxit)
}else if (optim.method=="lm"){
    optim.arguments=nls.lm.control(ftol=1e-25,ptol=1e-16,gtol=1e-15,maxiter=maxit)
}else if (optim.method=="nmkb"){
    optim.arguments=list(trace=NA,fnscale=NA,factr=NA,maxit=NA)
}else if (optim.method=="gensa"){
    optim.arguments=list(max.time=maxtime)
}



out.pars <- mclapply(1:nrow(in.pars),function(x){parm.bounds[["start"]] <- as.numeric(in.pars[x,2:13]);do.treated.fitting(data.matrix=all.data,parameter.vector=parm.bounds,optim.method=optim.method,optim.arguments=optim.arguments)},mc.cores=no.cores)

parameters=lapply(out.pars,"[[","par")
costs=lapply(out.pars,"[[","value")
messages=lapply(out.pars,"[[","message")

out.object=as_tibble(cbind(do.call("rbind",parameters),do.call("rbind",costs)))
out.object=bind_cols(out.object,tibble(message=do.call("rbind",messages)))
colnames(out.object) <- c(par.names,"cost","message")

if (optim.method=="lbfgsb"){
    out.object=bind_cols(out.object,as_tibble(optim.arguments)[rep(1,nrow(out.object)),])
}else if (optim.method=="lm"){
    out.object=bind_cols(out.object,as_tibble(optim.arguments[c("ftol","ptol","gtol","maxiter")])[rep(1,nrow(out.object)),])
}else if (optim.method=="nmkb"){
    out.object=bind_cols(out.object,as_tibble(optim.arguments)[rep(1,nrow(out.object)),])
}else if (optim.method=="gensa"){
    out.object=bind_cols(out.object,as_tibble(optim.arguments)[rep(1,nrow(out.object)),])
}


out.pars <- bind_cols(tibble(gene=rep(current.gene,nrow(out.object)),seed=rep(seed.no,nrow(out.object)),no.reps=rep(no.reps.2,nrow(out.object)),method=rep(optim.method,nrow(out.object))),as_tibble(out.object))

write_csv(out.pars,path=out.file,col_names=TRUE)


