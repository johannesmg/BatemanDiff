library(parallel)
library(lhs)
library(gsl)
library(tidyverse)
source("functions.R")
source("gradient_functions.R")
library(cowplot)
library(numDeriv)

args = commandArgs(trailingOnly=TRUE)
parameter.bounds.file=args[1]
fit.data.file=args[2]
current.gene=args[3]

par.names <- c("g1p","g2p","k01","n1","k02","n2","a","b","g1","g2","n","dt")
no.pars=length(par.names)
parm.bounds <- load.bounds.from.file(parameter.bounds.file)
names(parm.bounds) <- c("lo","hi","start")


read_csv(fit.data.file) %>%
    filter(gene==current.gene) ->
    data.in
data.in[,c("batch","t","conc","mny","sdy")] %>%
    as.matrix ->
    data.in

colnames(data.in)=c("batch","t","conc","expression","sd.expression")

bounds.frame=tibble(parameter=par.names,lo.observed=parm.bounds[["lo"]],hi.observed=parm.bounds[["hi"]])
starting.parameters=provide.starting.parameters(gene.name=current.gene,seed.no=18931,no.reps=1000,no.pars=12,par.names=par.names,bounds.frame=bounds.frame)

emp=as_tibble(do.call("rbind",(lapply(1:nrow(starting.parameters),function(x){cat(paste("...",sprintf("%.2f",x/nrow(starting.parameters)),sep=""));grad(cost.function.treated,x=as.numeric(starting.parameters[x,-1]),side=rep(1,12),data.matrix=data.in,method="Richardson")}))))
theo=as_tibble(do.call("rbind",(lapply(1:nrow(starting.parameters),function(x){cat(paste("...",sprintf("%.2f",x/nrow(starting.parameters)),sep=""));analytical.gradient.cost.function(as.numeric(starting.parameters[x,-1]),data.in)}))))

colnames(emp)=par.names
emp$type="numerical"
emp$run=1:nrow(emp)
colnames(theo)=par.names
theo$type="analytical"
theo$run=1:nrow(theo)

emp=gather(emp,"parameter","value",g1p:dt)
theo=gather(theo,"parameter","value",g1p:dt)

starting.parameters$run=1:nrow(starting.parameters)

bind_rows(emp,theo) %>%
    left_join(starting.parameters) %>%
    mutate(value=log10(abs(value))) %>%
    spread(type,value) %>%
    filter(!is.infinite(numerical)) %>%
    ggplot(aes(x=numerical,y=analytical))+
    geom_abline(slope=1,intercept=0)+
    geom_point(size=0.2)+
    facet_wrap(~parameter)+
    theme_minimal()+
    scale_colour_gradient(low="darkblue",high="gold")->
    p1

save_plot(plot=p1,filename=paste("gradient_test_",current.gene,".pdf",sep=""))

