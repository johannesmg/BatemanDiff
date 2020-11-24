predict.genes <- function(gene.vec,t.vec,conc.vec,parameter.frame){
    model.pred.points <- as_data_frame(expand.grid(gene=gene.vec,t=t.vec,conc=conc.vec,stringsAsFactors=FALSE))
    model.pred.points %>%
    left_join(parameter.frame,by="gene") %>%
        group_by(gene) %>%
        mutate(model.response=treatment.model(t=t,conc=conc,parms=c(g1p[1],g2p[1],k01[1],n1[1],k02[1],n2[1],a[1],b[1],g1[1],g2[1],n[1],dt[1]))) %>%
        mutate(y=model.response) %>%
        dplyr::select(gene,t,conc,y) ->
        model.prediction
}

collect.kinetic.dose.data <- function(data.fr){
    data.fr %>%
        filter(conc %in% c(0,0.6)) %>%
        group_by(gene,t) %>%
        mutate(y=mny,fc=mny-mny[conc==0]) %>%
        dplyr::select(gene,t,conc,y,fc,sdy) ->
        data.kinetic

    data.fr %>%
        filter(!(conc ==0.6)) %>%
        group_by(gene,t) %>%
        mutate(fc=mny,y=mny+mny[conc==0]) %>%
        dplyr::select(gene,t,conc,y,fc,sdy) %>%
        filter(conc!=0) ->
        data.dose

    data.all <- bind_rows(data.kinetic,data.dose)
}    



combine.data.and.prediction <- function(gene.vec,t.vec,conc.vec,parameter.frame,data.fr){
    data.all <- collect.kinetic.dose.data(data.fr=data.fr)
    model.prediction <- predict.genes(gene.vec,t.vec=t.vec,conc.vec=conc.vec,parameter.frame=parameter.frame)
    model.prediction %>%
        group_by(gene,t) %>%
        mutate(fc=y-y[conc==0]) ->
        model.prediction
        
    model.prediction$sdy <- NA
    model.and.data <- bind_rows(model=model.prediction,data=data.all,.id="type")
}    



library(tidyverse)
library(gsl)
source("functions.R")


args = commandArgs(trailingOnly=TRUE)
parameter.file=args[1]
fit.data.file=args[2]
out.file=args[3]

parameters=read_csv(parameter.file)
fit.data <- read_csv(fit.data.file)
genes.to.test <- unique(parameters$gene)
model.and.data <- combine.data.and.prediction(gene.vec=genes.to.test,t.vec=unique(fit.data$t),conc.vec=unique(fit.data$conc),parameter.frame=parameters,data.fr=fit.data)

model.and.data %>%
    ungroup %>%
    filter(type=="data") %>%
    select(gene,t,conc,sdy) ->
    sd.data

model.and.data %>%
    ungroup %>%
    filter(t %in% unique(fit.data$t)) %>%
    filter(conc %in% c(0,0.6)) %>%
    select(type,gene,t,conc,y)%>%
    spread(type,y) %>%
    mutate(df=data-model) %>%
    filter(!is.na(df)) %>%
    left_join(sd.data) %>%
    mutate(df=df/sdy) %>%
    group_by(t,conc) %>%
    summarise(scaling.factor=sd(df)) ->
    scaling.time

model.and.data %>%
    ungroup %>%
    filter(t %in% unique(fit.data$t)) %>%
    filter(t==144,conc!=0.6,conc!=0) %>%
    select(type,gene,t,conc,fc)%>%
    spread(type,fc) %>%
    mutate(df=data-model) %>%
    filter(!is.na(df)) %>%
    left_join(sd.data) %>%
    mutate(df=df/sdy) %>%
    group_by(t,conc) %>%
    summarise(scaling.factor=sd(df)) ->
    scaling.conc

scaling.all=bind_rows(scaling.time,scaling.conc)


write_csv(scaling.all,path=paste(out.file,".scaling",sep=""))


fit.data %>%
    left_join(scaling.all) %>%
    mutate(sdy=sdy*scaling.factor) %>%
    select(-scaling.factor) ->
    fit.data.out

write_csv(fit.data.out,path=out.file)
