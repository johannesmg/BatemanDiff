treatment.model <- function(t,conc,parms){
    g1p <- parms[1]
    g2p <- parms[2]
    g1 <- parms[9]
    g2 <- parms[10]
    k01 <- parms[3]
    n1 <- parms[4]
    k02 <- parms[5]
    n2 <- parms[6]
    gg1 <- create.gamma.curve(conc,g1p,g1,k01,n1)
    gg2 <- create.gamma.curve(conc,g2p,g2,k02,n2)
    bateman.model(t,gg1,gg2,parms[7:12])
}


bateman.model <- function(t,g1,g2,parms){
    A=parms[1]
    B=parms[2]
    n=parms[5]
    t=t+parms[6]
    #make use of Kummer's transformation to prevent
                                        #machine size infinities
    g1.greater <- g1>g2
    g1.smaller <- !g1.greater
    rv <- rep(0,length(g1))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- A+B*exp(-g2[g1.greater]*t[g1.greater])*t[g1.greater]^n*g1[g1.greater]^n*hyperg_1F1(a=n,b=1+n,x=-(g1[g1.greater]-g2[g1.greater])*t[g1.greater])/gamma(n+1)
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- A+B*exp(-g1[g1.smaller]*t[g1.smaller])*t[g1.smaller]^n*g1[g1.smaller]^n*hyperg_1F1(a=1,b=1+n,x=(g1[g1.smaller]-g2[g1.smaller])*t[g1.smaller])/gamma(n+1)
    }
#    log2(rv)
    log2(rv+1)
}


create.gamma.curve <- function(conc,gp,g,k,n){
                                        # this should only happen for large n
                                        # and large n's should not be taken...
    if (k^n>0 & max(conc)>0){
        rv <- (gp-g)*conc^n/(k^n+conc^n)+g
    }else{
        rv <- rep(g,length(conc))
        names(rv) <- names(conc)
    }
    rv
}



load.bounds.from.file <- function(file.name){
    p0 <- read.csv(file=file.name,stringsAsFactors=FALSE)
    rownames(p0) <- p0[,"X"]
    p0 <- p0[,-1]
    list(lo.pars=as.numeric(p0["lo.pars",]),hi.pars=as.numeric(p0["hi.pars",]),start.pars=p0["start.pars",])
}


cost.function.untreated <- function(pars,data.matrix,fn.string){
    responses <- bateman.model(data.matrix$t,rep(pars[3],length(data.matrix$t)),rep(pars[4],length(data.matrix$t)),pars)
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    sum((measured.dep-responses)^2/(measured.sd^2))
}


do.untreated.fitting <- function(data.matrix,max.time,par.names){
    parameter.vector=list()
    parameter.vector[["lo"]] <- as.numeric(data.matrix[1,paste(par.names,"lo",sep="_")])
    parameter.vector[["hi"]] <- as.numeric(data.matrix[1,paste(par.names,"hi",sep="_")])
    parameter.vector[["start"]] <- as.numeric(data.matrix[1,par.names])
    seed.no=as.numeric(data.matrix[1,"seed"])
    set.seed(seed.no)
    yy <- GenSA(par=parameter.vector[["start"]],fn=cost.function.untreated,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],data.matrix=data.matrix,control=list(max.time=max.time))
    df.0 <- data.frame(t(as.matrix(c(yy$par,yy$value))))
    colnames(df.0) <- c(par.names,"cost")
    df.0$gene=as.character(data.matrix[1,"gene"])
    df.0$seed=seed.no
    df.0[,c("gene","seed",par.names,"cost")]
}


cost.function.treated <- function(pars,data.matrix,fn.string){
    responses <- treatment.model(data.matrix[,"t"],data.matrix[,"conc"],pars)
    w.conc=which(data.matrix[,"t"]==144 & data.matrix[,"conc"]>0 & data.matrix[,"batch"]==3)
    w.utr=which(data.matrix[,"t"]==144 & data.matrix[,"batch"]==1)
    responses[w.conc] <- responses[w.conc]-responses[w.utr]
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    sum((measured.dep-responses)^2/(measured.sd^2))
}

cost.function.residual.vector <- function(pars,data.matrix,fn.string){
    responses <- treatment.model(data.matrix[,"t"],data.matrix[,"conc"],pars)
    w.conc=which(data.matrix[,"t"]==144 & data.matrix[,"conc"]>0 & data.matrix[,"batch"]==3)
    w.utr=which(data.matrix[,"t"]==144 & data.matrix[,"batch"]==1)
    responses[w.conc] <- responses[w.conc]-responses[w.utr]
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    (measured.dep-responses)/measured.sd
}


do.treated.fitting <- function(data.matrix,parameter.vector,optim.method,optim.arguments){
    if (optim.method=="lbfgsb"){
        yy=optim(par=parameter.vector[["start"]],fn=cost.function.treated,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],gr=function(x,...){-analytical.gradient.cost.function(x,...)},data.matrix=data.matrix,method="L-BFGS-B",control=optim.arguments)
    }else if (optim.method=="lm"){
        yy=nls.lm(par=parameter.vector[["start"]],fn=cost.function.residual.vector,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],jac=function(x,...){-jacobian.cost.function(x,...)},data.matrix=data.matrix,control=optim.arguments)
        yy$value=yy$deviance
    }else if (optim.method=="gensa"){
        yy <- GenSA(par=parameter.vector[["start"]],fn=cost.function.treated,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],data.matrix=data.matrix,control=optim.arguments)
        yy$message=NA
    }else if (optim.method=="nmkb"){
        yy <- tryCatch(nmkb(par=mod.par(parameter.vector[["start"]],parameter.vector),fn=pl.interface,lower=parameter.vector[["lo"]],upper=parameter.vector[["hi"]],data.matrix=data.matrix),error=function(e){list(par=parameter.vector[["start"]],value=cost.function.treated(parameter.vector[["start"]],data.matrix=data.matrix),message="NA")})
    }else stop(cat("Optim method must be one of \"lbfgsb\", \"lm\", or \"nmbk\""))
    yy
}



provide.starting.parameters <- function(gene.name,seed.no,no.reps,no.pars,par.names,bounds.frame,random.gen="randomLHS"){
    fran <- get(random.gen)
    set.seed(seed.no)
    random.pars <- fran(no.reps,no.pars)
    random.pars <- lapply(1:no.pars,function(y){10^(log10(bounds.frame$lo.observed[y])+(log10(bounds.frame$hi.observed[y])-log10(bounds.frame$lo.observed)[y])*random.pars[,y])})
    random.pars <- as_tibble(do.call("cbind",random.pars))
    colnames(random.pars) <- par.names
    bind_cols(tibble(gene=rep(gene.name,nrow(random.pars))),random.pars)
}    



