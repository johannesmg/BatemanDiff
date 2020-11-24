

hyperg_1F1.kummer.product.1 <- function(t,g1,g2,n){
    g1.greater <- g1>g2
    g1.smaller <- !g1.greater
    rv <- rep(0,length(g1))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- exp(-(g2[g1.greater]*t[g1.greater]))*hyperg_1F1(a=n,b=1+n,x=-(g1[g1.greater]-g2[g1.greater])*t[g1.greater])
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- exp(-(g1[g1.smaller]*t[g1.smaller]))*hyperg_1F1(a=1,b=1+n,x=(g1[g1.smaller]-g2[g1.smaller])*t[g1.smaller])
    }
    rv
}


hyperg_1F1.kummer.product.2 <- function(t,g1,g2,n){
    g1.greater <- g1>g2
    g1.smaller <- !g1.greater
    rv <- rep(0,length(g1))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- exp(-g2[g1.greater]*t[g1.greater])*g1[g1.greater]^n*hyperg_1F1(a=n,b=1+n,x=-(g1[g1.greater]-g2[g1.greater])*t[g1.greater])
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- exp(-g1[g1.smaller]*t[g1.smaller])*g1[g1.smaller]^n*hyperg_1F1(a=1,b=1+n,x=(g1[g1.smaller]-g2[g1.smaller])*t[g1.smaller])
    }
    rv
}

hyperg_1F1.kummer.minimal <- function(a,b,x){
    g1.greater <- x<0
    g1.smaller <- !g1.greater
    rv <- rep(0,length(x))
    if (sum(g1.greater)>0){
        rv[g1.greater] <- hyperg_1F1(a=a,b=b,x=(x[g1.greater]))
    }
    if (sum(g1.smaller)>0){
        rv[g1.smaller] <- exp(x[g1.smaller])*hyperg_1F1(a=b-a,b=b,x=-(x[g1.smaller]))
    }
    rv
}

integrand=function(a,b,c,d,t,x){
    t^(a-1)*(1-t)^(b-a-1)*hyperg_1F1.kummer.minimal(d,c,x*t)
}


hyperg_2F2 <- function(p,q,g1,g2,t){
    a=p[1]
    d=p[2]
    b=q[1]
    c=q[2]
    n=p[1]
    z1=(-g1 + g2)*t
    z0=exp(-g2*t)*1/(n*gamma(1 + n))*t^n*g1^n
    ww=which(z0<1e-200)
    rix=1:length(z1)
    rix=setdiff(rix,ww)
    rv=rep(0,length(z1))
    if (length(rix)>0){
        rv[rix]=sapply(rix,function(z){gamma(b)/(gamma(a)*gamma(b-a))*integrate(function(y){z0[z]*integrand(a,b,c,d,y,z1[z])},0,1,stop.on.error=FALSE)$value})
    }
    rv
}



EulerGamma=-digamma(1)

HarmonicNumber=function(n){
    integrate(function(t){(1-t^n)/(1-t)},0,1)$value
}



derivn=function(t,g1,g2,n){
    t3=1 + EulerGamma*n - n*HarmonicNumber(n) + n*log(g1*t)
    t4=1/(n*gamma(1 + n))*t^n*g1^n
    t2=hyperg_1F1.kummer.product.1(t,g1,g2,n)
    t1=-hyperg_2F2(c(n,n),c(1 + n,1 + n),g1,g2,t)
    rv <- t1+t4*t2*t3
                                        #limit for t->0 is 0
    rv[t==0]=0
    rv
}


derivg1=function(t,g1,g2,n){
    equal.rates=g1==g2
                                        #if (g2==g1){
    if (length(equal.rates)>0){
        g2[equal.rates]=g2[equal.rates]+1e-8
    }
    g1^( - 1)*(g1 - g2)^(-1)*t^n*((g1^n)*exp(-g1*t)*g1 - g2*hyperg_1F1.kummer.product.2(t,g1,g2,n))/gamma(n)
    }

derivg2=function(t,g1,g2,n){
    equal.rates=g1==g2
    if (length(equal.rates)>0){
        g2[equal.rates]=g2[equal.rates]+1e-8
    }
    (g1 - g2)^(-1)*t^n*(g1^(n)*exp(-g1*t)*(-1) + (n+(-g1+g2)*t)*n^(-1)*hyperg_1F1.kummer.product.2(t,g1,g2,n))/gamma(n)
    }

derivdt=function(t,g1,g2,n){
     t^(n - 1)*(g1^n*exp(-g1*t)*1-n^(-1)*g2*t*hyperg_1F1.kummer.product.2(t,g1,g2,n))/gamma(n)
}

derivhill_gp=function(c,gp,k,n,g){
    c^n/(c^n + k^n)
}

derivhill_n=function(c,gp,k,n,g){
    rv= -(((g - gp)*(c*k)^n*log(c/k))/(c^n + k^n)^2)
    rv[c==0]=0
    rv
}

derivhill_g=function(c,gp,k,n,g){
    1 - c^n/(c^n + k^n)
}

derivhill_k=function(c,gp,k,n,g){
    -((c^n*(-g + gp)*k^(-1+n)*n)/(c^n + k^n)^2)
}


analytical.gradient.treatment.model=function(t,conc,batch,pars){
    g1p <- pars[1]
    g2p <- pars[2]
    k01 <- pars[3]
    n1 <- pars[4]
    k02 <- pars[5]
    n2 <- pars[6]
    A <- pars[7]
    B <- pars[8]
    g1 <- pars[9]
    g2 <- pars[10]
    n <- pars[11]
    dt <- pars[12]
    gg1 <- create.gamma.curve(conc,g1p,g1,k01,n1)
    gg2 <- create.gamma.curve(conc,g2p,g2,k02,n2)
    w.conc=which(t==144 & conc>0 & batch==3)
    w.utr=which(t==144 & batch==1)
    t.mod=t+dt
    model.output=treatment.model(t,conc,pars)
                                        #take care of time shift in argument!
    pre.factor=1/log(2)*1/(exp(log(2)*model.output))
    pre.factor.utr=1/log(2)*1/(exp(log(2)*model.output[w.utr]))
                                        #g1p
    d.g1p=pre.factor*B*derivg1(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_gp(conc,g1p,k01,n1,g1)
    #g2p
    d.g2p=pre.factor*B*derivg2(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_gp(conc,g2p,k02,n2,g2)
    #k01
    d.k01=pre.factor*B*derivg1(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_k(conc,g1p,k01,n1,g1)
    #n1
    d.n1=pre.factor*B*derivg1(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_n(conc,g1p,k01,n1,g1)
    #k02
    d.k02=pre.factor*B*derivg2(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_k(conc,g2p,k02,n2,g2)
    #n2
    d.n2=pre.factor*B*derivg2(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_n(conc,g2p,k02,n2,g2)
    #A
    utr.correction=pre.factor.utr
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    # For the concentration series data (batch 3) we use fold changes so we have 
    # to take this into account when calculating the derivative
    d.A=pre.factor-utr.correction
    #B
#    utr.correction=pre.factor.utr
    utr.correction=pre.factor.utr*(exp(log(2)*model.output[w.utr])-A-1)/B
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    d.B=pre.factor*((exp(log(2)*model.output)-A-1)/B)-utr.correction
    #g1
    deriv.utr=derivg1(t=144+dt,g1=g1,g2=g2,n=n)
    utr.correction=pre.factor.utr*B*deriv.utr
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    d.g1=pre.factor*B*derivg1(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_g(conc,g1p,k01,n1,g1)-utr.correction
                                        #g2
    deriv.utr=derivg2(t=144+dt,g1=g1,g2=g2,n=n)
    utr.correction=pre.factor.utr*B*deriv.utr
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    d.g2=pre.factor*B*derivg2(t=t.mod,g1=gg1,g2=gg2,n=n)*derivhill_g(conc,g2p,k02,n2,g2)-utr.correction
    #n
    deriv.utr=derivn(t=144+dt,g1=g1,g2=g2,n=n)
    utr.correction=pre.factor.utr*B*deriv.utr
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    d.n=pre.factor*B*derivn(t=t.mod,g1=gg1,g2=gg2,n=n)-utr.correction
                                        #dt
    deriv.utr=derivdt(t=144+dt,g1=g1,g2=g2,n=n)
    utr.correction=pre.factor.utr*B*deriv.utr
    utr.correction=utr.correction*as.numeric(1:length(t) %in% w.conc)
    d.dt=pre.factor*B*derivdt(t=t.mod,g1=gg1,g2=gg2,n=n)-utr.correction
    cbind(d.g1p,d.g2p,d.k01,d.n1,d.k02,d.n2,d.A,d.B,d.g1,d.g2,d.n,d.dt)
}

numerical.gradient <- function(pars.in,dx){
    ix=which(data.in[,"batch"]%in%c(1,2,3))
    comp.grad=c()
    for (i in 1:12){
        pars.alt=pars.in
        pars.alt[i]=pars.alt[i]+dx
        comp.grad[i]=(cost.function.treated(pars.in,data.in[ix,,drop=FALSE])-cost.function.treated(pars.alt,data.in[ix,,drop=FALSE]))/dx
    }
    comp.grad
}

analytical.gradient.cost.function <- function(pars,data.matrix){
    responses <- treatment.model(data.matrix[,"t"],data.matrix[,"conc"],pars)
    w.conc=which(data.matrix[,"t"]==144 & data.matrix[,"conc"]>0 & data.matrix[,"batch"]==3)
    w.utr=which(data.matrix[,"t"]==144 & data.matrix[,"batch"]==1)
    responses[w.conc] <- responses[w.conc]-responses[w.utr]
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    term.1=2*(measured.dep-responses)/(measured.sd^2)
    term.2=analytical.gradient.treatment.model(data.matrix[,"t"],data.matrix[,"conc"],data.matrix[,"batch"],pars)
    as.numeric(matrix(term.1,nrow=1) %*% term.2)
}



jacobian.cost.function <- function(pars,data.matrix){
    responses <- treatment.model(data.matrix[,"t"],data.matrix[,"conc"],pars)
    w.conc=which(data.matrix[,"t"]==144 & data.matrix[,"conc"]>0 & data.matrix[,"batch"]==3)
    w.utr=which(data.matrix[,"t"]==144 & data.matrix[,"batch"]==1)
    responses[w.conc] <- responses[w.conc]-responses[w.utr]
    measured.dep <- data.matrix[,"expression"]
    measured.sd <- data.matrix[,"sd.expression"]
    term.1=1/(measured.sd)
    term.2=analytical.gradient.treatment.model(data.matrix[,"t"],data.matrix[,"conc"],data.matrix[,"batch"],pars)
    term.1*term.2
}

