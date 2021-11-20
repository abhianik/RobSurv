# Simulate a random sample from Parametric proportiona hazard (PPH) model with agiven baseline (Exponential or Weibull)


simPPH <- function(n,beta,gamma,baseline,Covmean = rep(0, length(beta)), Covsigma = diag(length(beta)), rateCens=0.01){
  
  library(mvtnorm)
  
    p=length(beta)
    s=length(gamma)
     
    if (baseline=="Exponential"){
      gamma1=gamma[1]
      gamma2=1}
    else if (baseline=="Weibull"){
      gamma1=gamma[1]
      gamma2=gamma[2]}
    else {print("Error: Incorrect Specification of Baseline Hazard!")
      return()}
    
    #Generate covariate and time
    covariate=rmvnorm(n, Covmean, Covsigma)
    a=as.vector(exp(covariate%*%as.vector(beta)))
    u=runif(n,0,1)
    t=(((-log(u))/a)^(1/gamma2))/gamma1
    
    #generate censoring and censored data
    cen=rexp(n,rate=rateCens)
    x=pmin(t,cen)
    delta=as.numeric(t<=cen)
    table(delta)
    data=data.frame(cbind(x),cbind(delta),covariate)
    
    #Output 
    return(data)
}
