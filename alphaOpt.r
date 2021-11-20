# Function to find optimal alpha, by the procedure as described in Nady et al. (2020). 

alphaOpt <- function(times,delta,covar,alpha_pilot,beta_ini,gamma_ini,baseline,Tol=0.00001,Range=seq(0,1,0.025)){
  
  if (baseline=="Exponential"){source("Helper_functions_exp.r")}
  else if (baseline=="Weibull"){source("Helper_functions_weibull.r")}
  else {print("Error: Incorrect Specification of Baseline Hazard!")
    return()}
  
  N=length(Range)
  p=length(beta_ini)+length(gamma_ini)
  amse=rep(0,N)
  out_pilot=mdpdePPH(times,delta,covar,alpha_pilot,beta_ini,gamma_ini,baseline=baseline,Tol=Tol)
  beta_pilot<-out_pilot$beta[,1]
  gamma_pilot<-out_pilot$gamma[,1]
  theta_pilot<-c(beta_pilot,gamma_pilot)

  out_mle=mdpdePPH(times,delta,covar,0.000000001,beta_ini,gamma_ini,baseline=baseline,Tol=Tol)
  beta_mle<-out_mle$beta[,1]
  gamma_mle<-out_mle$gamma[,1]
  theta_mle<-c(beta_mle,gamma_mle)
  se_mle<-c(out_mle$beta[,2],out_mle$gamma[,1])
  
  amse[1]=sum((theta_pilot-theta_mle)^2)+sum(se_mle^2)

  for(i in 2:length(Range))
  {
    out<-mdpdePPH(times,delta,covar,Range[i],beta_ini,gamma_ini,baseline=baseline,Tol=Tol)
    beta_dpd=out$beta[,1]
    gamma_dpd=out$gamma[,1]
    theta_dpd<-c(beta_dpd,gamma_dpd)
    se_dpd<-c(out$beta[,2],out$gamma[,1])
    amse[i]=sum((theta_pilot-theta_dpd)^2)+sum(se_dpd^2)
    print(i)
}

  amse
  opt_alpha=Range[which(amse==min(amse))]
  return(opt_alpha)
}
