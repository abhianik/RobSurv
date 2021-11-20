# Compute the MDPDE, perform tests for the regression coefficiencts, and compute the DIC under the PPH model (Nandy et al., 2020)

mdpdePPH <- function(times,delta,covar,alpha,beta_ini,gamma_ini,baseline,Tol=0.00001,beta0=matrix(rep(0,length(beta_ini)),ncol=1)){
  
  
  library(psych)
  
  if(alpha<0)
    {print("Error: alpha must be a non-negative number!")}
  else if (alpha==0)
  {source("mlePPH.r")
   Out<-mlePPH(times,delta,covar,beta_ini,gamma_ini,baseline,Tol,beta0)
  return(Out)
  }
  else {
  n=nrow(covar)
  p=length(beta_ini)
  s=length(gamma_ini)
  theta_ini=c(beta_ini,gamma_ini)
  
  if (baseline=="Exponential"){
    source("Helper_functions_exp.r")
    lower=0.000000000001}
  else if (baseline=="Weibull"){
    source("Helper_functions_weibull.r")
    lower= c(0.00000000001,alpha/(alpha+3))}
  else {print("Error: Incorrect Specification of Baseline Hazard!")
        return()}
  

  dif1=1
  while(dif1>Tol)
  {
    #Minimization using optim

    h_n=function(gammas){d_n(beta_ini,gammas,alpha,times,delta,covar)}
    # opt=optim(gamma_ini, h_n,method="BFGS")
    opt=optim(gamma_ini, h_n,method="L-BFGS-B",lower =lower)
    # # print(paste("Convergence message in optim = ",opt$convergence))
    gamma_new=opt$par
      
      
    #Estimating beta using N-R
   
    dif2=1
    while(dif2>Tol)
    {
      beta_new<-beta_ini-(solve(beta_heis_func_dpd(beta_ini,gamma_new,alpha,times,delta,covar))%*%beta_grad_func_dpd(beta_ini,gamma_new,alpha,times,delta,covar))
      dif2=sqrt(sum((round(beta_new,5)-round(beta_ini,5))^2))
      beta_ini<-as.vector(beta_new)
    }
    theta_new=c(beta_new,gamma_new)
    dif1=sqrt(sum((theta_new-theta_ini)^2))

    theta_ini=theta_new
    beta_ini=beta_new
    gamma_ini=gamma_new
  }
  
gamma_est=gamma_new
beta_est=beta_new


## Compute asymptotic Standard Error in 'se' 
J=J_n(beta_est,gamma_est,alpha,times,delta,covar)
K=K_n(beta_est,gamma_est,alpha,times,delta,covar)
se=sqrt(diag(solve(J)%*%K%*%solve(J))/n)
se_beta=se[1:p]
  if(s>1){se_gamma=se[(p+1):(p+s)]} 
  else {se_gamma=se[p+1]}


## Perform Wald-type test for beta=beta0 (individually) & record P-values in 'PV'
W<-((beta_est-beta0)/se_beta)^2   #Test statistics
PV=1-pchisq(W,1) #p-values


## Compute DIC
dic=(d_n(beta_est,gamma_est,alpha,times,delta,covar)/n)+(((alpha+1)/n)*tr(solve(J)%*%K))+(1/(alpha))
  

#Table for beta and Gamma
Betas<-cbind(beta_est,se_beta,PV)
colnames(Betas)<-c("Estimate", "Standard Error", "p-value")

Gammas<-cbind(gamma_est,se_gamma)
colnames(Gammas)<-c("Estimate", "Standard Error")

#Output 
Out<-list(beta=Betas,gamma=Gammas,DIC=dic)
return(Out)
}
}
