# Compute the MLE, perform tests for the regression coefficiencts, and compute the DIC under the PPH model (Nandy et al., 2020)

mlePPH <- function(times,delta,covar,beta_ini,gamma_ini,baseline,Tol=0.00001,beta0=matrix(rep(0,length(beta_ini)),ncol=1)){
  
  n=nrow(covar)
  p=length(beta_ini)
  s=length(gamma_ini)
  theta_ini=c(beta_ini,gamma_ini)
  if (baseline=="Exponential"){source("Helper_functions_exp.r")}
  else if (baseline=="Weibull"){source("Helper_functions_weibull.r")}
  else {print("Error: Incorrect Specification of Baseline Hazard!")
    return()}
  
  
  dif1=1
  while(dif1>Tol)
  {
    #Minimization using optim
    
    h_n=function(egammas){-log_likelihood(beta_ini,egammas,times,delta,covar)}
    opt=optim((gamma_ini), h_n, method="BFGS")
    # print(paste("Convergence message in optim = ",opt$convergence))
    gamma_new=exp(opt$par)

    
    
    #Estimating beta using N-R
    
    dif2=1
    while(dif2>Tol)
    {
      beta_new=beta_ini-(solve(beta_heis_func_mle(beta_ini,gamma_new,times,delta,covar))%*%beta_grad_func_mle(beta_ini,gamma_new,times,delta,covar))
      dif2=sqrt(sum((round(beta_new,5)-round(beta_ini,5))^2))
      beta_ini=beta_new
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
  J=J_n(beta_est,gamma_est,0,times,delta,covar)
  K=K_n(beta_est,gamma_est,0,times,delta,covar)
  se=sqrt(diag(solve(J)%*%K%*%solve(J))/n)
  se_beta=se[1:p]
  if(s>1){se_gamma=se[(p+1):(p+s)]} else {se_gamma=se[p+1]}
  
  
  ## Perform Wald-type test for beta=beta0 (individually) & record P-values in 'PV'
  W<-((beta_est-beta0)/se_beta)^2   #Test statistics
  PV=1-pchisq(W,1) #p-values
  
  
  ## Compute DIC
  dic=(-log_likelihood(beta_est,gamma_est,times,delta,covar)/n)+((1/n)*tr(solve(J)%*%K))
  
  
  #Table for beta and Gamma
  Betas<-cbind(beta_est,se_beta,PV)
  colnames(Betas)<-c("Estimate", "Standard Error", "p-value")
  
  Gammas<-cbind(gamma_est,se_gamma)
  colnames(Gammas)<-c("Estimate", "Standard Error")
  
  #Output 
  Out<-list(beta=Betas,gamma=Gammas,DIC=dic)
  return(Out)
}
