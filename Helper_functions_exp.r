
## DPD function  --------

d_n=function(beta,gamma,alpha,times,delta,covar){

  beta=as.vector(beta)
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  pt1=exp(-gamma*times[which(delta==1)]*a1)
  pt2=exp(-gamma*times[which(delta==0)]*a2)
  
  
  d1=(gamma^alpha)*(((a1^alpha)/(1+alpha))-((a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                    -((1+alpha)*(a1^alpha)*(pt1^alpha)/alpha))
  d2=(pt2^(1+alpha))-((1+alpha)*(pt2^alpha)/alpha)
  
  d=sum(d1)+sum(d2)
  return(d)
}


## gradient function for DPD for updating beta using N-R ----

beta_grad_func_dpd=function(beta,gamma,alpha,times,delta,covar){

  
  beta=as.vector(beta)
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  t1=times[which(delta==1)]
  t2=times[which(delta==0)]
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  pt1=exp(-gamma*t1*a1)
  pt2=exp(-gamma*t2*a2)
  
  g1=(gamma^alpha)*((alpha*(a1^alpha)/(1+alpha))
                    -((alpha-((1+alpha)*gamma*t1*a1)/(1+alpha))*(a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                    -((1+alpha)*(1-(gamma*t1*a1))*(a1^alpha)*(pt1^alpha)))
  
  g2=(1+alpha)*gamma*t2*a2*((pt2^alpha)-(pt2^(1+alpha)))  

  gg1=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    gg1=gg1+(g1[i]*uncen_covar[i,])
  }
  gg2=rep(0,ncol(covar))
  for(i in 1:nrow(cen_covar))
  {
    gg2=gg2+(g2[i]*cen_covar[i,])
  }
  
  g=as.vector(gg1)+as.vector(gg2)
  return(as.vector(g))
}

#hessian function for DPD for updating beta using N-R ------------

beta_heis_func_dpd=function(beta,gamma,alpha,times,delta,covar){

  beta=as.vector(beta)
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  t1=times[which(delta==1)]
  t2=times[which(delta==0)]
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  pt1=exp(-gamma*t1*a1)
  pt2=exp(-gamma*t2*a2)
  
  
  h1=(gamma^alpha)*((alpha*alpha*(a1^alpha)/(1+alpha)) 
                    - (((alpha-(1+alpha)*gamma*t1*a1)^2)*(a1^alpha)*(pt1^(1+alpha))/(1+alpha)) 
                    + (gamma*t1*((a1*pt1)^(1+alpha))) 
                    - (alpha*(1+alpha)*((1-(gamma*t1*a1))^2)*((a1*pt1)^alpha))
                    + ((1+alpha)*gamma*t1*(a1^(1+alpha))*(pt1^alpha)))
  
  h2 = (((1+alpha)*gamma*t2*a2*((pt2^alpha)-(pt2^(1+alpha))))
      + ((1+alpha)*gamma*gamma*t2*t2*a2*a2*((pt2^alpha)-(pt2^(1+alpha)))))
  
  hh1=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    hh1=hh1 + h1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))
  }
  hh2=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(cen_covar))
  {
    hh2=hh2 + h2[i]*(cen_covar[i,]%*%t(cen_covar[i,]))
  }
  
  h=as.matrix(hh1)+as.matrix(hh2)
  return(as.matrix(h))
}


#K-Matrix calculation -----------

K_n=function(beta,gamma,alpha,times,delta,covar){

  
  beta=as.vector(beta)
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  t1=times[which(delta==1)]
  t2=times[which(delta==0)]
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  pt1=exp(-gamma*t1*a1)
  pt2=exp(-gamma*t2*a2)
  n=nrow(covar)
  
  g1=(gamma^alpha)*((alpha*(a1^alpha)/(1+alpha))
                    -((alpha-((1+alpha)*gamma*t1*a1)/(1+alpha))*(a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                    -((1+alpha)*(1-(gamma*t1*a1))*(a1^alpha)*(pt1^alpha)))
  
  g2=(1+alpha)*gamma*t2*a2*((pt2^alpha)-(pt2^(1+alpha)))  
  
  gg1=( alpha*(gamma^(alpha-1))*((alpha*(a1^alpha)/(1+alpha))
                              -((alpha-((1+alpha)*gamma*t1*a1)/(1+alpha))*(a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                              -((1+alpha)*(1-(gamma*t1*a1))*(a1^alpha)*(pt1^alpha)))
      + ((gamma^alpha)*t1*(a1^(alpha+1))*((pt1^(1+alpha))-((1+alpha)*(pt1^alpha)))))
  
  gg2=(1+alpha)*t2*a2*((pt2^alpha)-(pt2^(1+alpha)))          
  
  
  K11=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    K11=K11 + g1[i]*g1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))
  }
  for(i in 1:nrow(cen_covar))
  {
    K11=K11 + g2[i]*g2[i]*(cen_covar[i,]%*%t(cen_covar[i,]))
  }
  
  K12=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    K12=K12+ (g1[i]*gg1[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    K12=K12+(g2[i]*gg2[i]*cen_covar[i,])
  }
  
  K22=sum(gg1^2)+sum(gg2^2)

  K1=cbind(K11, as.vector(K12))
  K2=c(K12,K22)
  K=rbind(K1,K2)
  K=K/n
  return(as.matrix(K))
}

#J-matrix calculation ---------------

J_n=function(beta,gamma,alpha,times,delta,covar){

  
  beta=as.vector(beta)
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  t1=times[which(delta==1)]
  t2=times[which(delta==0)]
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  pt1=exp(-gamma*t1*a1)
  pt2=exp(-gamma*t2*a2)
  n=nrow(covar)
  
  h1=(gamma^alpha)*((alpha*alpha*(a1^alpha)/(1+alpha)) 
                    - (((alpha-(1+alpha)*gamma*t1*a1)^2)*(a1^alpha)*(pt1^(1+alpha))/(1+alpha)) 
                    + (gamma*t1*((a1*pt1)^(1+alpha))) 
                    - (alpha*(1+alpha)*((1-(gamma*t1*a1))^2)*((a1*pt1)^alpha))
                    + ((1+alpha)*gamma*t1*(a1^(1+alpha))*(pt1^alpha)))
  
  h2 = (((1+alpha)*gamma*t2*a2*((pt2^alpha)-(pt2^(1+alpha))))
      + ((1+alpha)*gamma*gamma*t2*t2*a2*a2*((pt2^alpha)-(pt2^(1+alpha)))))
  
  hg1=(alpha*(alpha-1)*(gamma^(alpha-2))*((alpha*(a1^alpha)/(1+alpha))
                               -((alpha-((1+alpha)*gamma*t1*a1)/(1+alpha))*(a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                               -((1+alpha)*(1-(gamma*t1*a1))*(a1^alpha)*(pt1^alpha)))
      + (alpha*(gamma^(alpha-1))*t1*(a1^(alpha+1))*((pt1^(1+alpha))-((1+alpha)*(pt1^alpha))))
      + ((gamma^alpha)*t1*t1*(a1^(alpha+2))*(1+alpha)*(-(pt1^(1+alpha))+(alpha*(pt1^alpha)))))
  
  hg2=(1+alpha)*t2*t2*a2*a2*(-(alpha*(pt2^alpha))+((1+alpha)*(pt2^(1+alpha))))
  
  hc1=(alpha*(gamma^(alpha-1))*((alpha*(a1^alpha)/(1+alpha))
                    -((alpha-((1+alpha)*gamma*t1*a1)/(1+alpha))*(a1^alpha)*(pt1^(1+alpha))/(1+alpha))
                    -((1+alpha)*(1-(gamma*t1*a1))*(a1^alpha)*(pt1^alpha))) 
     + (gamma^alpha)*((t1*(a1^(1+alpha))*(pt1^(1+alpha)))+(t1*a1*(alpha-((1+alpha)*gamma*t1*a1))*(a1^alpha)*(pt1^(1+alpha)))
                      +((1+alpha)*t1*(a1^(1+alpha))*(pt1^alpha))
                      +(alpha*(1+alpha)*(1-(gamma*t1*a1))*t1*(a1^(1+alpha))*(pt1^alpha))))
      
    
  hc2=(((1+alpha)*t2*a2*((pt2^alpha)-(pt2^(1+alpha)))) 
       + (alpha*(1+alpha)*gamma*t2*t2*a2*a2*(-(alpha*(pt2^alpha))+((1+alpha)*(pt2^(1+alpha)))))) 
  
  J11=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    J11=J11 + h1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))
  }
  for(i in 1:nrow(cen_covar))
  {
    J11=J11 + h2[i]*(cen_covar[i,]%*%t(cen_covar[i,]))
  }
  
  J12=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    J12=J12+ (hc1[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    J12=J12+(hc2[i]*cen_covar[i,])
  }
  
  J22=sum(hg1)+sum(hg2)
  
  J1=cbind(J11, as.vector(J12))
  J2=c(J12,J22)
  J=rbind(J1,J2)
  J=-J/n
  return(as.matrix(J))
}



#### MLE gradient function --------------

beta_grad_func_mle=function(beta,gamma,times,delta,covar){
  beta=as.vector(beta)
  return(beta_grad_func_dpd(beta,gamma,0,times,delta,covar))
}



## Likelihood hessian function

beta_heis_func_mle=function(beta,gamma,times,delta,covar){
  beta=as.vector(beta)
  return(beta_heis_func_dpd(beta,gamma,0,times,delta,covar))
}


## log-likelihood function
log_likelihood=function(beta,exp_gamma,times,delta,covar)
{
  beta=as.vector(beta)
  gamma=as.numeric(exp_gamma)
  
  uncen_covar=as.matrix(covar[which(delta==1),])
  len=length(delta[delta==1])
  ll=(len*exp_gamma)+sum(uncen_covar%*%beta)-(exp(exp_gamma)*sum(times*exp(covar%*%beta)))
  return(ll)
}
