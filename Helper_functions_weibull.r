
## DPD function  --------

library(pracma)

gammaincApprx=function(x,a){
  if(a>0){
  u1=seq(from=a, to=0,by=-1)
  I=(0:(length(u1)-1))
  u=cumprod(u1)
  Ap=u/(x^I)
  g=(x^(a-1))*(exp(-x))*sum(Ap)}
  else{g<-0}
  return(g)
}



gammaincM=function(x,a){
  g=rep(0,length(x))
  # print(x)
  for (i in 1:length(x)){
    # if(is.nan(x[i])){g[i]<-Inf}
    # else{ g[i]<-gammainc(x[i],a)[1]}
    g[i]<-gammaincApprx(x[i],a)
  }
  return(as.vector(g))
}

d_n=function(beta,gamma,alpha,times,delta,covar){
  
  beta=as.vector(beta)
  gamma1=gamma[1]
  gamma2=gamma[2]
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  Lam1=((gamma1*times[which(delta==1)])^gamma2)*a1
  Lam2=((gamma1*times[which(delta==0)])^gamma2)*a2
  pt1=exp(-Lam1)
  pt2=exp(-Lam2)
  l1=gamma2*(gamma1^gamma2)*((times[which(delta==1)])^(gamma2-1))*a1
  
  # print(gamma2)
  s=alpha+1-(alpha/gamma2)
  C1=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s))/((1+alpha)^s)
  d1=C1-((1+alpha)*(l1^alpha)*(pt1^alpha)/alpha)
  d2=(pt2^(1+alpha))-((1+alpha)*(pt2^alpha)/alpha)
  
  d=sum(d1)+sum(d2)
  return(d)
}


## gradient function for DPD for updating beta using N-R ----

beta_grad_func_dpd=function(beta,gamma,alpha,times,delta,covar){
  
  
  beta=as.vector(beta)
  gamma1=gamma[1]
  gamma2=gamma[2]
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  Lam1=((gamma1*times[which(delta==1)])^gamma2)*a1
  Lam2=((gamma1*times[which(delta==0)])^gamma2)*a2
  pt1=exp(-Lam1)
  pt2=exp(-Lam2)
  l1=gamma2*(gamma1^gamma2)*((times[which(delta==1)])^(gamma2-1))*a1
  
  s=alpha+1-(alpha/gamma2)
  # print(s)
  C1=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s))/((1+alpha)^s)
  C2=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+1))/(((1+alpha)^(s+1))*(gamma1^gamma2)*a1)
  
  g1=(1+alpha)*(C1-C2-(((l1*pt1)^alpha)*(1-Lam1)))
  
  g2=(1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha)))  
  
  gg1=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    if(!is.na(g1[i])){
      gg1=gg1+(g1[i]*uncen_covar[i,])}
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
  gamma1=gamma[1]
  gamma2=gamma[2]
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  Lam1=((gamma1*times[which(delta==1)])^gamma2)*a1
  Lam2=((gamma1*times[which(delta==0)])^gamma2)*a2
  pt1=exp(-Lam1)
  pt2=exp(-Lam2)
  l1=gamma2*(gamma1^gamma2)*((times[which(delta==1)])^(gamma2-1))*a1
  
  s=alpha+1-(alpha/gamma2)
  # print(s)
  C1=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s))/((1+alpha)^s)
  C2=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+1))/(((1+alpha)^(s+1))*(gamma1^gamma2)*a1)
  C3=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+2))/(((1+alpha)^(s+2))*(gamma1^(2*gamma2))*a1*a1)
  
  
  h1=(1+alpha)*(((1+alpha)*C1)-((3+2*alpha)*C2)+((2+alpha)*C3)-(alpha*((l1*pt1)^alpha))
                +((2*alpha+1)*((l1*pt1)^alpha)*Lam1)-((alpha+1)*((l1*pt1)^alpha)*Lam1*Lam1))
  
  h2 = (((1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha))))
        + ((1+alpha)*Lam2*Lam2*(((1+alpha)*(pt2^alpha))-(alpha*(pt2^(1+alpha))))))
  
  hh1=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    if(!is.na(h1[i])){
      hh1=hh1 + h1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))}
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
  gamma1=gamma[1]
  gamma2=gamma[2]
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  Lam1=((gamma1*times[which(delta==1)])^gamma2)*a1
  Lam2=((gamma1*times[which(delta==0)])^gamma2)*a2
  pt1=exp(-Lam1)
  pt2=exp(-Lam2)
  n=nrow(covar)
  l1=gamma2*(gamma1^gamma2)*((times[which(delta==1)])^(gamma2-1))*a1
  
  s=alpha+1-(alpha/gamma2)
  C1=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s))/((1+alpha)^s)
  C2=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+1))/(((1+alpha)^(s+1))*(gamma1^gamma2)*a1)
  
  g1=(1+alpha)*(C1-C2-(((l1*pt1)^alpha)*(1-Lam1)))
  
  g2=(1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha)))  
  
  
  gg11=g1*gamma2/gamma1
  
  gg12=(1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha)))*(gamma2/gamma1)          
  
  gg21=  (1+alpha)*((C1/gamma2)+(log(gamma1*times[which(delta==1)])*(C1-C2))-(((l1*pt1)^alpha)*((1/gamma2)+(log(gamma1*times[which(delta==1)])*(1-Lam1)))))
  
  gg22=(1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha)))*(log(gamma1*times[which(delta==0)]))          
  
  
  K11=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    K11=K11 + g1[i]*g1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))
  }
  for(i in 1:nrow(cen_covar))
  {
    K11=K11 + g2[i]*g2[i]*(cen_covar[i,]%*%t(cen_covar[i,]))
  }
  
  K121=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    K121=K121+ (g1[i]*gg11[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    K121=K121+(g2[i]*gg12[i]*cen_covar[i,])
  }
  
  K122=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    K122=K122+ (g1[i]*gg21[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    K122=K122+(g2[i]*gg22[i]*cen_covar[i,])
  }

  K221=sum(gg11^2)+sum(gg12^2)
  K222=sum(gg21^2)+sum(gg22^2)
  K220=sum(gg11*gg21)+sum(gg12*gg22)
  
  K1=cbind(K11, as.vector(K121),as.vector(K122))
  K21=c(K121,K221,K220)
  K22=c(K122,K220,K222)
  K=rbind(K1,K21,K22)
  K=K/n
  return(as.matrix(K))
}

#J-matrix calculation ---------------

J_n=function(beta,gamma,alpha,times,delta,covar){
  
  
  beta=as.vector(beta)
  gamma1=gamma[1]
  gamma2=gamma[2]
  uncen_covar=as.matrix(covar[which(delta==1),])
  cen_covar=as.matrix(covar[which(delta==0),])
  a1=exp(beta%*%t(uncen_covar))
  a2=exp(beta%*%t(cen_covar))
  Lam1=((gamma1*times[which(delta==1)])^gamma2)*a1
  Lam2=((gamma1*times[which(delta==0)])^gamma2)*a2
  pt1=exp(-Lam1)
  pt2=exp(-Lam2)
  n=nrow(covar)
  l1=gamma2*(gamma1^gamma2)*((times[which(delta==1)])^(gamma2-1))*a1
  
  s=alpha+1-(alpha/gamma2)
  C1=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s))/((1+alpha)^s)
  C2=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+1))/(((1+alpha)^(s+1))*(gamma1^gamma2)*a1)
  C3=((gamma1*gamma2)^alpha)*(a1^(-alpha/gamma2))*(gammaincM((1+alpha)*Lam1,s+2))/(((1+alpha)^(s+2))*(gamma1^(2*gamma2))*a1*a1)
  C0=log(gamma1*times[which(delta==0)])
  C01=log(gamma1*times[which(delta==1)])*times[which(delta==1)]/2
  
  h1=(1+alpha)*(((1+alpha)*C1)-((3+2*alpha)*C2)+((2+alpha)*C3)-(alpha*((l1*pt1)^alpha))
                +((2*alpha+1)*((l1*pt1)^alpha)*Lam1)-((alpha+1)*((l1*pt1)^alpha)*Lam1*Lam1))
  
  h2 = (((1+alpha)*Lam2*((pt2^alpha)-(pt2^(1+alpha))))
        + ((1+alpha)*Lam2*Lam2*(((1+alpha)*(pt2^alpha))-(alpha*(pt2^(1+alpha))))))
  
  
  hg11=(1+alpha)*(gamma2/(gamma1*gamma1))*(((gamma2*(1+alpha)-1)*C1)-((gamma2*(3+2*alpha)-1)*C2)+((1+alpha)*gamma2*C3)
                      -(((l1*pt1)^alpha)*((alpha*gamma2-1)-((gamma2*(1+2*alpha)-1)*Lam1)+(alpha*gamma2*Lam1*Lam1))))
  
  
  hg12=(1+alpha)*((Lam2*Lam2*((gamma2/gamma1)^2)*(((1+alpha)*(pt2^(1+alpha)))-(alpha*(pt2^alpha))))
                  +(((pt2^alpha)-(pt2^(1+alpha)))*Lam2*(gamma2/gamma1)*((gamma2-1)/gamma1)))

  hg21=((1+alpha)*(alpha/(gamma2*gamma2)) + (2*(1+alpha)*(1+alpha)*C01*(C1-C2)/gamma2)
        + ((1+alpha)*(1+alpha)*C01*(C1-C2-C2+C3))
        -(((l1*pt1)^alpha)*(1+alpha)*(((alpha-1)/(gamma2*gamma2))+(2*alpha*C01*(1-Lam1)/gamma2)
                                      +(alpha*C01*C01*(1-Lam1)*(1-Lam1)))))
  
  
  hg22=(1+alpha)*C0*C0*((Lam2*Lam2*(((1+alpha)*(pt2^(1+alpha)))-(alpha*(pt2^alpha))))
                        +(((pt2^alpha)-(pt2^(1+alpha)))*Lam2))
  
  hgc1=(1+alpha)*(((1+alpha)*(C1-C2)/gamma1)+((1+alpha)*gamma2*C01*C1/gamma1)
                  -((3+2*alpha)*gamma2*C01*C2/gamma1)-((1+alpha)*gamma2*C01*C3/gamma1)
        + (((l1*pt1)^alpha)*(1/gamma1)*((alpha*(1-Lam1))+(gamma2*C01*(alpha-((2*alpha+1)*Lam1)+(alpha*Lam1*Lam1))))))

  hgc2=(1+alpha)*((Lam2*Lam2*((gamma2/gamma1)*C0)*(((1+alpha)*(pt2^(1+alpha)))-(alpha*(pt2^alpha))))
                  +(((pt2^alpha)-(pt2^(1+alpha)))*Lam2*(1/gamma1)*((gamma2*C0)+1)))
  
  hc11=(1+alpha)*(gamma2/gamma1)*(((1+alpha)*C1)-((3+2*alpha)*C2)+((1+alpha)*C3)
                                  -(((l1*pt1)^alpha)*(alpha-((2*alpha+1)*Lam1)+(alpha*Lam1*Lam1))))
  
  
  hc12=(1+alpha)*((Lam2*Lam2*((gamma2/gamma1))*(((1+alpha)*(pt2^(1+alpha)))-(alpha*(pt2^alpha))))
                  +(((pt2^alpha)-(pt2^(1+alpha)))*Lam2*(gamma2/gamma1)))
  
  hc21=(1+alpha)*(((1+alpha)*(C1-C2)/gamma1)+((1+alpha)*C01*C1)-((3+2*alpha)*C01*C2)+((1+alpha)*C01*C3)
                  + (((l1*pt1)^alpha)*((alpha*(1-Lam1)/gamma2)+(C01*(alpha-((2*alpha+1)*Lam1)+(alpha*Lam1*Lam1))))))
  
  hc22=(1+alpha)*((Lam2*Lam2*(C0)*(((1+alpha)*(pt2^(1+alpha)))-(alpha*(pt2^alpha))))
                  +(((pt2^alpha)-(pt2^(1+alpha)))*Lam2*C0))
  J11=matrix(0,ncol(covar),ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    J11=J11 + h1[i]*(uncen_covar[i,]%*%t(uncen_covar[i,]))
  }
  for(i in 1:nrow(cen_covar))
  {
    J11=J11 + h2[i]*(cen_covar[i,]%*%t(cen_covar[i,]))
  }
  
  J121=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    J121=J121+ (hc11[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    J121=J121+(hc12[i]*cen_covar[i,])
  }
  J122=rep(0,ncol(covar))
  for(i in 1:nrow(uncen_covar))
  {
    J122=J122+ (hc21[i]*uncen_covar[i,])
  }
  for(i in 1:nrow(cen_covar))
  {
    J122=J122+(hc22[i]*cen_covar[i,])
  }
  
  J221=sum(hg11)+sum(hg12)
  J222=sum(hg21)+sum(hg22)
  J220=sum(hgc1)+sum(hgc2)
  
  # print(J121)
  
  J1=cbind(J11, as.vector(J121),as.vector(J122))
  J21=c(J121,J221,J220)
  J22=c(J122,J220,J222)
  J=rbind(J1,J21,J22)
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
  gamma2=exp(exp_gamma[2])
  gamma1=exp(exp_gamma[1])
  
  uncen_covar=as.matrix(covar[which(delta==1),])
  len=length(delta[delta==1])
  L=(len*log(gamma2*(gamma1^gamma2)))+sum(uncen_covar%*%beta)+((gamma2-1)*sum(log(times[which(delta==1)])))-sum((gamma1*times)^(gamma2)*exp(covar%*%beta))
  return(L)
}
