library(coxrobust)
library(psych)
source("mdpdePPH.r")
source("alphaOpt.r")


### Example 1 - smallcell  data : https://www.rdocumentation.org/packages/emplik/versions/1.0-3/topics/smallcell -----

# Get data ---
library("emplik")
data(smallcell)
data_new=smallcell
covar=as.matrix(data_new[,1:2])
times=data_new[,3]
delta=data_new[,4]


# Semiparametric methods

result=coxr(Surv(survival, indicator) ~  arm+entry, data = data_new)
result


# Exponential Baseline #########################################


beta_ini=result$coefficients
gamma_ini=sum(delta)/sum(times)
M<-mad(times)


## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini,gamma_ini*M,baseline="Exponential",Tol=0.1)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+1,ncol = ln)
SE<-matrix(0,nrow=p+1,ncol = ln)
Pval<-matrix(0,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)


for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini*M,baseline="Exponential",Tol=0.00001)
  Est[,i]<-c(out$beta[,1],out$gamma[1])
  SE[,i]<-c(out$beta[,2],out$gamma[2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
}
  
Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
covarM=as.matrix(covar[,-1])
beta_iniM=beta_ini[-1]
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.001,Range=seq(0,1,0.01))
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.00001)
alpha_opt
out



# Weibull baseline ##########################################

beta_ini=result$coefficients
result1=survreg(Surv(survival, indicator) ~1, data = data_new,dist="weibull")
s=result1$scale
r=result1$coefficients
gamma_ini=c(exp(-r),1/s)
M<-mad(times)
gamma_ini=c(M*sum(delta)/sum(times),1/s)


## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini,gamma_ini,baseline="Weibull",Tol=0.1)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
alpha=c(0.0001,0.05,0.1,0.2,0.3,0.4,0.5)
# alpha<-0.5
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+2,ncol = ln)
SE<-matrix(0,nrow=p+2,ncol = ln)
Pval<-matrix(0,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)


for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini,baseline="Weibull",Tol=0.1)
  Est[,i]<-c(out$beta[,1],out$gamma[,1])
  SE[,i]<-c(out$beta[,2],out$gamma[,2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
  print(i)
}

Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
covarM=as.matrix(covar[,-1])
beta_iniM=beta_ini[-1]
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.1,Range=seq(0,.4,0.01))
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.1)
alpha_opt
out




# Example 2: - veteran lung cancer data (outlier = case 44): -----------------------
# http://astrostatistics.psu.edu/datasets/R/html/survival/html/veteran.html

library("survival")
library("fastDummies")
data=veteran
# write.csv(data,"Veteran.csv")
data_new=dummy_cols(data,select_columns = "celltype",remove_first_dummy = FALSE)
covar=as.matrix(data_new[,c(5,10,11,12)])
times=data_new[,3]
delta=data_new[,4]


# Semiparametric methods

result=coxr(Surv(time, status) ~  karno + celltype_smallcell+celltype_adeno+celltype_large, data = data_new)
result


# Exponential Baseline #########################################


beta_ini=result$coefficients
gamma_ini=sum(delta)/sum(times)
M<-mad(times)
M<-1

## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini/2,gamma_ini*M,baseline="Exponential",Tol=.91)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+1,ncol = ln)
SE<-matrix(0,nrow=p+1,ncol = ln)
Pval<-matrix(0,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)

for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini*M,baseline="Exponential",Tol=0.04)
  Est[,i]<-c(out$beta[,1],out$gamma[1])
  SE[,i]<-c(out$beta[,2],out$gamma[2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
}

Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
covarM=as.matrix(covar)
beta_iniM=beta_ini
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.04,Range=seq(0,1,0.01))
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.04)
alpha_opt
out


# Weibull baseline ##########################################

beta_ini=result$coefficients
result1=survreg(Surv(time, status) ~1, data = data_new,dist="weibull")
s=result1$scale
r=result1$coefficients
gamma_ini=c(exp(-r),1/s)
M<-mad(times)
M<-1
gamma_ini=c(M*sum(delta)/sum(times),1.2)


## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini,gamma_ini,baseline="Weibull",Tol=0.1)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.49)
# alpha<-0.5
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+2,ncol = ln)
SE<-matrix(0,nrow=p+2,ncol = ln)
Pval<-matrix(0,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)


for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini,baseline="Weibull",Tol=0.001)
  Est[,i]<-c(out$beta[,1],out$gamma[,1])
  SE[,i]<-c(out$beta[,2],out$gamma[,2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
  print(i)
}

Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
covarM=as.matrix(covar)
beta_iniM=beta_ini
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.1,Range=seq(0,.5,0.025))
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.1)
alpha_opt
out



# Example 3: - Criminal recidivism data: -----------------------

data=read.csv(file.choose(),header=T)
library(fastDummies)
data_new=dummy_cols(data,select_columns = c("fin","wexp"),remove_first_dummy = FALSE)
covar=data_new[,c(65,67,10)]
covar=as.matrix(covar)
delta=data_new[,3]
times=data_new[,2]


# Semiparametric methods

result=coxr(Surv(week,arrest)~fin_yes+wexp_yes+prio,data=data_new)
result


# Exponential Baseline #########################################


beta_ini=result$coefficients
gamma_ini=sum(delta)/sum(times)
M<-sd(times)
M<-10


## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini,gamma_ini*M,baseline="Exponential",Tol=.00001)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7)
# alpha<-0.9
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+1,ncol = ln)
SE<-matrix(0,nrow=p+1,ncol = ln)
Pval<-matrix(0,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)

for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini*M,baseline="Exponential",Tol=0.01)
  Est[,i]<-c(out$beta[,1],out$gamma[1])
  SE[,i]<-c(out$beta[,2],out$gamma[2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
}

Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
Range=seq(0,.7,0.01)
covarM=as.matrix(covar[,-1])
beta_iniM=beta_ini[-1]
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.001,Range)
alpha_opt
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini*M,baseline="Exponential",Tol=0.000001)
out


# Weibull baseline ##########################################

beta_ini=result$coefficients
result1=survreg(Surv(week, arrest) ~1, data = data_new,dist="weibull")
s=result1$scale
r=result1$coefficients
gamma_ini=c(exp(-r),1/s)
M<-sd(times)
M<-10
gamma_ini=c(M*sum(delta)/sum(times),1/s)


## MLE 

out<-mdpdePPH(times/M,delta,covar,0,beta_ini,gamma_ini,baseline="Weibull",Tol=0.1)
out

## MDPDE 
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)
alpha=c(0.000000001,0.05,0.1,0.2,0.3,0.4,0.5,0.7)
# alpha<-0.5
ln<-length(alpha)
p<-ncol(covar)
Est<-matrix(0,nrow=p+2,ncol = ln)
SE<-matrix(0,nrow=p+2,ncol = ln)
Pval<-matrix(NA,nrow=p,ncol = ln)
DIC<-matrix(0,nrow=1,ncol = ln)


for (i in 1:ln){
  out<-mdpdePPH(times/M,delta,covar,alpha[i],beta_ini,gamma_ini,baseline="Weibull",Tol=0.01)
  Est[,i]<-c(out$beta[,1],out$gamma[,1])
  SE[,i]<-c(out$beta[,2],out$gamma[,2])
  Pval[,i]<-out$beta[,3]
  DIC[i]<-out$DIC
  print(i)
}

Est
SE
Pval
DIC

## optimal alpha
alpha_pilot=0.5
Range=seq(0,.7,0.01)
covarM=as.matrix(covar[,1])
beta_iniM=beta_ini[1]
alpha_opt<-alphaOpt(times/M,delta,covarM,alpha_pilot,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.1,Range)
out<-mdpdePPH(times/M,delta,covarM,alpha_opt,beta_iniM,gamma_ini,baseline="Weibull",Tol=0.01)
alpha_opt
out






