library(survival)
#install.packages('SurvRegCensCov',dependencies = TRUE,repos = 'http://cran.rstudio.com')
library(SurvRegCensCov)
#install.packages('eha',dependencies = TRUE,repos = 'http://cran.rstudio.com')
library(eha)
#install.packages('writexl',dependencies = TRUE,repos = 'http://cran.rstudio.com')
library(writexl)

#Parameters: alpha is the shape parameter
#         : lambda is the scale parameter
#         : theta is the censoring parameter

np=function(n,lambda,b1,b2,bt,alpha,theta){
  #Generate covariates
  x1=rnorm(n,0,1)
  x2=rbinom(n,1,0.5)
  
  #Pontential switching time for every subject in the study. Every subject has the potential to move from unexposed to exposed
  a0=1
  a1=1
  a2=1
  ts=-log(runif(n))/(lambda*exp(a0+a1*x1+a2*x2))
  
  
  #Exponetial survival times
  u=runif(n)
  Ti=ifelse(-log(u)<lambda*exp(b1*x1+b2*x2)*ts,-log(u)/(lambda*exp(b1*x1+b2*x2)),
            (-log(u)-lambda*exp(b1*x1+b2*x2)*ts+lambda*exp(b1*x1+b2*x2+bt)*ts)/(lambda*exp(b1*x1+b2*x2+bt)))
  
  #censoring times
  Ci=runif(n,0,theta)
  
  #observed times
  Zi=pmin(Ti,Ci)
  di=as.numeric(Ti<=Ci)
  
  
  #The time at which the time-dependent covariate moves from unexposed(zt=0) to exposed (zt=1)
  t0=pmin(Zi,ts)
  
  
  #exposure status:whether a subject has changes statuses or not.
  ext=-log(u)>=lambda*exp(b1*x1+b2*x2)*ts & Ci>ts
  
  #time-varying covariate
  zt=as.numeric(Ti>ts)
  
  data.frame(id=1:n,x1=x1,x2=x2,di=di,Ti=Ti,zt=zt,Ci=Ci, ts=ts, Zi=Zi, t0=t0, ext=ext)
  
}
#data.e=np(50,1,1,1,-0.5,2)
#data.e
#1-sum(data.e$di)/50



runs=10000



p=data.frame(matrix(nrow = runs, ncol =1))

#Data storage for aft weibull model
aft.wei=data.frame(matrix(nrow = runs, ncol = 5+5+5+5+5))
colnames(aft.wei)=c("lambda","alpha","b1", "b2","bt","lambda_se","alpha_se","se_b1","se_b2","se_bt","lambda_cil","alpha_cil","cil_b1","cil_b2", "cil_bt","lambda_ciu","alpha_ciu" ,"ciu_b1","ciu_b2", "ciu_bt","lambda_coverage","alpha_coverage" ,"coverage_b1","coverage_b2","coverage_bt" )





set.seed(101)
for(i in 1:runs){
  
  #10%
  data.e=np(50,1,1,1,-0.5,1,14)
  p[i,]=1-sum(data.e$di)/50
  
  #data.e=np(100,1,1,1,-0.5,1,14)
  #p[i,]=1-sum(data.e$di)/100
  
  #data.e=np(250,1,1,1,-0.5,1,14)
  #p[i,]=1-sum(data.e$di)/250
  
  #data.e=np(1000,1,1,1,-0.5,1,14)
  #p[i,]=1-sum(data.e$di)/1000
  
  
  
  
  
  ##########################################################################################  
  #30%
  #data.e=np(50,1,1,1,-0.5,1,4)
  #p[i,]=1-sum(data.e$di)/50
  
  #data.e=np(100,1,1,1,-0.5,1,4)
  #p[i,]=1-sum(data.e$di)/100
  
  #data.e=np(250,1,1,1,-0.5,1,4)
  #p[i,]=1-sum(data.e$di)/250
  
  #data.e=np(1000,1,1,1,-0.5,1,4)
  #p[i,]=1-sum(data.e$di)/1000
  
  
  
  #############################################################################################
  #45%
  #data.e=np(50,1,1,1,-0.5,1,1.75)
  #p[i,]=1-sum(data.e$di)/50
  
  #data.e=np(100,1,1,1,-0.5,1,1.75)
  #p[i,]=1-sum(data.e$di)/100
  
  #data.e=np(250,1,1,1,-0.5,1,1.75)
  #p[i,]=1-sum(data.e$di)/250
  
  #data.e=np(1000,1,1,1,-0.5,1,1.75)
  #p[i,]=1-sum(data.e$di)/1000
  
  
  

  
  #Fit AFT weibull model
  fit.aft.wei=survreg(Surv(Zi, di) ~ x1+x2+zt, data=data.e, dist = "weibull")
  wei=ConvertWeibull(fit.aft.wei)
  aft.wei[i,1:5]=wei$vars[,1] #coefficients
  aft.wei[i,6:10]=wei$vars[,2] #standard errors
  aft.wei[i,11:15]=aft.wei[i,1:5]+qnorm(.025) *aft.wei[i,6:10]
  aft.wei[i,16:20]=aft.wei[i,1:5]+qnorm(.975) *aft.wei[i,6:10] 
  aft.wei[i,21]=as.integer((aft.wei[i,11]<=1)&(aft.wei[i,16]>=1))
  aft.wei[i,22]=as.integer((aft.wei[i,12]<=1)&(aft.wei[i,17]>=1))
  aft.wei[i,23:24]=as.integer((aft.wei[i,13:14]<=1)&(aft.wei[i,18:19]>=1))
  aft.wei[i,25]=as.integer((aft.wei[i,15]<=-0.5)&(aft.wei[i,20]>=-0.5))
  
  
}

#Censoring percentage
colMeans(p)


####AFT weibull model results########
#Estimates
Estimates=colMeans(aft.wei[,1:5])
Estimates
#Bias
true_values=c(1,0.1,1,1,-0.5)
Bias=Estimates-true_values
Bias
#Estimated standard error
Est_Se=colMeans(aft.wei[,6:10])
Est_Se
#Empirical standard error
Emp_Se=sqrt(diag(var(aft.wei[,1:5])))
Emp_Se
#True_matrix
True_matrix=matrix(true_values,ncol = 5,nrow = runs, byrow = T)
#MSE
mse=colMeans((True_matrix- aft.wei[,1:5])^2)
mse
#Coverage probability
cp=colMeans(aft.wei[,21:25]) 
cp 



