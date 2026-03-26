lss.var<-function(beta,x,y,z,delta,sigma,varlb,eps,bandmu){
  n<-nrow(y)
  e<-(y[,1]-x%*%beta)/sqrt(sigma)
  mu<-x%*%beta
  eresvar<-lss.eres(e, z, y[,2],eps)
  variance<- (y[,2]*(y[,1]-mu)^2+(1-y[,2])*(eresvar*sqrt(sigma))^2)
  o<-order(mu)
  mui<-mu[o]
  variancei<-variance[o]
  dummy<- 1:n
  npoly<-1
  m<-n
  muout<-mui
  svarmiss<-lwlsld(n,bandmu,npoly,mui,variancei,m,muout)
  svar<-missvar(n,p,varlb,svarmiss)
  sigma<-svar
  sigma[dummy[o]]<-svar
  sigma
}



lwlsld<-function(n,bw,npoly,xin,yin,m,xou){
  nmax<-n
  aux1<-rep(0,nmax)
  aux2<-rep(0,nmax)
  aux3<-rep(0,nmax)
  xmat<-matrix(0,nrow=nmax,ncol=npoly+1)
  w<-rep(0,nmax)
  uu<-rep(0,m)
  you<-rep(0,m)
  aa<-1
  for(j in 1:m){
    uu<-xou[j]
    ul<-uu-aa*bw
    uh<-uu+aa*bw
    count<-0
    for(i in 1:n){
      tf<-0
      if(ul<=xin[i] && xin[i]<=uh){tf<-1}
      if(xin[i]>uh){break}
      if(tf==1){
        count<-count+1
        aux1[count]<-(xin[i]-uu)/bw
        aux2[count]<-yin[i]
        aux3[count]<-1}
    }
    if(count<3){
      if(npoly==1){
        if(count==1){you[j]<-aux2[1]}
        if(count==2){
          if(aux1[1]>=0){you[j]<-aux2[1]}
          if(aux1[2]<=0){you[j]<-aux2[2]}
          if(aux1[1]<0 && aux1[2]>0){
            xh<-aux1[2]-aux1[1]
            if(aux2[2]>aux2[1]){you[j]<--aux1[1]/xh*(aux2[2]-aux2[1])+aux2[1]}
            else{
              if(aux2[2]<aux2[1]){you[j]<-aux2[2]/xh*(aux2[1]-aux2[2])+aux2[2]}
              else{you[j]<-(aux2[1]+aux2[2])/2}
            }}
        } }}
    else{
      for(i in 2:count){
        if(aux1[i]-aux1[i-1]<1e-4){stop=1}
        else{stop=0
        break}}
      if(stop==1){you[j]=sum(aux2)/count}
      else{
        for(i in 1:count){
          w[i]<-1-(aux1[i]^2)
          if(w[i]<0){w[i]<-0}
          w[i]<-3/4*w[i]/bw}
        aux1<-aux1*bw
        for(i in 1:count){
          xmat[i,1]<-1
          for(k in 2:(npoly+1)){
            xmat[i,k]<-aux1[i]^(k-1)}}
        wxmat<-sqrt(diag(w[1:count]))%*%xmat[1:count,]
        waux2<-sqrt(diag(w[1:count]))%*%aux2[1:count]
        a<-t(wxmat)%*%wxmat
        coeff<-solve(a)%*%t(wxmat)%*%waux2
        you[j]<-coeff[npoly]}}
  }
  return(you)
}

missvar<-function(n,p,varlb,variancei){
  for(i in 1:n){
    j<-1
    while(variancei[i]==0){
      if(i<n/2){variancei[i]=variancei[j+i]
      if(j+i==n){break}}
      else{variancei[i]=variancei[i-j]
      if(i-j==0){break}}
      j<-j+1
    }
  }
  for(i in 1:n){
    if(variancei[i]<varlb){variancei[i]=varlb}}
  return(variancei)
}


lss<-function(zdummy,x, y, maxiter,tolerance,eps){
  nobs <- nrow(y)
  if(is.vector(x))
    nvar <- 1
  else
    nvar <- ncol(x)
  beta <- lss.betag(x, y[,1], y[,2], zdummy)
  niter=0
  if(is.vector(x))
    xbar <- mean(x)
  else
    xbar <- apply(x, 2, mean)
  xm <- x - rep(xbar, rep(nobs, nvar))
  xinv <- solve(t(xm) %*% xm)
  while(niter < maxiter)
  {
    niter <- niter + 1
    betaprev <- beta
    e <- y[,1] - x %*% beta
    eres <- lss.eres(e, zdummy, y[,2], eps)
    yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% beta)
    ybar <- mean(yhat)
    beta <- xinv %*% (t(xm) %*% (yhat - ybar))
    beta0 <- sum(yhat- x %*% beta)/nobs
    bb <- abs(beta)
    bb[bb<0.01] <- 0.01
    mm <- max(abs(beta - betaprev) / bb)
    if(mm < tolerance)
    {
      break
    }
  }
  beta <- matrix(c(beta0, beta),ncol=1)
  object <- list(beta=beta)
  object
}


lssnew<-function(zdummy,x, y, maxiter,tolerance,eps){
  nobs <- nrow(y)
  beta <- lss.betag(x, y[,1], y[,2], zdummy)
  niter=0
  if(is.vector(x))
    xbar <- mean(x)
  else
    xbar <- apply(x, 2, mean)
  xinv <- solve(t(x) %*% x)
  while(niter < maxiter)
  {
    niter <- niter + 1
    betaprev <- beta
    e <- y[,1] - x %*% beta
    eres <- lss.eres(e, zdummy, y[,2], eps)
    yhat <- y[,2] * y[,1] + (1 - y[,2]) * (eres + x %*% beta)
    ybar <- mean(yhat)
    beta <- xinv %*% (t(x) %*% (yhat ))
    beta0 <- sum(yhat- x %*% beta)/nobs
    bb <- abs(beta)
    bb[bb<0.01] <- 0.01
    mm <- max(abs(beta - betaprev) / bb)
    if(mm < tolerance)
    {
      break
    }
  }
  beta <- matrix(c(beta0, beta),ncol=1)
  object <- list(beta=beta)
  object
}


weight<-function(y,delta,x,maxiter=100,tolerance,varlb,eps,bandmu){
  y <- cbind(y,delta)
  x <- as.matrix(x)
  if(all(x[, 1] == 1))
    x <- x[, -1]
  if(ncol(as.matrix(y)) != 2)
    stop("Response must be a right-deltaed survival object!")
  nobs <- nrow(y)
  if(is.vector(x))
    nvar <- 1
  else
    nvar <- ncol(x)
  zdummy <- matrix(rep(1,nobs), ncol=1)
  initialbeta <- lss(zdummy,x, y, maxiter,tolerance,eps)
  beta <- initialbeta$beta
  sigma <-1
  sigma <- lss.var(beta,cbind(1,x),y,zdummy,delta,sigma,varlb,eps,bandmu)
  niter=0
  ynew<-y
  newbeta0<-0
  betaflag<-matrix(0,nrow=maxiter,ncol=(nvar+1))
  while(niter < maxiter)
  {
    niter <- niter + 1
    print("niter")
    print(niter)
    betaprev <- beta
    if(nlevels(as.factor(sigma))>2){
      xnew <- cbind(1/sqrt(sigma), x/sqrt(sigma))
      ynew[,1] <- y[,1]/sqrt(sigma)
      nvarnew <- ncol(xnew)}
    else {nvarnew <- nvar
    xnew <- x
    ynew <- y}
    beta<-lssnew(zdummy,xnew,ynew,maxiter,tolerance,eps)$beta
    if(nvarnew==nvar){
      newbeta0<-0
      betaflag[niter,2:(nvar+1)]<-beta}
    else{ newbeta0<-beta[1]
    beta<-beta[-1]
    betaflag[niter,]<-beta}
    sigma <- lss.var(beta, cbind(1,x), y, zdummy, delta, sigma,varlb,eps,bandmu)
    print("beta")
    print(beta)
    if(niter>1){
      bb <-abs(betaflag[niter,2:(nvar+1)])
      bb[bb<0.01] <- 0.01
      mm <-max(abs(betaflag[niter,2:(nvar+1)] - betaflag[(niter-1),2:(nvar+1)])/bb)
      if(mm<tolerance){break}
      else{
        stop<-0
        if(niter>9){
          for(jj in (niter-2):5){
            mm <- max(abs(betaflag[niter,2:(nvar+1)] - betaflag[jj,2:(nvar+1)]) / bb)
            if(mm<tolerance){stop<-1
            break}}
          if(stop==1){break}}}}
  }
  object<-list(beta=beta)
  object
}




weightlss<-function(y,delta, x, maxiter=100, eps=.Machine$double.eps^(2/3),
                    tolerance=0.1, varlb=0.01, bandmu)
{
  library(survival)
  library(quantreg)
  library(lss)
  if(all(x[, 1] == 1)){p<-dim(x)[2]}
  else{
    p<-dim(x)[2]+1}
  ndata<-nrow(x)
  nboot<-100
  A.boot<-matrix(0, nrow=nboot, ncol=p)
  data<-cbind(y,delta,x)
  beta<-weight(y,delta,x,maxiter=100,tolerance,varlb,eps,bandmu)
  for(ii in 1:nboot){
    print("boot")
    print(ii)
    print(ndata)
    seed<-10*ii+6
    position<-floor(runif(ndata,0,1)*ndata)+1
    datanew<-data[position,]
    out.boot<-weight(datanew[,1],datanew[,2],datanew[,3:(p+2)] ,
                     maxiter=100,tolerance,varlb,eps,bandmu)
    A.boot[ii, ]<-out.boot$beta
  }
  B.boot<-apply(A.boot,2,var)
  object<-list(beta=beta, B.boot=B.boot)
  object
}


np=function(n,lambda,b1,b2,bt,theta){
  #Generate covariates
  x1=rnorm(n,0,1)
  x2=rbinom(n,1,0.5)
  
  #Time of exposure
  a0=1
  a1=1
  a2=1
  t0=-log(runif(n))/(lambda*exp(a0+a1*x1+a2*x2))
  
  
  #Exponetial survival times
  u=runif(n)
  t.e=ifelse(-log(u)<lambda*exp(b1*x1+b2*x2)*t0,-log(u)/(lambda*exp(b1*x1+b2*x2)),
             (-log(u)-lambda*exp(b1*x1+b2*x2)*t0+lambda*exp(b1*x1+b2*x2+bt)*t0)/(lambda*exp(b1*x1+b2*x2+bt)))
  
  #censoring times
  c=runif(n,0,theta)
  
  #observed times
  time.e=pmin(t.e,c)
  status.e=as.numeric(t.e<=c)
  
  #time-varying covariate
  zt=-log(u)>=lambda*exp(b1*x1+b2*x2)*t0 & c>t0
  
  data.frame(id=1:n,x1=x1,x2=x2,zt=zt,time.e=time.e,status.e=status.e)
  
}


#install.packages('writexl',dependencies = TRUE,repos = 'http://cran.rstudio.com')
library(writexl)


library(parallel)
nworkers = detectCores()
cl = makeCluster(nworkers) 
clusterSetRNGStream(cl,iseed = 101)
clusterEvalQ(cl,library(quantreg))
clusterEvalQ(cl,library(survival))
clusterEvalQ(cl,library(SurvRegCensCov))
clusterEvalQ(cl,library(Metrics))
clusterEvalQ(cl,library(lss))




#set.seed(1)
clusterExport(cl,c("lss.var","lwlsld","missvar","lss","lssnew","weight","weightlss","np"))

results=parLapply(cl,1:1500,function(i){
  
  set.seed(i)
  
  data.e=np(50,1,1,1,-0.5,12)
  p=1-sum(data.e$status.e)/50
  
  
  y=log(data.e$time.e)
  delta=data.e$status.e
  x=cbind(1,data.e$x1,data.e$x2,data.e$zt)
  uit=weightlss(y,delta,x,bandmu =50^(-1/5))
  betas=t(-uit$beta$beta)
  se=sqrt( t(uit$B.boot))
  ci_l=betas+qnorm(.025)*se
  ci_u=betas+qnorm(.975)*se
  cp_b0=as.integer((ci_l[1]<=0)&(ci_u[1]>=0))
  cp_b1=as.integer((ci_l[2]<=1)&(ci_u[2]>=1))
  cp_b2=as.integer((ci_l[3]<=1)&(ci_u[3]>=1))
  cp_b3=as.integer((ci_l[4]<=-0.5)&(ci_u[4]>=-0.5))
  out=cbind(betas,se,ci_l,ci_u,cp_b0,cp_b1,cp_b2,cp_b3,p)
  
  #write.table(out,file=paste("C:/Users/u28122161/Desktop/Iketle programs/Output/uit.",i,"txt",sep=""),append=FALSE,row.names=FALSE,col.names=FALSE)
  
  return(out)
  
  
})

df <- data.frame(matrix(unlist(results), nrow=length(results), byrow=T),stringsAsFactors=FALSE)
colnames(df) <- c("b0", "b1","b2","bt","se_b0","se_b1","se_b2","se_bt","cil_b0","cil_b1","cil_b2","cil_bt","ciu_b0","ciu_b1","ciu_b2","ciu_bt","cv_b0","cv_b1","cv_b2","cv_bt","p")
write_xlsx(df,"C:\\Users\\User\\Desktop\\Iketle\\Exponential nonproportional\\10 percent\\waft_np50.xlsx")


++#Censoring percentage
  mean(df$p)
#Estimates
Estimates=colMeans(df[,2:4])
Estimates
#Bias
true_values=c(1,1,-0.5)
Bias=Estimates-true_values
Bias
#Estimated standard error
Est_Se=colMeans(df[,6:8])
Est_Se
#Empirical standard error
Emp_Se=sqrt(diag(var(df[,2:4])))
Emp_Se
True_matrix=matrix(true_values,ncol = 3,nrow = 1500, byrow = T)
#True_matrix
mse=colMeans((True_matrix- df[,2:4])^2)
mse
#Coverage probability
cp=colMeans(df[,18:20]) 
cp
