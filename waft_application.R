lss.var<-function(beta,x,y,z,delta,adjust,sigma){
  n<-nrow(y)
  e<-y[,1]-x%*%beta-adjust
  mu<-x%*%beta+adjust
  eresvar<-lss.eres(e, z, y[,2],eps)
  #variance<- (y[,2]*(y[,1]-mu)^2+(1-y[,2])*(eresvar)^2)
  variance<- (y[,2]*(y[,1]-mu)^2+(1-y[,2])*(eresvar*sqrt(sigma))^2)
  
  o<-order(mu)
  mui<-mu[o]
  variancei<-variance[o]
  
  dummy<- 1:nobs
  ####a.3.smooth
  npoly<-1
  m<-n
  muout<-mui
  
  svarmiss<-lwlsld(n,bandmu,npoly,mui,variancei,m,muout)
  svar<-missvar(n,p,varlb,svarmiss)
  sigma<-svar	
  sigma[dummy[o]]<-svar	
  #maxvar<-2*var(y[,1])
  
  #sigma[sigma>maxvar]<-maxvar
  object<-list(variance=variance, sigma=sigma)
  object
}

############################################################################################
##############################
####################################################################################################
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
        }	}}
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
  
  return(you)}		 

#########################################################

##from chiou
bwgrid<-function(factor,n,nbw,gridin){
  gridout<-rep(0,nbw)
  o<-order(gridin)
  gridin<-gridin[o]
  bwmin<-gridin[1]
  bwmax<-gridin[n]
  range<-bwmax-bwmin
  
  gridout[1]<-factor*range/n
  fact<-n/factor
  fact<-fact^(1/(nbw-1))
  for(i in 2:nbw){
    gridout[i]<-gridout[i-1]*fact}
  
  return(gridout)}

###################################################################################################
##varince lower bound

varlowbd<-function(n,data,vartyp){
  if(vartyp==1){
    #		varlb=min(data[,1])
    #		if(varlb==0){varlb=0.5
    varlb<-0.01}
  if(vartyp==2){
    varlb=0.01*(1-0.01)}
  
  return(varlb)}



#####################################################################################################
##replace the zero variance

missvar<-function(n,p,varlb,variancei){
  n<-nobs
  
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
  return(variancei)}

##################################################################################################

lss.eres <- function(e,z,delta,eps,sh=FALSE)
{
  delta[e==max(e)]<-1
  nobs=length(e)
  ord <- order(e)
  ei <- e[ord]
  zi <- z[ord]
  deltai <- delta[ord]
  tie <- c(diff(ei)>eps,1)
  tie1 <- c(1,diff(ei)>eps)
  dummy <- 1:nobs
  repeats <- diff(c(dummy[c(TRUE,diff(ei)>eps)],nobs+1))
  Ni <- rev(cumsum(rev(zi)))
  di=cumsum(zi*deltai)
  di=di[tie>eps]
  di=c(di[1],diff(di))
  ieb <- 1 - di / Ni[tie1>eps]
  Shat <- cumprod(ieb)
  if(sh)
  {
    return(Shat)	
  }
  Shat <- rep(Shat,repeats)
  edif <- c(diff(ei),0)
  ehat <- rev(cumsum(rev(edif * Shat)))
  ehat[Shat<eps] <- 0
  Shat[Shat<eps] <- 1
  ehat <- ehat/Shat + ei
  eres <- ehat
  eres[dummy[ord]] <- ehat
  eres
}

lss.betag<-function(x,y,delta,z)
{
  if(is.vector(x))
    row=1
  else
    row=ncol(x)
  col=ncol(z)
  
  
  betagm<-matrix(0,ncol=col,nrow=row)
  
  ynew<-1000*(length(y))^2
  if(is.vector(x))
  {
    n1<-length(x)
    n2<- 1
  }
  else
  {
    dimnum<-dim(x)
    n1<-dimnum[1]
    n2<-dimnum[2]
  }
  
  yy0<-rep(y,rep(n1,n1))
  delta1<-rep(delta,rep(n1,n1))
  yy1<-rep(y,n1)
  yy2<-delta1*(yy0-yy1)
  
  xx0<-matrix(rep(as.vector(x),rep(n1,n1*n2)),nrow=n1*n1)
  xx1<-t(matrix(rep(as.vector(t(x)),n1),nrow=n2))
  xx2<-xx0-xx1
  
  for(i in 1:col)
  {
    zz=rep(z[,i],rep(n1,n1))*rep(z[,i],n1)
    xxdif<-xx2*zz*delta1
    xnew<-apply(xxdif,2,sum)
    xnew<-rbind(xxdif,-xnew)
    yynew<-c(yy2*zz,ynew)
    
    if(is.na(LETTERS[c(NA,2)][1])) # if running in R
      fit <- rq(yynew ~ xnew - 1, tau = 0.5)
    else
      fit <- l1fit(xnew, yynew, intercept=FALSE)
    betagm[,i] <- fit$coef
    
  }
  
  betagm
}


weightlss<-function(y,censor, x,  maxiter)
  
{
  tolerance=0.001
  
  
  y <- cbind(y,censor)
  
  x <- as.matrix(x)
  
  if(all(x[, 1] == 1))
    x <- x[, -1]	
  
  if(ncol(as.matrix(y)) != 2)
    stop("Response must be a right-censored survival object!")
  
  nobs <- nrow(y)
  if(is.vector(x))
    nvar <- 1
  else
    nvar <- ncol(x)
  
  zdummy <- matrix(rep(1,nobs), ncol=1)
  
  initialbeta <- lss(zdummy,x, y, maxiter)
  beta <- initialbeta$beta
  sigma<-rep(1,ndata)
  sigma <- lss.var(beta,cbind(1,x),y,zdummy,delta,0,sigma)$sigma
  
  niter=0
  ynew<-y
  newbeta0<-0
  betaflag<-matrix(0,nrow=maxiter,ncol=(nvar+1))
  
  while(niter < maxiter)
  {
    niter <- niter + 1
    print('niter')
    print(niter)
    betaprev <- beta
    
    if(nlevels(as.factor(sigma))>2){
      
      xnew <- cbind(1/sqrt(sigma), x/sqrt(sigma))
      ynew[,1] <- y[,1]/sqrt(sigma)
      nvarnew <- ncol(xnew)}
    else 	{nvarnew <- nvar
    xnew <- x
    ynew <- y}
    
    beta<-lssnew(zdummy,xnew,ynew,maxiter)$beta
    if(nvarnew==nvar){
      newbeta0<-0
      betaflag[niter,2:(nvar+1)]<-beta}
    else{ newbeta0<-beta[1]
    beta<-beta[-1]
    betaflag[niter,]<-beta}
    
    sigmaout <- lss.var(beta, cbind(1,x), y, zdummy, delta, 0,sigma)
    
    sigma <- sigmaout$sigma
    
    #obvar<- sigmaout$variance
    print('beta')
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
  mu<-cbind(1,x)%*%beta
  e<-y[,1]-cbind(1,x)%*%beta
  eresvar<-lss.eres(e, zdummy, y[,2],eps)
  variance<- (y[,2]*(y[,1]-mu)^2+(1-y[,2])*(eresvar)^2)
  obvar<-variance
  object<-list(beta=beta, mu=mu, sigma=sigma, obvar=obvar)
  object
}


lss<-function(zdummy,x, y, maxiter)
{
  tolerance=0.001
  
  if(all(x[, 1] == 1))
    x <- x[, -1]	
  
  
  if(ncol(as.matrix(y)) != 2)
    stop("Response must be a right-censored survival object!")
  
  nobs <- nrow(y)
  if(is.vector(x))
    nvar <- 1
  else
    nvar <- ncol(x)
  #xy<-cbind(y,x)
  #xysub<-xy[xy[,2]==1,]
  #beta<-lm(xysub[,1]~xysub[,3:(nvar+2)])$coeff[2:(nvar+1)]
  beta <- lss.betag(x, y[,1], y[,2], zdummy)
  #print('betag')
  #print(beta)
  
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
    #print('bjbeta')
    #print(c(beta0,beta))
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

lssnew<-function(zdummy,x, y, maxiter)
{
  tolerance=0.001
  
  if(all(x[, 1] == 1))
    x <- x[, -1]	
  
  
  if(ncol(as.matrix(y)) != 2)
    stop("Response must be a right-censored survival object!")
  
  nobs <- nrow(y)
  if(is.vector(x))
    nvar <- 1
  else
    nvar <- ncol(x)
  xy<-cbind(y,x)
  xysub<-xy[xy[,2]==1,]
  #beta<-lm(xysub[,1]~xysub[,3:(nvar+2)])$coeff[2:(nvar+1)]
  beta <- lss.betag(x, y[,1], y[,2], zdummy)
  #print('betag')
  #print(beta)
  
  niter=0
  if(is.vector(x))
    xbar <- mean(x)
  else
    xbar <- apply(x, 2, mean)
  xm <- x - rep(xbar, rep(nobs, nvar))
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
    #print('bjbeta')
    #print(c(beta0,beta))
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

###############################################################

rgev<-function(n, loc = 0, scale = 1, shape = 0)
{
  if(min(scale) < 0) stop("invalid scale")
  if(length(shape) != 1) stop("invalid shape")
  if(shape == 0) return(loc - scale * log(rexp(n)))
  else return(loc + scale * (rexp(n)^(-shape) - 1)/shape)
}
##################################################################################
qgev<-function(p, loc = 0, scale = 1, shape = 0, lower.tail = TRUE)
{
  if(min(p, na.rm = TRUE) <= 0 || max(p, na.rm = TRUE) >=1)
    stop("`p' must contain probabilities in (0,1)")
  if(min(scale) < 0) stop("invalid scale")
  if(length(shape) != 1) stop("invalid shape")
  if(!lower.tail) p <- 1 - p
  if(shape == 0) return(loc - scale * log(-log(p)))
  else return(loc + scale * ((-log(p))^(-shape) - 1)/shape)
}
#######################################################################
##################################################################################################3
gevMoments =
  function(xi = 0, mu = 0, beta = 1)
  {   
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Compute true statistics for Generalized Extreme Value distribution
    
    # Value:
    #   Returns true mean for xi < 1 and variance for xi < 1/2
    #   of GEV distribution, otherwise NaN is returned
    
    # FUNCTION:
    
    # MEAN: Returns for x >= 1 NaN:
    g = c(1, 0, NaN)
    xinv = 1/ ( xi + sign(abs(xi)) - 1 )
    
    # For xi = the result is eulers constant
    euler = 0.57721566490153286060651209008240243104
    xi0 = c(0, beta*euler, 0)
    
    # Supress warning for NaN's from Gamma Function:
    warn.old <- getOption("warn")
    options(warn = -1)
    gevMean = mu + beta * xinv * (gamma(1-xi)-1) * g[sign(xi-1)+2] +
      xi0[(sign(xi)+2)]
    options(warn = warn.old)
    
    # VAR: Returns for x >= 1 NaN:
    g = c(1, 0, NaN)
    xinv = 1/ ( xi + sign(abs(xi)) - 1 )
    xi0 = c(0, (beta*pi)^2 / 6, 0)
    
    # Supress warning for NaN's from Gamma Function:
    warn.old <- getOption("warn")
    options(warn = -1)
    gevVar = (beta*xinv)^2 * (gamma(1-2*xi) - gamma(1-xi)^2 ) *
      g[sign(2*xi-1)+2] + xi0[(sign(xi)+2)]
    options(warn = warn.old)
    
    # Result:
    param = c(xi = xi, mu = mu, beta = beta)
    ans = list(param = param, mean = gevMean, var = gevVar)
    
    # Return Value:
    ans
  }


###############################################################################################################################bias and var

biva<-function(beta,lbeta,true,nrun){
  bias<-rep(0,lbeta)
  var<-rep(0,lbeta)
  for(j in 1:lbeta){
    bias[j]<-0
    var[j]<-0
    for(i in 1:nrun){
      bias[j]<-bias[j]+(beta[i,j]-true[j])
      var[j]<-var[j]+(beta[i,j]-true[j])^2}
    bias[j]<-bias[j]/nrun
    var[j]<-sqrt(var[j]/nrun)}
  object<-list(bias=bias, var=var)
  object}


#return(bias,var)} 
#################################################################################### 95% coverage prob
cp<-function(beta,lbeta,truebeta,var){
  count<-rep(0,lbeta)
  percent<-rep(0,lbeta)
  for(j in 1:lbeta){
    ccount<-rep(0,nrun)
    ubound<-beta[,j]+1.96*sqrt(var[,j])
    lbound<-beta[,j]-1.96*sqrt(var[,j])
    for(i in 1:nrun){
      #print(i)
      if(truebeta[j]<ubound[i] & truebeta[j]>lbound[i]){count[j]<-count[j]+1
      
      ccount[i]<-1}}
    percent[j]<-count[j]/nrun}
  return(percent)}


########################################################################################################################################################################################################
library(survival)
library(quantreg)
library(lss2)

varlb<-0.01
set.seed(1)
eps <- .Machine$double.eps^(2/3)

###########################################################################
data("cancer", package = "survival")
attach(veteran)

veteran$celllarge=ifelse(veteran$celltype=='large',1,0)
veteran$celladeno=ifelse(veteran$celltype=='adeno',1,0)
veteran$cellsmall=ifelse(veteran$celltype=='smallcell',1,0)
veteran$cellsquamous=ifelse(veteran$celltype=='squamous',1,0)


View(veteran)



data=cbind(log(veteran$time),veteran$status)
x=cbind(1,veteran$trt, veteran$celllarge, veteran$celladeno, veteran$cellsmall, veteran$karno, veteran$diagtime, veteran$age, veteran$prior)
data=cbind(data,x)
############################################################################

ndata<-dim(data)[1]
tn<-ndata
p<-dim(data)[2]-2

y<-data[,1]
censor<-data[,2]

nboot<-100
A.boot<-matrix(0,nrow=nboot,ncol=(dim(x)[2]))

#########################################################################
nobs<-ndata
#bandmuvec<-c(ndata^(-1/5))
ncvgroup<- floor(ndata/10)
###################### for cv
wholedata<-data
mse<-rep(0,ncvgroup)
tmse<-rep(0,9)
for(j in 1:9){
  print('j')
  print(j)
  for(i in 1:10){
    print(c('i',i))
    data<-wholedata
    bandmu<-ndata^(-1/5)
    
    validata<-data[(((i-1)*ncvgroup+1):min(tn,(i*ncvgroup))),]
    traindata<-data[-(((i-1)*ncvgroup+1):min(tn,(i*ncvgroup))),]
    
    data<-traindata
    ndata<-dim(data)[1]
    nobs<-ndata
    #########################for cv end
    
    out<-weightlss(data[,1],data[,2],data[,3:(p+2)],3) 
    beta<-out$beta
    sigma<-out$sigma
    mu<-out$mu
    obvar<-out$obvar
    plot(mu,obvar, xlab='mean',ylab='observed variance',ylim=c(0,4))
    points(mu,sigma,pch=2)
    
    ##########################for cv
    nvalid<-dim(validata)[1]
    zdummy <- matrix(rep(1,nobs), ncol=1)
    
    mse[i]<-sum(validata[,2]*(validata[,1]-validata[,3:(p+2)] %*% beta)^2)
  }
  tmse[j]<-sum(mse)
}
###########################for cv end

for(ii in 1:nboot){
  print('boot')
  print(ii)
  seed<-10*ii+106
  position<-floor(runif(ndata,0,1)*ndata)+1
  datanew<-data[position,]
  
  out.boot<-weightlss(datanew[,1],datanew[,2],datanew[,3:(p+2)] , 50)
  A.boot[ii, ]<-out.boot$beta
  
}

B.boot<-apply(A.boot,2,var)




A.boot
beta
B.boot

se=sqrt(B.boot)
se
ci_l=beta + qnorm(0.025)*se
ci_l
ci_u=beta + qnorm(0.975)*se
ci_u
