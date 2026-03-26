# Load libraries
#install.packages('',dependencies = TRUE,repos = 'http://cran.rstudio.com')
library(survival)
library(survminer)
library(dplyr)
library(mdhglm)
library(MASS)
library(SurvRegCensCov)
library(eha)
library(writexl)
library(ushr)
library(ggplot2)


#data
data("cancer", package = "survival")
attach(veteran)

veteran$celllarge=ifelse(veteran$celltype=='large',1,0)
veteran$celladeno=ifelse(veteran$celltype=='adeno',1,0)
veteran$cellsmall=ifelse(veteran$celltype=='smallcell',1,0)
veteran$cellsquamous=ifelse(veteran$celltype=='squamous',1,0)


####################################################################################################################

#Fit Cox model
vet.cox=coxph(Surv(time,status) ~ trt +celllarge + celladeno + cellsmall + karno + diagtime + age + prior, data=veteran)
summary(vet.cox)
vet.cox.est=vet.cox$coefficients
vet.cox.est
vet.cox.se=sqrt(diag(vet.cox$var))
vet.cox.se
ci_l=vet.cox.est+qnorm(0.025)*vet.cox.se
ci_l
ci_u=vet.cox.est+qnorm(0.975)*vet.cox.se
ci_u

###########################################################################################################################
#Test PH
vet.ph=cox.zph(vet.cox)
vet.ph

#windows.options(width=10, height=10)
par(mfrow=c(1, 1))

plot(vet.ph,var = 2)
abline(h = coef(vet.cox)[2], col = "red", lwd = 2)
plot(vet.ph,var = 3)
abline(h = coef(vet.cox)[3], col = "red", lwd = 2)
plot(vet.ph,var = 5)
abline(h = coef(vet.cox)[5], col = "red", lwd = 2)


largesurv=survfit(Surv(time, status) ~ celllarge, data = veteran)
par(mfrow=c(1,1))
plot(largesurv,fun="cloglog",
     xlab =  "log of time",
     ylab = "log-log Survival",
     col=c("blue", "red"))
legend("bottomright",legend =c("celllarge=0", "celllarge=1"), col=c("blue","red"),lwd=2 )

adenosurv=survfit(Surv(time, status) ~ celladeno, data = veteran)
par(mfrow=c(1,1))
plot(adenosurv,fun="cloglog",
     xlab =  "log of time",
     ylab = "log-log Survival",
     col=c("blue", "red"))
legend("bottomright",legend =c("celladeno=0", "celladeno=1"), col=c("blue","red"),lwd=2 )


quantile(veteran$karno)
kar=cut(veteran$karno, breaks=c(10,40,60,75,99), labels=c("kar1","kar2", "kar3", "kar4"))
karsurv=survfit(Surv(time, status) ~ kar, data = veteran)
par(mfrow=c(1,1))
plot(karsurv,fun="cloglog",
     xlab =  "log of time",
     ylab = "log-log Survival",
     col=c("blue", "red", "green","black"))
legend("bottomright",legend =c("kar1", "kar2", "kar3", "kar4"), col=c("blue","red","green","black"),lwd=4 )




######################################################################################################################################
veteran$karcategory=as.numeric(with(veteran, kar))
veteran$K= with(veteran, interaction(celllarge,celladeno, karcategory))






#####################################################################################################################################
id=1:nrow(veteran)
veteran=cbind(id,veteran)
#vet1=tmerge(veteran, veteran,id=id, endtime=event(time,status))
vet.ext=tmerge(veteran, veteran,id=id, endtime=event(time,status), kr=tdc(karno), ceel=tdc(as.numeric(celltype)))
vet.ext



######################################################AFT Weibull################################################################


aftwei=survreg(Surv(time,status)~ trt +celllarge + celladeno + cellsmall + karno + diagtime + age + prior, data=veteran, dist = "weibull")
wei.aft=ConvertWeibull(aftwei)
wei.aft
ci_l=wei.aft$vars[3:10,1]+qnorm(.025)*wei.aft$vars[3:10,2]
ci_l
ci_u=wei.aft$vars[3:10,1]+qnorm(.975)*wei.aft$vars[3:10,2]
ci_u
weiAIC=AIC(aftwei)
weiAIC
weiBIC=BIC(aftwei)
weiBIC

par(mfrow=c(1,1))
kaplan=survfit(Surv(time, status) ~ 1,data=veteran)
kaplan.surv=kaplan$surv
kaplan.time=kaplan$time
log.log.kmestimates=log(-log(kaplan.surv))
log.time=log(kaplan.time)
plot(log.log.kmestimates ~ log.time)
lm.model=lm(log.log.kmestimates[1:100] ~ log.time[1:100])
abline(lm.model)



