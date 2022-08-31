## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  warning = FALSE,
  fig.width=12, 
  fig.height=8,
  collapse = TRUE,
  comment = "#>"
)

## ----echo = F, results = 'hide'-----------------------------------------------
library(ggthemes)
col.blindfree = colorblind_pal()(8)
mycol.GCM = col.blindfree[1:4]
mycol.RCM = col.blindfree[4:8]

# export figure
export.figure = FALSE

## -----------------------------------------------------------------------------
library(QUALYPSO)

# Projections of mean winter temperature for 20 simulations of 129 years are contained in the matrix Y 
#[20 x 129]
dim(Y)

# Corresponding years are provided in the vector X_time_vec
X_time_vec

# GCM and RCM models GCM et RCM corresponding to the 20 simulations are given in scenAvail (data.frame 
# with two columns "GCM" et "RCM" and 20 lines corresponding to the 20 simulations). These 20 simulations
# correspond to 4 GCMs downscaled with 5 RCMs
dim(scenAvail)
scenAvail
apply(scenAvail,2,unique)

# We define a vector of future years to reduce the dimension of the SS-ANOVA treatment
Xfut_time = seq(from=1999,to=2099,by=5)

# This list contains the different objects processed by the ANOVA methods
lInput = list(Y=Y,vecYears=X_time_vec,Xref=1999,Xfut=Xfut_time,scenAvail=scenAvail)

# in listOption, we can specify the type of change (absolute "abs" or relative "rel")
listOption = list(typeChangeVariable="abs")

# The call to QUALYPSO extracts the climate change responses and applies a simple
# ANOVA for each future year
QUALYPSOOUT = QUALYPSO(Y=lInput$Y,X=lInput$vecYears,scenAvail=lInput$scenAvail,
                       Xref=lInput$Xref,Xfut=lInput$Xfut,
                       listOption=listOption)
                       
# phiStar is matrix 20 x 21 (n=20 simulations x nFut=21 future time steps) containing the climate change responses
phiStar = QUALYPSOOUT$CLIMATEESPONSE$phiStar

# vec.GCM contains the vector of the 4 different GCMs 
vec.GCM = unique(scenAvail$GCM)

# Figure 1
par(mar=c(4.5,4.5,0.5,0.5))
plot(-1,-1,xlim=range(lInput$Xfut),ylim=range(phiStar),cex.axis=1.3,
     xlab="Years",ylab="Temperature change [degC]",cex.lab=1.7)
for(i in 1:nrow(scenAvail)){
  i.GCM = which(scenAvail$GCM[i]==vec.GCM)
  lines(lInput$Xfut,phiStar[i,],lwd=2,col=mycol.GCM[i.GCM])
}
legend('topleft',bty='n',cex=1.3,lty=1,lwd=2,col=mycol.GCM,legend=vec.GCM)

## ----echo = T, results = 'hide'-----------------------------------------------
# we load the package qualypsoss
library(qualypsoss,exclude=c("Y","scenAvail"))

# the different options listed here specify: absolute change; 1,000 MCMC draws for the burn-in period,
#5,000 MCMC draws to represent the posterior distributions, MCMC draws are returned in the output
#(returnMCMC=TRUE), a unique climate change response is used (uniqueFit=TRUE) with the smoothing parameter
#spar=1, we consider an AR1 process for the autocorrelation of the residual errors, and a linear evolution
#of their standard deviations.
listOption = list(typeChangeVariable="abs",nBurn=1000,nKeep=5000,returnMCMC=TRUE,
                  uniqueFit=TRUE,spar=1,type.temporal.dep="AR1",type.hetero="linear")

# Call to the main function of the package qualypsoss, see ?QUALYPSOSS for further information on the
#inputs and possible options 
QUALYPSOSSOUT = QUALYPSOSS(ClimateProjections=t(lInput$Y), 
                           scenAvail=lInput$scenAvail,
                           vecYears = lInput$vecYears,
                           predContUnique = lInput$Xfut,
                           iCpredCont = which(lInput$vecYears==lInput$Xref), # 1999
                           iCpredContUnique = which(lInput$Xfut==lInput$Xref), # 1999
                           listOption=listOption)

## -----------------------------------------------------------------------------
# the functions plotQUALYPSOeffect and plotQUALYPSOSSeffect directly display the evolution of the main
#effects from the two ANOVA approaches
par(mar=c(2.5,2.5,3,0.5),mfrow=c(2,2))

plotQUALYPSOeffect(QUALYPSOOUT,nameEff = "GCM",col = mycol.GCM,lim = c(-0.6,0.6),
                   xlab="",ylab="",main="(a) ANOVA-TI / Main GCM effects",cex.lab=1.3)
grid()
plotQUALYPSOSSeffect(QUALYPSOSSOUT,iEff=1,col = mycol.GCM,lim = c(-0.6,0.6),
                     xlab="",ylab="",main="(b) SS-ANOVA / Main GCM effects",cex.lab=1.3)
grid()
plotQUALYPSOeffect(QUALYPSOOUT,nameEff = "RCM",col = mycol.RCM,lim = c(-0.6,0.55),
                   xlab="",ylab="",main="(c) ANOVA-TI / Main RCM effects",cex.lab=1.3)
grid()
plotQUALYPSOSSeffect(QUALYPSOSSOUT,iEff=2,col = mycol.RCM,lim = c(-0.6,0.55),
                     xlab="",ylab="",main="(d) SS-ANOVA / Main RCM effects",cex.lab=1.3)
grid()

## -----------------------------------------------------------------------------
# QUALYPSOSSOUT$MCMC contains the MCMC draws of the different inferred quantities. We plot here the density
#of the posterior distributions using these draws.
par(mar=c(4.5,4.5,0.5,0.5),mfrow=c(1,3))
plot(density(QUALYPSOSSOUT$MCMC$LAMBDA[,3]),main="",cex.lab=1.5,xlim=c(0,0.00035),
     xlab=expression(lambda),ylab="Posterior density",lwd=2)
lines(density(QUALYPSOSSOUT$MCMC$LAMBDA[,2]),lty=2,lwd=2)
lines(density(QUALYPSOSSOUT$MCMC$LAMBDA[,1]),lty=3,lwd=2)
legend("top",lwd=2,lty=c(2,3,1),bty="n",cex=1.5,
       legend=c(expression(lambda[1]),expression(lambda[2]),expression(lambda[3])))
plot(density(QUALYPSOSSOUT$MCMC$RESIDUALVAR),main="",cex.lab=1.5,lwd=2,
     xlab=expression(delta[RV]),ylab="Posterior density")
plot(density(QUALYPSOSSOUT$MCMC$RHO),main="",cex.lab=1.5,lwd=2,
     xlab=expression(rho),ylab="Posterior density")

## -----------------------------------------------------------------------------
# we extract here the standard deviations of the residual errors (square root of the variances)
sigRes.QUA = sqrt(QUALYPSOOUT$RESIDUALVAR$MEAN)
sigRes.SSinf = sqrt(QUALYPSOSSOUT$BAYES$RESIDUALVAR[2,])
sigRes.SSmed = sqrt(QUALYPSOSSOUT$BAYES$RESIDUALVAR[4,])
sigRes.SSsup = sqrt(QUALYPSOSSOUT$BAYES$RESIDUALVAR[6,])

par(mar=c(4.5,4.5,0.5,0.5),mfrow=c(1,1))
plot(lInput$Xfut,sigRes.QUA,ylim=c(0,0.22),type="l",lwd=2,cex.lab=1.3,
     xlab="Years",ylab="Residual var. (standard deviation)")
polygon(x = c(lInput$Xfut,rev(lInput$Xfut)),y=c(sigRes.SSinf,rev(sigRes.SSsup)),
        col=adjustcolor("black", alpha.f = 0.2), lty = 0)
lines(lInput$Xfut,sigRes.SSmed,ylim=c(0,0.15),lwd=2,lty=2)
legend("topleft", bty = "n", fill = c(NA, NA, "grey"), lwd = c(2, 2, NA), lty = c(1, 2, NA), 
       border = c(NA,NA,"black"), col = c("black","black", NA), 
       legend = c("ANOVA-TI","SS-ANOVA (Median)","SS-ANOVA (95%CI)"))

