### Guillaume Evin
### 22/10/2019, Grenoble
###  IRSTEA
### guillaume.evin@irstea.fr
### 
### apply QE-ANOVA Bayesian to climate projections with possible missing RCM/GCM couples, using smoothing splines
###
### REFERENCES
###
### Evin, G., B. Hingray, J. Blanchet, N. Eckert, S. Morin, and D. Verfaillie.
### Partitioning Uncertainty Components of an Incomplete Ensemble of Climate Projections Using Data Augmentation,
### Journal of Climate. https://doi.org/10.1175/JCLI-D-18-0606.1.
###
### Cheng, Chin-I., and Paul L. Speckman. 2012. Bayesian Smoothing Spline Analysis of Variance. Computational Statistics 
### & Data Analysis 56 (12): https://doi.org/10.1016/j.csda.2012.05.020.


#==============================================================================
#' QUALYPSOSS.check.option
#'
#' Check if input options provided in \code{\link{QUALYPSOSS}} are valid and assigned default values if missing.
#'
#' @param listOption list of options
#'
#' @return List containing the complete set of options.
#'
#' @author Guillaume Evin
QUALYPSOSS.check.option = function(listOption){
  if(is.null(listOption)){
    listOption = list()
  }
  
  # lambdaClimateResponse
  if('lambdaClimateResponse' %in% names(listOption)){
    lambdaClimateResponse = listOption[['lambdaClimateResponse']]
    if(!(is.numeric(lambdaClimateResponse)&lambdaClimateResponse>0)) stop('lambdaClimateResponse must be numeric and positive')
  }else{
    listOption[['lambdaClimateResponse']] = 0.1
  }
  
  # lambdaHyperParANOVA
  if('lambdaHyperParANOVA' %in% names(listOption)){
    lambdaHyperParANOVA = listOption[['lambdaHyperParANOVA']]
    if(!(is.numeric(lambdaHyperParANOVA)&lambdaHyperParANOVA>0)) stop('lambdaHyperParANOVA must be numeric and positive')
  }else{
    listOption[['lambdaHyperParANOVA']] = 0.000001
  }
  
  # typeChangeVariable
  if('typeChangeVariable' %in% names(listOption)){
    typeChangeVariable = listOption[['typeChangeVariable']]
    if(!(typeChangeVariable%in%c('abs','rel'))) stop("typeChangeVariable must be equal to 'abs' or 'rel'")
  }else{
    listOption[['typeChangeVariable']] = 'abs'
  }
  
  # nBurn
  if('nBurn' %in% names(listOption)){
    nBurn = listOption[['nBurn']]
    if(!(is.numeric(nBurn)&(nBurn>=0))) stop('wrong value for nBurn')
  }else{
    listOption[['nBurn']] = 1000
  }
  
  # nKeep
  if('nKeep' %in% names(listOption)){
    nKeep = listOption[['nKeep']]
    if(!(is.numeric(nKeep)&(nKeep>=0))) stop('wrong value for nKeep')
  }else{
    listOption[['nKeep']] = 2000
  }
  
  # nMCMC
  listOption$nMCMC = listOption$nKeep+listOption$nBurn
  listOption$vecKeepMCMC = (listOption$nBurn+1):listOption$nMCMC
  
  # nCluster
  if('nCluster' %in% names(listOption)){
    nCluster = listOption[['nCluster']]
    if(!(is.numeric(nCluster)&(nCluster>=0))) stop('wrong value for nCluster')
  }else{
    listOption[['nCluster']] = 1
  }
  
  # quantileCompress
  if('quantileCompress' %in% names(listOption)){
    quantileCompress = listOption[['quantileCompress']]
    if(!(is.numeric(quantileCompress))) stop('wrong value for quantileCompress')
  }else{
    listOption[['quantileCompress']] = c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)
  }
  
  # uniqueFit
  if('uniqueFit' %in% names(listOption)){
    uniqueFit = listOption[['uniqueFit']]
    if(!(is.logical(uniqueFit))) stop('wrong value for uniqueFit')
  }else{
    listOption[['uniqueFit']] = TRUE
  }
  
  # returnMCMC
  if('returnMCMC' %in% names(listOption)){
    returnMCMC = listOption[['returnMCMC']]
    if(!(is.logical(returnMCMC))) stop('wrong value for returnMCMC')
  }else{
    listOption[['returnMCMC']] = FALSE
  }
  
  # returnOnlyCR
  if('returnOnlyCR' %in% names(listOption)){
    returnOnlyCR = listOption[['returnOnlyCR']]
    if(!(is.logical(returnOnlyCR))) stop('wrong value for returnOnlyCR')
  }else{
    listOption[['returnOnlyCR']] = FALSE
  }
  
  # Version
  listOption[['version']] = 'v1.0.0'
  
  return(listOption)
}


#==============================================================================
#' QUALYPSOSS.process.scenario
#'
#' compute different objects used for the application of Smoothing-Splines ANOVA (SS-ANOVA),
#' these objects being processed outputs of the scenario characteristics
#'
#' @param scenAvail matrix of scenario characteristics \code{nS} x \code{nK}.
#' @param predContUnique (optional) unique values of continuous predictors.
#' 
#' @return list containing diverse information aboutwith the following fields:
#' \itemize{
#'   \item \strong{scenAvail}: Record first argument of the function,
#'   \item \strong{predContUnique}: Record second argument of the function,  
#'   \item \strong{XFull}: data.frame with all possible combinations of predictors (continuous AND discrete),
#'   \item \strong{nFull}: number of rows of \code{XFull},
#'   \item \strong{nK}: Number of columns of \code{ScenAvail} (i.e. number of discrete predictors),   
#'   \item \strong{predDiscreteUnique}: List containing possible values for each discrete predictor.
#' }
#' 
#' @author Guillaume Evin
QUALYPSOSS.process.scenario = function(scenAvail,predContUnique){
  #==== dimensions
  nK = ncol(scenAvail)
  nS = nrow(scenAvail)
  
  #==== unique predictors
  predDiscreteUnique = list()
  if(is.null(colnames(scenAvail))){
    predNames = 1:nK
  }else{
    predNames = colnames(scenAvail)
  }
  for(iCol in 1:nK) predDiscreteUnique[[predNames[iCol]]] = as.character(unique(scenAvail[,iCol]))
  
  # we expand over all possible values of the continuous predictor
  Xunique = predDiscreteUnique
  Xunique$PredCont = predContUnique
  
  #==== all possible combinations of discrete predictors
  XFull = expand.grid(Xunique)
  nFull = nrow(XFull)
  
  
  #==== scenario characeristics
  lScen = list(XFull=XFull,nFull=nFull,nK=nK,predContUnique=predContUnique,
               predDiscreteUnique=predDiscreteUnique,scenAvail=scenAvail)
}


#==============================================================================
# Bernoulli polynomials for reproducing kernels
B2scaled = function(x){
  return((x^2-x+1/6)/2)
}
B4scaled = function(x){
  return((x^4-2*x^3+x^2-1/30)/24)
}


#==============================================================================
#' reproducing.kernel
#'
#' see par 2.3 in Cheng and Speckman
#' 
#' @param x vector of predictors (continuous or discrete) 
#' @param y vector of predictors (continuous or discrete) 
#' @param type 'continuous' or 'discrete'
#' @param typeRK type of reproducing kernels: c('Cheng','Gu','Gaussian')
#' 
#' @return matrix n x n 
#' 
#' @author Guillaume Evin
reproducing.kernel = function(x,y=NULL,type,typeRK="Cheng"){
  # "Cheng": par2.3: Cheng, Chin-I., and Paul L. Speckman. 2016. Bayes Factors for Smoothing Spline ANOVA.
  # Bayesian Analysis. https://doi.org/10.1214/15-BA974.
  
  # "Gu": Eq. 4.2. Gu, Chong, and Grace Wahba. 1993. Smoothing Spline ANOVA with Component-Wise Bayesian Confidence Intervals 
  # Journal of Computational and Graphical Statistics. https://doi.org/10.2307/1390957.
  
  # square reproducing kernel
  if(is.null(y)) y=x
  
  # length of the vectors
  m = length(x)
  n = length(y)
  
  if(type=='continuous'){
    # normalize in [0,1]
    if(typeRK%in%c('Cheng','Gu')){
      min.xy = min(c(x,y))
      max.xy = max(c(x,y))
      range.xy = max.xy-min.xy
      x = (x-min.xy)/range.xy
      y = (y-min.xy)/range.xy
    }else if(typeRK=='Gaussian'){
      x = (x-mean(x))/sd(x)
      y = (y-mean(y))/sd(y)
    }
  }else if(type=='discrete'){
    # number of effects for this predictor
    nEffect = length(unique(c(x,y)))
  }else{stop('unknown type')} 
  
  # initialise matrix of reproducing kernel
  RK = matrix(nrow=m,ncol=n)
  
  for(i in 1:m){
    pi = x[i]
    
    if(type=='continuous'){
      if(typeRK=='Cheng'){
        RK[i,] = pmin(pi,y)^2*((3*pmax(pi,y)-pmin(pi,y))/6)
      }else if(typeRK=='Gaussian'){
        RK[i,] = exp(-(y-pi)^2/2)
      }else if(typeRK=='Gu'){
        RK[i,] = B2scaled(pi)*B2scaled(y)-B4scaled(abs(y-pi))
      }
    }else if(type=='discrete'){
      isEq = y==pi
      RK[i,isEq] = (nEffect-1)/nEffect
      RK[i,!isEq] = -1/nEffect
    }else{stop('unknown type')}
  }
  
  return(RK)
}


#==============================================================================
#' get.spectral.decomp
#'
#' compute different objects used for the application of Smoothing-Splines ANOVA (SS-ANOVA)
#'
#' @param SIGMA reproducing kernel
#'
#' @return list with the following fields:
#' \itemize{
#'   \item \strong{Q}: Matrix of eigen vectors n x r,
#'   \item \strong{D}: Vector of nonzero eigen values (size r),
#'   \item \strong{r}: Number of nonzero eigen values (scalar).
#' }
#'
#' @author Guillaume Evin
get.spectral.decomp = function(SIGMA){
  # decompose reproducing kernel matrix with a spectral decomposition
  eigen.out = eigen(SIGMA)
  
  # nonzero eigen values
  zz = eigen.out$values>0
  
  # eigen vectors
  Q = eigen.out$vectors[,zz]
  
  # eigen values
  D = eigen.out$values[zz]
  
  # number of nonzero eigen values
  r = length(D)
  
  # return list with all objects
  return(list(Q=Q, D=D, r=r))
}


#==============================================================================
#  extract.climate.response
#'
#' Extract climate response for one time series z
#' @param ClimateProjections matrix of climate projections
#' @param predCont matrix of continuous predictor corresponding to the climate projections
#' @param predContUnique vector of predictors for which we need fitted climate reponses
#' @param nMCMC number of MCMC samples
#' @param lam fixed smoothing parameter lambda
#' @param uniqueFit logical value indicating if only one fit is applied
#' @param parSmooth smoothing parameter \code{spar} in \code{\link[stats]{smooth.spline}}: varies in [0,1]
#' @param listCR list of objects for the extraction of the climate response
#'
#' @return list with the following fields:
#' \itemize{
#'   \item \strong{phi}: MCMC draws of climate response
#'   \item \strong{eta}: MCMC draws of deviation from the climate response
#'   \item \strong{sigma2}: MCMC draws of sigma2
#'   \item \strong{beta}: MCMC draws of beta
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin
extract.climate.response = function(ClimateProjections,predCont,predContUnique,nMCMC,lam,uniqueFit,parSmooth=1,listCR=NULL){
  # dimensions
  nT = nrow(ClimateProjections)
  nS = ncol(ClimateProjections)
  nP = length(predContUnique)
  
  # prepare outputs
  sigma2.MCMC = matrix(nrow=nMCMC,ncol=nS)
  phi.MCMC = array(dim=c(nMCMC,nS,nP))
  eta.MCMC = array(dim=c(nMCMC,nS,nT))
  
  
  # return some objects for faster computation on grids (there are the)
  if(is.null(listCR)){
    loadListCR = FALSE
    listCR = list()
  }else{
    loadListCR = TRUE
  }
  
  for(iS in 1:nS){
    # extract response and predictor for this scenario
    zRaw = ClimateProjections[,iS]
    predContRaw = predCont[,iS]
    # non missing responses
    nz = !is.na(zRaw)&!is.na(predContRaw)
    # select non missing predictors and responses
    z = zRaw[nz]
    predContS = predContRaw[nz]
    n = length(z)
    
    # avoid scaling issues in predictors
    min.xy = min(c(predContS,predContUnique))
    max.xy = max(c(predContS,predContUnique))
    range.xy = max.xy-min.xy
    predContS.scale = (predContS-min.xy)/range.xy
    predContUnique.scale = (predContUnique-min.xy)/range.xy
    
    if(uniqueFit){
      #=================== CUBIC SMOOTH SPLINE ===================
      # find predictor
      smooth.spline.out<-stats::smooth.spline(predContS,z,spar=parSmooth)
      #  fitted available responses
      phi = smooth.spline.out$y
      # fitted responses at unnown points
      phiNP = predict(smooth.spline.out, predContUnique)$y
      # residuals for available responses
      eta = z - phi
      
      #===== Store results
      sigma2.MCMC[,iS] = rep(var(eta),nMCMC)
      phi.MCMC[,iS,] = matrix(nrow=nMCMC,ncol=nP,rep(phiNP,nMCMC),byrow = T)
      eta.MCMC[,iS,nz] = matrix(nrow=nMCMC,ncol=sum(nz),rep(eta,nMCMC),byrow = T)
    }else{
      #=================== HYPARAMETERS + MCMC ===================
      if(loadListCR){
        X = listCR[[iS]]$X
        XXinvX = listCR[[iS]]$XXinvX
        XXinv = listCR[[iS]]$XXinv
        matPredNP = listCR[[iS]]$matPredNP
        XNP = listCR[[iS]]$XNP
        Q = listCR[[iS]]$Q
        Qt = listCR[[iS]]$Qt
        D = listCR[[iS]]$D
        r = listCR[[iS]]$r
      }else{
        # linear parametric terms
        X = cbind(rep(1,sum(nz)),predContS.scale)
        XXinv = MASS::ginv(t(X)%*%X)
        XXinvX = XXinv%*%t(X)
        
        # Reproducing kernel for available responses
        # only one predictor: assume equally time spaced climate response
        RK = reproducing.kernel(x=predContS.scale,type="continuous")
        RKinv = MASS::ginv(RK)
        RK.DECOMP = get.spectral.decomp(RK)
        Q = RK.DECOMP$Q
        Qt = t(Q)
        D = RK.DECOMP$D
        r = RK.DECOMP$r
        
        # Reproducing kernel for new predictors
        RKNP = reproducing.kernel(x=predContUnique.scale,y=predContS.scale,type="continuous")
        matPredNP = RKNP%*%RKinv
        XNP = cbind(rep(1,nP),predContUnique.scale)
        
        # retur objects for faster computation on grids (there are the)
        listCR[[iS]] = list(X= X, XXinvX = XXinvX, XXinv=XXinv, matPredNP=matPredNP, XNP=XNP,
                            Q=Q, Qt=Qt, D=D, r=r)
      }
      
      
      #===================== i.MCMC = 0 ===========================
      # start with a cubic spline
      smooth.spline.out<-stats::smooth.spline(predContS,z,spar=parSmooth)
      
      #  fitted available responses
      phi = smooth.spline.out$y
      
      # sigma2: invariance prior
      sigma2 = mean(((z - phi)^2))
      
      # corresponding smoothing spline effect g
      g = phi
      
      
      #================= i.MCMC = 1..n.MCMC =======================
      for(iMCMC in 1:nMCMC){
        #===== beta
        mu.beta = XXinvX%*%(z-g)
        sig.beta = XXinv*sigma2
        beta = mvtnorm::rmvnorm(n=1, mean = mu.beta, sigma = sig.beta)
        XBeta = X%*%t(beta)
        
        #===== nu
        Cinv = 1/(1 + lam/D)
        mu = diag(Cinv)%*%(Qt%*%(z-XBeta))
        nu = mu + rnorm(n=r,sd=sqrt(sigma2*Cinv))
        
        #===== g and phi
        g = Q%*%nu
        phi = XBeta + g
        
        #===== for new predictors
        phiNP = XNP%*%t(beta) + matPredNP%*%g
        
        #===== sigma2
        diff = sum((z - phi)^2)
        nuStar = sum(nu^2/D)
        sigma2 = 1/rgamma(n=1, shape=(n+r)/2, rate=(diff/2+lam*nuStar/2))
        
        #===== Store MCMC draws
        sigma2.MCMC[iMCMC,iS] = sigma2
        phi.MCMC[iMCMC,iS,] = phiNP
        eta.MCMC[iMCMC,iS,nz] = z - phi
      }
    }
  }
  
  # return results
  return(list(phi=phi.MCMC,eta=eta.MCMC,sigma2=sigma2.MCMC,listCR=listCR))
}


#=========================================================================
#' QUALYPSOSS.get.RK
#'
#' Get reproducing kernel for each discrete predictor
#' @param X matrix of predictors
#' @param nK number of discrete predictors
#' @param nCluster number of clusters used to compute the reproducing kernels
#'
#' @return strong{RK}: list containing the reproducing kernels, obtained using spectral decomposition
#' 
#' @export
#' 
#' @author Guillaume Evin  
QUALYPSOSS.get.RK = function(X,nK,nCluster){
  L = nK+1
  
  # main reproducing kernels
  SIGMA.CONTINUOUS = reproducing.kernel(x=X$PredCont, type="continuous")
  SIGMA.PRED = list()
  for(i in 1:nK) SIGMA.PRED[[i]] = reproducing.kernel(x=X[,i], type="discrete")
  
  # predictors
  SIGMA.LIST = list()
  for(i in 1:nK) SIGMA.LIST[[i]] = SIGMA.PRED[[i]]*SIGMA.CONTINUOUS
  SIGMA.LIST[[L]] = SIGMA.CONTINUOUS
  
  # spectral decomposition using parallelisation
  cl <- parallel::makeCluster(nCluster)
  doParallel::registerDoParallel(cl)
  RK = foreach(l=1:L,.export=c("get.spectral.decomp")) %dopar% {
    return(get.spectral.decomp(SIGMA.LIST[[l]]))
  }
  parallel::stopCluster(cl)
  for(l in 1:L){
    RK[[l]]$Qt = t(RK[[l]]$Q)
  }
  
  # return
  return(RK)
}


#=========================================================================
#' QUALYPSOSS.ANOVA
#'
#' SSANOVA decomposition of the ensemble of climate change responses using a Bayesian approach.
#' The different fields of the returned list contain \code{n} samples from the posterior distributions
#' of the different inferred quantities.
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param yMCMC array \code{nMCMC} x \code{nFull} of climate change responses
#' @param RK large object containing the reproducing kernels, returned by \code{\link{QUALYPSOSS.get.RK}}
#'
#' @return list containing diverse information aboutwith the following fields:
#' \itemize{
#'   \item \strong{g}: Smooth effects \code{g}: array \code{n} x \code{nFull} x \code{L} where
#'   \code{nFull} is the number of possible combinations of predictors (discrete AND continuous),
#'   \item \strong{lambda}: Smoothing parameters: matrix \code{n} x \code{L},  
#'   \item \strong{sigma2}: Residual variance: vector of length \code{n},
#'   \item \strong{MCMC.list}: list containing previous objects, for records (according to the option \code{returnMCMC}).
#' }
#' 
#' @author Guillaume Evin  
QUALYPSOSS.ANOVA = function(lOpt, yMCMC, RK){
  #====== dimensions
  L = length(RK)
  nFull = nrow(RK[[1]]$Q)
  nMCMC = lOpt$nMCMC
  MCMC.list = list()
  
  #====== dimension of the reproducing kernels
  r = vector(length=L)
  for(l in 1:L){
    r[l] = RK[[l]]$r
  }
  
  
  #====== missing responses
  isMiss = is.na(yMCMC[1,])
  nMiss = sum(isMiss)
  
  
  #====== initialise matrix and arrays
  lambda.MCMC = matrix(nrow=nMCMC,ncol=L)
  g.MCMC = array(dim=c(nMCMC,nFull,L))
  sigma2.MCMC = vector(length=nMCMC)
  nu = nu.Cinv = nu.mu = list()
  nuStar = vector(length=L)
  g = matrix(nrow=nFull,ncol=L)
  
  # Hyperparameters for the predictors
  b = rep(lOpt$lambdaHyperParANOVA,L)
  
  
  #====== i.MCMC = 0
  # sigma2: empirical variance
  sigma2 = 1
  
  # lambda: gamma prior: Gamma(1/2,2bl)
  lam = rgamma(n = L, shape = 1/2, scale = 1+b)
  
  # nu: normal prior: N(0,sigma^2/lambda_l D_l)
  for(l in 1:L) nu[[l]] = rnorm(n=RK[[l]]$r,sd=sqrt((sigma2/lam[l])*RK[[l]]$D))
  
  # g: retrieve g
  for(l in 1:L) g[,l] = RK[[l]]$Q%*%nu[[l]]
  
  # sum of main effects
  gSum = apply(g,1,sum)
  
  # climate change response
  y = yMCMC[1,]
  y[isMiss] = gSum[isMiss] + rnorm(n=nMiss, sd = sqrt(sigma2))
  
  
  #====== i.MCMC = 1..n.MCMC
  pb <- txtProgressBar(min = 1, max = nMCMC, style = 3)
  for(iMC in 1:nMCMC){
    setTxtProgressBar(pb, iMC)
    y[!isMiss] = yMCMC[iMC,!isMiss]
    
    # nuStar = nu'%*%D^{-1}%*%nu
    for(l in 1:L) nuStar[l] = sum(nu[[l]]^2/RK[[l]]$D[[l]])
    
    # sigma2
    diff = (y - gSum)
    sigma2 = 1/rgamma(n=1, shape=(nFull+sum(r))/2, rate=(sum(diff^2)/2+sum(lam*nuStar)/2))
    
    # lambda 
    lam = rgamma(n=L, shape=(r+1)/2, scale=(2*sigma2*b)/(b*nuStar+sigma2))
    
    #===== nu: Eq. 24
    for(l in 1:L){
      nu.Cinv = 1/(1 + lam[l]/RK[[l]]$D)
      diff = as.matrix(y-rowSums(g[,-l]))
      nu.mu = nu.Cinv*(RK[[l]]$Qt%*%diff)
      nu[[l]] = nu.mu + rnorm(n=r[l])*sqrt(sigma2*nu.Cinv)                       
    }
    
    # g
    for(l in 1:L) g[,l] = RK[[l]]$Q%*%nu[[l]]
    
    # sum of main effects
    gSum = rowSums(g)
    
    # completed obs: normal conditional posterior distribution
    y[isMiss] = gSum[isMiss] + rnorm(n=nMiss, sd = sqrt(sigma2))
    
    # store values
    lambda.MCMC[iMC,] = lam
    g.MCMC[iMC,,] = g
    sigma2.MCMC[iMC] = sigma2
  }
  close(pb)
  
  # store all MCMC draws if required
  if(lOpt$returnMCMC){
    MCMC.list = list(lambda=lambda.MCMC, g=g.MCMC, sigma2=sigma2.MCMC)
  }else{
    MCMC.list = NULL
  }
  
  
  # remove burn-in period
  g.sel = g.MCMC[lOpt$vecKeepMCMC,,]
  lambda.sel = lambda.MCMC[lOpt$vecKeepMCMC,]
  sigma2.sel = sigma2.MCMC[lOpt$vecKeepMCMC]
  
  # return
  return(list(g=g.sel,lambda=lambda.sel,sigma2=sigma2.sel,MCMC.list=MCMC.list))
}


#==============================================================================
#' QUALYPSOSS
#'
#' @param ClimateProjections matrix \code{nT} x \code{nS} of climate projections where \code{nT} is the number of values for the continuous predictor
#' (years, global temperature) and \code{nS} the number of scenarios.
#' @param scenAvail matrix of scenario characteristics \code{nS} x \code{nK} where \code{nK} is the number of discrete predictors.
#' @param vecYears (optional) vector of years of length \code{nT} (by default, a vector \code{1:nT}).
#' @param predCont (optional) matrix \code{nT} x \code{nS} of continuous predictors.
#' @param predContUnique (optional) vector of length \code{nP} corresponding to the continuous predictor for which we want to obtain the prediction.
#' @param iCpredCont (optional) index in \code{1:nT} indicating the reference period (reference period) for the computation of change variables.
#' @param iCpredContUnique (optional) index in \code{1:nP} indicating the reference continuous predictor for the computation of change variables.
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{lambdaClimateResponse}: smoothing parameter > 0 for the extraction of the climate response.
#'   \item \strong{lambdaHyperParANOVA}: hyperparameter \eqn{b} for the \eqn{\lambda} parameter related to each predictor \eqn{g}.
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative).
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000). If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: number of kept samples (default: 2000). If \code{nKeep} is too small, MCMC samples might not be represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{nCluster}: number of clusters used for the parallelization (default: 1). When \code{nCluster} is greater than one, parallelization is used to
#'   apply \code{QUALYPSOSS} over multiple time steps or grid points simultaneously.
#'   \item \strong{quantileCompress}: vector of probabilities (in [0,1]) for which we compute the quantiles from the posterior distributions
#'    \code{quantileCompress = c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)} by default.
#'   \item \strong{uniqueFit}: logical, if \code{FALSE} (default), climate responses are fitted using Bayesian smoothing splines, otherwise,if \code{TRUE},
#'    a unique cubic smoothing spline is fitted for each run, using the function \link[stats]{smooth.spline}.
#'   \item \strong{returnMCMC}: logical, if \code{TRUE}, the list \code{MCMC} contains MCMC chains.
#'   \item \strong{returnOnlyCR}: logical, if \code{TRUE} (default), only Climate Responses are fitted and returned.
#' }
#' @param RK Reproducing kernels: list
#' 
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{MEAN}: list containing the mean estimate of different quantities: \code{ResidualVariability} (residual variability),
#'   \code{InternalVariability} (internal variability), \code{lambda} (smoothing parameters), \code{grandMean} (grand mean for all time
#'   steps), \code{effect} (list with one item per discrete predictor \code{i}, containing matrices \code{nT} x \code{nEffi}, 
#'   where \code{nEffi} is the number of possible values for the discrete predictor \code{i}).
#'   \item \strong{QUANT}: list containing quantiles of different estimated quantities, listed in \strong{MEAN}.
#'   \item \strong{DECOMPVAR}: list with the contribution of all components to the total uncertainty, provided in \code{TotalVar} for 
#'   all time steps. In addition, for each discrete predictor, \code{ContribEffect} provides the relative contribution of possible
#'   discrete value (e.g. the contribution of one RCM to the uncertainty due to RCMs).
#'   \item \strong{MCMC.list}: list containing the MCMC chains (not returned by default).
#'   \item \strong{climateResponse}: list containing different objects related to the extraction of the climate response.
#'   phiStar (\eqn{\phi^*}) is an array \code{nQ} x \code{nS} x \code{nP} containing climate change responses, where \code{nQ} is the
#'   number of returned quantiles, \code{nS} is the number of scenarios and \code{nP} is the length of \code{predContUnique} (e.g. number
#'   of future years). 
#'   Similarly, etaStar (\eqn{\eta^*}) contains the deviation from the climate change response. 
#'   phi (\eqn{\phi}) contains the climate responses and eta (\eqn{\eta}) contains the deviations from the climate responses.
#'   \item \strong{listCR}: list containing objects created during the extraction of the climate responses (to be used as an 
#'   argument in \code{\link{QUALYPSOSSlight}})
#'   \item \strong{ClimateProjections}: argument of the call to the function, for records.
#'   \item \strong{predCont}: (optional) argument of the call to the function, for records.
#'   \item \strong{predContUnique}:  (optional) argument of the call to the function, for records.
#'   \item \strong{predDiscreteUnique}: list of possible values taken by the discrete predictors given in \code{scenAvail}.
#'   \item \strong{listOption}: list of options
#'   \item \strong{listScenario}: list of scenario characteristics (obtained from \code{\link{QUALYPSOSS.process.scenario}})
#'   \item \strong{RK}: list containing the reproducing kernels (to be used as an argument in \code{\link{QUALYPSOSSlight}})
#' }
#'
#' @examples
#' ##########################################################################
#' # SYNTHETIC SCENARIOS
#' ##########################################################################
#' # create nS=3 fictive climate scenarios with 2 GCMs and 2 RCMs, for a period of nY=20 years
#' n=20
#' t=1:n/n
#'
#' # GCM effects (sums to 0 for each t)
#' effGCM1 = t*2
#' effGCM2 = t*-2
#'
#' # RCM effects (sums to 0 for each t)
#' effRCM1 = t*1
#' effRCM2 = t*-1
#'
#' # These climate scenarios are a sum of effects and a random gaussian noise
#' scenGCM1RCM1 = effGCM1 + effRCM1 + rnorm(n=n,sd=0.5)
#' scenGCM1RCM2 = effGCM1 + effRCM2 + rnorm(n=n,sd=0.5)
#' scenGCM2RCM1 = effGCM2 + effRCM1 + rnorm(n=n,sd=0.5)
#' ClimateProjections = cbind(scenGCM1RCM1,scenGCM1RCM2,scenGCM2RCM1)
#'
#' # Here, scenAvail indicates that the first scenario is obtained with the combination of the
#' # GCM "GCM1" and RCM "RCM1", the second scenario is obtained with the combination of
#' # the GCM "GCM1" and RCM "RCM2" and the third scenario is obtained with the combination
#' # of the GCM "GCM2" and RCM "RCM1".
#' scenAvail = data.frame(GCM=c('GCM1','GCM1','GCM2'),RCM=c('RCM1','RCM2','RCM1'))
#' 
#' listOption = list(nBurn=20,nKeep=30)
#' QUALYPSOSSOUT = QUALYPSOSS(ClimateProjections=ClimateProjections,scenAvail=scenAvail,
#' listOption=listOption)
#'
#' # QUALYPSOSSOUT output contains many different information about climate projections uncertainties,
#' # which can be plotted using the following functions.
#' 
#' # plotQUALYPSOSSClimateResponse draws the climate responses, for all simulation chains, 
#' # in comparison to the raw climate responses.
#' plotQUALYPSOSSClimateResponse(QUALYPSOSSOUT)
#' 
#' # plotQUALYPSOSSClimateChangeResponse draws the climate change responses, for all simulation chains.
#' plotQUALYPSOSSClimateChangeResponse(QUALYPSOSSOUT)
#' 
#' # plotQUALYPSOSSeffect draws the estimated effects, for a discrete predictor specified by iEff,
#' # as a function of the continuous predictor.
#' plotQUALYPSOSSeffect(QUALYPSOSSOUT, iEff = 1)
#' plotQUALYPSOSSeffect(QUALYPSOSSOUT, iEff = 2)
#' 
#' # plotQUALYPSOSSgrandmean draws the estimated grand mean, as a function of the continuous predictor.
#' plotQUALYPSOSSgrandmean(QUALYPSOSSOUT)
#' 
#' # plotQUALYPSOSSTotalVarianceByScenario draws the total uncertainty and the mean effect, 
#' # for one discrete predictor, usually a RCP scenario (e.g. it provides an illustration of the
#' # future evolution and associated uncertainties for one RCP scenario).
#' plotQUALYPSOSSTotalVarianceByScenario(QUALYPSOSSOUT,nameScenario = "GCM1",iEff = 1)
#' 
#' # plotQUALYPSOSSTotalVarianceDecomposition draws the decomposition of the total variance responses,
#' # as a function of the continuous predictor.
#' plotQUALYPSOSSTotalVarianceDecomposition(QUALYPSOSSOUT)
#' 
#' @export
#'
#' @author Guillaume Evin
QUALYPSOSS = function(ClimateProjections,scenAvail,vecYears=NULL,predCont=NULL,predContUnique=NULL,iCpredCont=NULL,iCpredContUnique=NULL,listOption=NULL,RK=NULL){
  # check options
  lOpt = QUALYPSOSS.check.option(listOption)
  
  # Dimensions
  nT = nrow(ClimateProjections)
  nS = ncol(ClimateProjections)
  nK = ncol(scenAvail)
  
  # check dimensions
  if(nS!=nrow(scenAvail)) stop("The number of columns of 'ClimateProjections' must match the number of rows of 'scenAvail'")
  
  # vecYears
  if(!is.null(vecYears)){
    if(!all(is.numeric(vecYears))) stop('vecYears must be a vector of years')
  }else{
    vecYears = 1:nT
  }
  
  # predContList
  if(!is.null(predCont)){
    # check lengths
    if(any(dim(ClimateProjections)!=dim(predCont))) stop('QUALYPSOSS: ClimateProjections and predCont must have the same dimensions')
  }else{
    predCont = replicate(n=nS,vecYears)
  }
  
  # unique continuous predictor
  if(is.null(predContUnique)){
    nP = 100
    predContUnique = seq(from=min(predCont),to=max(predCont),length.out=nP)
  }else{
    nP = length(predContUnique)
  }
  
  # iCpredContUnique in predContUnique
  if(!is.null(iCpredContUnique)){
    if(!any(iCpredContUnique==1:nP)) stop('iCpredContUnique must be an index in vector predContUnique')
  }else{
    iCpredContUnique = 1
  }
  
  # iCpredCont in predCont
  if(!is.null(iCpredCont)){
    if(!any(iCpredCont==1:nT)) stop('iCpredCont must be an index indicating a reference row in predCont')
  }else{
    iCpredCont = 1
  }
  
  # process scenarios
  lScen = QUALYPSOSS.process.scenario(scenAvail,predContUnique)
  
  
  
  #=================== EXTRACT CLIMATE CHANGE RESPONSE ===================
  climResponse = extract.climate.response(ClimateProjections=ClimateProjections,
                                          predCont=predCont,
                                          predContUnique=predContUnique,
                                          nMCMC = lOpt$nMCMC,
                                          lam=lOpt$lambdaClimateResponse,
                                          uniqueFit = lOpt$uniqueFit)
  phi = climResponse$phi
  eta = climResponse$eta
  
  
  
  #=================== COMPUTE CHANGE VARIABLES ===================
  # ANOVA is applied from the reference year until the end (assume no nas)
  trimEta = iCpredCont:nT
  nEta = length(trimEta)
  
  # compute absolute or relative climate change responses for all
  # the nMCMC fitted climate responses
  typeChangeVariable = lOpt$typeChangeVariable
  phiStar = array(dim=c(lOpt$nMCMC,nS,nP))
  etaStar = array(dim=c(lOpt$nMCMC,nS,nEta))
  
  # loop over the runs
  for(i in 1:nS){
    if(typeChangeVariable=='abs'){
      phiStar[,i,] = phi[,i,] - phi[,i,iCpredContUnique]
      etaStar[,i,] = eta[,i,trimEta]
    }else if(typeChangeVariable=='rel'){
      phiStar[,i,] = phi[,i,]/phi[,i,iCpredContUnique]-1
      etaStar[,i,] = eta[,i,trimEta]/phi[,i,iCpredContUnique]
    }
  }
  
  # if relative changes are computed on zeros, projections are meaningless
  if(any(is.nan(phiStar))) return(NULL)
  
  # Climate Response
  phiQ = apply(phi[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress)
  etaQ = apply(eta[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress,na.rm=TRUE)
  phiStarQ = apply(phiStar[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress)
  etaStarQ = apply(etaStar[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress,na.rm=TRUE)
  climateResponse=list(phiStar=phiStarQ,etaStar=etaStarQ,phi=phiQ,eta=etaQ)
  
  if(lOpt$returnOnlyCR){
    return(list(climateResponse=climateResponse,ClimateProjections=ClimateProjections,
                predCont=predCont,listOption=lOpt,
                predContUnique=predContUnique,listScenario=lScen))
  }
  
  #=================== ANOVA ===================
  
  # matrix of predictors, such that columns of yMCMC corresponds
  # rows of lScen$Xfull
  yMCMC = matrix(nrow=lOpt$nMCMC,ncol=lScen$nFull)
  
  # collapse scenario characteristics
  vXFull <- apply(lScen$XFull[,1:nK], 1, paste, collapse='.')
  vscenAvail <- apply(scenAvail, 1, paste, collapse='.')
  
  for(i in 1:nS){
    zz = which(vXFull==vscenAvail[i])
    phiStarI = phiStar[,i,]
    yMCMC[,zz] = phiStarI
  }
  
  # reproducing kernels (computationally intensive)
  returnRK = is.null(RK)
  if(is.null(RK)) RK = QUALYPSOSS.get.RK(X=lScen$XFull,nK=nK,nCluster=lOpt$nCluster)
  
  # Apply smoothing-spline ANOVA model
  QUALYPSOSS.ANOVA.OUT = QUALYPSOSS.ANOVA(lOpt, yMCMC, RK)
  g = QUALYPSOSS.ANOVA.OUT$g
  lambda = QUALYPSOSS.ANOVA.OUT$lambda
  
  # store all MCMC draws if required
  if(lOpt$returnMCMC){
    MCMC.list = list(ANOVA = QUALYPSOSS.ANOVA.OUT$MCMC.list, 
                     CLIMATERESPONSE = list(phi=phi,eta=eta,phiStar=phiStar,etaStar=etaStar,sigma2=climResponse$sigma2))
  }else{
    MCMC.list = NULL
  }
  
  #===================  PROCESS RESULTS   =================== 
  # quantiles from the posterior distributions
  qq = lOpt$quantileCompress
  nqq = length(qq)
  
  # continuous predictors
  XPred = lScen$XFull
  
  # ¨discrete predictors
  predDiscreteUnique = lScen$predDiscreteUnique
  nD = unlist(lapply(predDiscreteUnique,length))
  
  #### initialize lists
  MEAN = QUANT = DECOMPVAR = list()
  
  #### Residual variability
  DECOMPVAR$ResidualVariability = mean(QUALYPSOSS.ANOVA.OUT$sigma2)
  MEAN$ResidualVariability = mean(QUALYPSOSS.ANOVA.OUT$sigma2)
  QUANT$ResidualVariability = quantile(QUALYPSOSS.ANOVA.OUT$sigma2,probs = qq)
  
  #### Internal variability
  InternalVariability = apply(apply(etaStar[lOpt$vecKeepMCMC,,],c(1,2),var,na.rm=TRUE),1,mean)
  DECOMPVAR$InternalVariability = mean(InternalVariability)
  MEAN$InternalVariability = mean(InternalVariability)
  QUANT$InternalVariability = quantile(InternalVariability,probs = qq)
  
  #### lambda
  MEAN$lambda = colMeans(lambda)
  QUANT$lambda = apply(lambda,2,quantile,probs = qq)
  
  
  #### Main change for the continuous predictor
  MEAN$grandMean = vector(length=nP)
  QUANT$grandMean = matrix(nrow=nP,ncol=nqq)
  MEAN$effect = QUANT$effect = MEAN$effectWithMean = MEAN$effectWithMean = list()
  # DECOMPVAR$effect is the variance related to each effect
  DECOMPVAR$effect = matrix(nrow=nP,ncol=nK)
  # DECOMPVAR$ContribEffect is the relative (%) contribution of each type of effect to the total variance
  DECOMPVAR$ContribEffect = list()
  
  # index continuous predictor in RK and g
  iPredCont = nK + 1
  for(i in 1:nP){
    ii = which(XPred$PredCont==predContUnique[i])
    MEAN$grandMean[i] = mean(g[,ii,iPredCont])
    QUANT$grandMean[i,] = quantile(g[,ii,iPredCont],probs=qq)
    
    for(iPred in 1:nK) DECOMPVAR$effect[i,iPred] = mean(g[,ii,iPred]^2)
  }
  DECOMPVAR$TotalVar = rowSums(DECOMPVAR$effect) + DECOMPVAR$ResidualVariability + DECOMPVAR$InternalVariability
  
  #### changes for the effects
  for(iPred in 1:nK){
    nEff = nD[iPred]
    DECOMPVAR$ContribEffect[[iPred]] = matrix(nrow=nP,ncol=nEff)
    MEAN$effect[[iPred]] = matrix(nrow=nP,ncol=nEff)
    QUANT$effect[[iPred]] = array(dim=c(nP,nEff,nqq))
    MEAN$effectWithMean[[iPred]] = matrix(nrow=nP,ncol=nEff)
    QUANT$effectWithMean[[iPred]] = array(dim=c(nP,nEff,nqq))
    
    for(iCont in 1:nP){
      for(iEff in 1:nEff){
        ii = which(XPred$PredCont==predContUnique[iCont]&XPred[,iPred]==predDiscreteUnique[[iPred]][iEff])
        MEAN$effect[[iPred]][iCont,iEff] = mean(g[,ii,iPred])
        QUANT$effect[[iPred]][iCont,iEff,] = quantile(g[,ii,iPred],probs=qq)
        MEAN$effectWithMean[[iPred]][iCont,iEff] = mean(g[,ii,iPred]+g[,ii,iPredCont])
        QUANT$effectWithMean[[iPred]][iCont,iEff,] = quantile(g[,ii,iPred]+g[,ii,iPredCont],probs=qq)
        DECOMPVAR$ContribEffect[[iPred]][iCont,iEff] = mean(g[,ii,iPred]^2)/(DECOMPVAR$effect[iCont,iPred]*nEff)
      }
    }
  }
  
  
  
  
  #===================  return results   =================== 
  
  # if returnRK is FALSE, we do not return RK (large object)
  if(!returnRK){
    RK = NULL
  }
  
  # return list of results
  return(list(MEAN=MEAN,QUANT=QUANT,DECOMPVAR=DECOMPVAR,MCMC.list=MCMC.list,
              climateResponse=climateResponse,listCR = climResponse$listCR,
              ClimateProjections=ClimateProjections,predCont=predCont,
              predContUnique=predContUnique,predDiscreteUnique=predDiscreteUnique,
              listOption=lOpt,listScenario=lScen,RK=RK))
}


#==============================================================================
#' plotQUALYPSOSSClimateResponse
#'
#' Plot climate responses.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param lim y-axis limits (default is NULL)
#' @param col color for the lines
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSClimateResponse = function(QUALYPSOSSOUT,lim=NULL,col=NULL,
                                       xlab="Years",ylab=expression(phi),...){
  # Continuous predictor
  predCont = QUALYPSOSSOUT$predCont
  predContUnique = QUALYPSOSSOUT$predContUnique
  
  # Available climate responses
  ClimateProj = QUALYPSOSSOUT$ClimateProjections
  
  # Available scenarios
  scenAvail = QUALYPSOSSOUT$listScenario$scenAvail
  nS = nrow(scenAvail)
  
  # climate responses
  iMedian = which(QUALYPSOSSOUT$listOption$quantileCompress==0.5)
  phi = QUALYPSOSSOUT$climateResponse$phi[iMedian,,]
  
  # process options
  if(is.null(lim)) lim = range(phi)
  if(is.null(col)) col = 1:nS
  
  # Figure
  plot(-100,-100,xlim=range(predContUnique),ylim=lim,xlab=xlab,ylab=ylab,...)
  for(i in 1:nS){
    # raw projections
    lines(predCont[,i],ClimateProj[,i],lwd=2,lty=2,col=i)
    # fitted climate response
    lines(predContUnique,phi[i,],col=col[i],lwd=2)
  }
}


#==============================================================================
#' plotQUALYPSOSSClimateChangeResponse
#'
#' Plot climate change responses.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param lim y-axis limits (default is NULL)
#' @param col color for the lines
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSClimateChangeResponse = function(QUALYPSOSSOUT,lim=NULL,col=NULL,
                                             xlab="Years",ylab=expression(phi^{star}),...){
  # Continuous predictor
  predContUnique = QUALYPSOSSOUT$predContUnique
  
  # Available scenarios
  scenAvail = QUALYPSOSSOUT$listScenario$scenAvail
  nS = nrow(scenAvail)
  
  # climate responses
  iMedian = which(QUALYPSOSSOUT$listOption$quantileCompress==0.5)
  phiStar = QUALYPSOSSOUT$climateResponse$phiStar[iMedian,,]
  
  # process options
  if(is.null(lim)) lim = range(phiStar)
  if(is.null(col)) col = 1:nS
  
  # Figure
  plot(-100,-100,xlim=range(predContUnique),ylim=lim,xlab=xlab,ylab=ylab,...)
  for(i in 1:nS){
    # fitted climate response
    lines(predContUnique,phiStar[i,],col=col[i],lwd=2)
  }
}


#==============================================================================
# plotQUALYPSOSSgetCI
#
# return confidence level given indices in QUALYPSOSSOUT$listOption$quantileCompress
plotQUALYPSOSSgetCI = function(QUALYPSOSSOUT,iBinf,iBsup){
  vecq = QUALYPSOSSOUT$listOption$quantileCompress
  return(round((vecq[iBsup]-vecq[iBinf])*100))
}

#==============================================================================
#' plotQUALYPSOSSgrandmean
#'
#' Plot prediction of grand mean ensemble. By default, we plot the credible interval corresponding to a probability 0.95.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param CIlevel probabilities for the credible intervals, default is equal to \code{c(0.025,0.975)}
#' @param lim y-axis limits (default is NULL)
#' @param col color for the overall mean and the credible interval
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSgrandmean = function(QUALYPSOSSOUT,CIlevel=c(0.025,0.975),lim=NULL,
                                 col='black',xlab="Continuous predictor",ylab="Grand mean",addLegend=T,
                                 ...){
  # Continuous predictor
  TTmix = QUALYPSOSSOUT$predContUnique
  Iord = order(TTmix)
  nTT = length(TTmix)
  TT = TTmix[Iord]
  
  # find index quantiles
  iMedian = which(QUALYPSOSSOUT$listOption$quantileCompress==0.5)
  iBinf = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[1])
  iBsup = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[2])
  
  # retrieve median and limits
  med = QUALYPSOSSOUT$QUANT$grandMean[,iMedian][Iord]
  binf = QUALYPSOSSOUT$QUANT$grandMean[,iBinf][Iord]
  bsup = QUALYPSOSSOUT$QUANT$grandMean[,iBsup][Iord]
  
  # colors polygon
  colPoly = adjustcolor(col,alpha.f=0.2)
  
  # initiate plot
  if(is.null(lim)) lim = range(c(binf,bsup),na.rm=TRUE)
  plot(-100,-100,xlim=range(TT),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)
  
  # add confidence interval
  polygon(c(TT,rev(TT)),c(binf,rev(bsup)),col=colPoly,lty=0)
  
  # add median
  lines(TT,med,lwd=3,col=col)
  
  # legend
  if(addLegend){
    pctCI = plotQUALYPSOSSgetCI(QUALYPSOSSOUT,iBinf,iBsup)
    legend('topleft',bty='n',fill=c(NA,colPoly),lwd=c(2,NA),lty=c(1,NA),
           border=c(NA,col),col=c(col,NA),legend=c('Median',paste0(pctCI,'%CI')))
  }
}


#==============================================================================
#' plotQUALYPSOSSeffect
#'
#' Plot prediction of ANOVA effects for one main effect. By default, we plot we plot the credible intervals corresponding to a probability 0.95.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param iEff index of the main effect to be plotted in \code{QUALYPSOSSOUT$listScenario$predDiscreteUnique}
#' @param includeMean if TRUE, the grand mean is added to the main effect in the plot
#' @param CIlevel probabilities for the credible intervals, default is equal to \code{c(0.025,0.975)}
#' @param lim y-axis limits (default is NULL)
#' @param col colors for each effect
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSeffect = function(QUALYPSOSSOUT,iEff,includeMean=FALSE,CIlevel=c(0.025,0.975),lim=NULL,
                              col=1:20,xlab="Continuous predictor",ylab="Effect",addLegend=TRUE,
                              ...){
  # Continuous predictor
  TTmix = QUALYPSOSSOUT$predContUnique
  Iord = order(TTmix)
  nTT = length(TTmix)
  TT = TTmix[Iord]
  
  # retrieve effects
  if(includeMean){
    QUANTEff = QUALYPSOSSOUT$QUANT$effectWithMean[[iEff]]
  }else{
    QUANTEff = QUALYPSOSSOUT$QUANT$effect[[iEff]]
  }
  nEff = dim(QUANTEff)[2]
  
  # find index quantiles
  iMedian = which(QUALYPSOSSOUT$listOption$quantileCompress==0.5)
  iBinf = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[1])
  iBsup = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[2])
  
  # retrieve median and limits
  med = QUANTEff[,,iMedian]
  binf = QUANTEff[,,iBinf]
  bsup = QUANTEff[,,iBsup]
  
  # initiate plot
  if(is.null(lim)) lim = range(c(binf,bsup),na.rm=TRUE)
  plot(-100,-100,xlim=range(TT),ylim=c(lim[1],lim[2]),xlab=xlab,ylab=ylab,...)
  
  for(i in 1:nEff){
    # colors polygon
    colPoly = adjustcolor(col[i],alpha.f=0.2)
    
    # add confidence interval
    polygon(c(TT,rev(TT)),c(binf[,i][Iord],rev(bsup[,i][Iord])),col=colPoly,lty=0)
    
    # add median
    lines(TT,med[,i][Iord],lwd=3,col=col[i])
  }
  
  # legend
  if(addLegend){
    pctCI = plotQUALYPSOSSgetCI(QUALYPSOSSOUT,iBinf,iBsup)
    legend('topleft',bty='n',fill=c(NA,'black'),lwd=c(2,NA),lty=c(1,NA),
           border=c(NA,'black'),col=c('black',NA),legend=c('Median',paste0(pctCI,'%CI')))
    
    legend('bottomleft',bty='n',lwd=2,lty=1,col=col,
           legend=QUALYPSOSSOUT$listScenario$predDiscreteUnique[[iEff]])
  }
}


#==============================================================================
#' plotQUALYPSOSSTotalVarianceDecomposition
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param col colors for each source of uncertainty, the first two colors corresponding to internal variability and residual variability, respectively
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSTotalVarianceDecomposition = function(QUALYPSOSSOUT,
                                                  col=c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1"),
                                                  xlab="Continuous predictor",ylab="% Total Variance",addLegend=TRUE,...){
  # Continuous predictor
  TTmix = QUALYPSOSSOUT$predContUnique
  Iord = order(TTmix)
  nTT = length(TTmix)
  TT = TTmix[Iord]
  
  # Variance decomposition
  VV = QUALYPSOSSOUT$DECOMPVAR
  
  # number of main effects
  nEff = ncol(VV$effect)
  
  # initialise matrix variances
  vNorm = matrix(nrow=nTT,ncol=nEff+2)
  
  # part of variance due to the different effects
  for(iEff in 1: nEff){
    vNorm[,iEff] = VV$eff[Iord,iEff]/VV$TotalVar[Iord]
  }
  
  # Residual variance
  vNorm[,nEff+1] = VV$ResidualVariability/VV$TotalVar[Iord]
  
  # Internal variability
  vNorm[,nEff+2] = VV$InternalVariability/VV$TotalVar[Iord]
  
  
  # figure
  col = col[1:(nEff+2)]
  cum=rep(0,nTT)
  plot(-1,-1,xlim=range(TT),ylim=c(0,1),xaxs="i",yaxs="i",las=1,
       xlab=xlab,ylab=ylab,...)
  for(i in 1:(nEff+2)){
    cumPrevious = cum
    cum = cum + vNorm[,i]
    polygon(c(TT,rev(TT)),c(cumPrevious,rev(cum)),col=rev(col)[i],lty=1)
  }
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)
  
  # legend
  if(addLegend){
    if(is.null(colnames(QUALYPSOSSOUT$listScenario$scenAvail))){
      namesEff = paste0("Eff",1:nEff)
    }else{
      namesEff = colnames(QUALYPSOSSOUT$listScenario$scenAvail)
    }
    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}


#==============================================================================
#' plotQUALYPSOSSTotalVarianceByScenario
#'
#' Plot fraction of total variance explained by each source of uncertainty.
#'
#' @param QUALYPSOSSOUT output from \code{\link{QUALYPSOSS}}
#' @param iEff index in \code{scenAvail} corresponding to the scenarios (e.g. RCP scenarios)
#' @param nameScenario name of the scenario to be plotted (as provided in \code{scenAvail})
#' @param probCI probability for the dredible interval, =0.9 by default
#' @param col colors for each source of uncertainty, the first two colors corresponding to internal variability and residual variability, respectively
#' @param ylim y-axis limits
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param addLegend if TRUE, a legend is added
#' @param ... additional arguments to be passed to \code{\link[graphics]{plot}}
#'
#' @export
#'
#' @author Guillaume Evin
plotQUALYPSOSSTotalVarianceByScenario = function(QUALYPSOSSOUT,iEff,nameScenario,probCI=0.9,col=NULL,ylim=NULL,
                                               xlab="Years",ylab="Change variable",addLegend=TRUE,...){
  # number of years
  TT = QUALYPSOSSOUT$predContUnique
  nTT = length(TT)
  
  # which scenario
  iScenario = which(QUALYPSOSSOUT$listScenario$predDiscreteUnique[[iEff]] == nameScenario)
  
  # mean prediction
  meanPred = QUALYPSOSSOUT$MEAN$effectWithMean[[iEff]][,iScenario]
  
  # Variance decomposition
  VV = QUALYPSOSSOUT$DECOMPVAR
  
  # remove effect corresponding to the scenarios from the total variance
  Veff = VV$effect[,-iEff]
  
  # concatenate variances
  Vbind = cbind(Veff,rep(VV$ResidualVariability,nTT),rep(VV$InternalVariability,nTT))
  nEff = ncol(Vbind)-2
  Vtot = rowSums(Vbind)
  Vnorm = Vbind/replicate(n = ncol(Vbind), Vtot)
  
  # reverse
  vNormRev = t(apply(t(Vnorm),2,rev))
  
  
  # compute the lower bound if the distribution is gaussian
  binf = qnorm(p = (1-probCI)/2, mean = meanPred, sd = sqrt(VV$TotalVar))
  bsup = qnorm(p = 0.5+probCI/2, mean = meanPred, sd = sqrt(VV$TotalVar))
  
  # figure
  if(is.null(col)){
    default.col = c("orange","yellow","cadetblue1","blue1","darkgreen","darkgoldenrod4","darkorchid1")
    col = default.col[1:(nEff+2)]
  }
  
  # obtain limits of the intervals, proportion corresponds to the part of the variance, lower and upper than the mean
  limIntInf = limIntSup = matrix(nrow=nTT,ncol=nEff+2)
  limIntInf[,1] = binf
  limIntSup[,1] = bsup
  for(i in 1:(nEff+1)){
    limIntInf[,i+1] = limIntInf[,i]+vNormRev[,i]*(meanPred-binf)
    limIntSup[,i+1] = limIntSup[,i]-vNormRev[,i]*(bsup-meanPred)
  }
  
  # figure
  if(is.null(ylim)) ylim = c(min(binf),max(bsup))
  plot(-1,-1,xlim=range(TT),ylim=ylim,xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",las=1,...)
  for(i in 1:(nEff+2)){
    polygon(c(TT,rev(TT)),c(limIntInf[,i],rev(limIntSup[,i])),col=col[i],lty=1)
  }
  lines(TT,meanPred,col="white",lwd=1)
  
  # add horizontal lines
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)
  
  # legend
  if(addLegend){
    if(is.null(colnames(QUALYPSOSSOUT$listScenario$scenAvail))){
      namesEff = paste0("Eff",1:nEff)
    }else{
      namesEff = colnames(QUALYPSOSSOUT$listScenario$scenAvail)[-iEff]
    }
    
    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}



#==============================================================================
#' QUALYPSOSSlight
#'
#' same as QUALYPSOSS, but less outputs are returned, and arguments are mandatory, in order to limit processing tasks which are repeated 
#' over a grid.
#'
#' @param ClimateProjections matrix \code{nT} x \code{nS} of climate projections where \code{nT} is the number of values for the continuous predictor
#' (years, global temperature) and \code{nS} the number of scenarios.
#' @param scenAvail matrix of scenario characteristics \code{nS} x \code{nK} where \code{nK} is the number of discrete predictors.
#' @param predCont (optional) matrix \code{nT} x \code{nS} of continuous predictors.
#' @param predContUnique (optional) vector of length \code{nP} corresponding to the continuous predictor for which we want to obtain the prediction.
#' @param iCpredCont (optional) index in \code{1:nT} indicating the reference period (reference period) for the computation of change variables.
#' @param iCpredContUnique (optional) index in \code{1:nP} indicating the reference continuous predictor for the computation of change variables.
#' @param listOption (optional) list of options
#' \itemize{
#'   \item \strong{lambdaClimateResponse}: smoothing parameter > 0 for the extraction of the climate response.
#'   \item \strong{lambdaHyperParANOVA}: hyperparameter \eqn{b} for the \eqn{\lambda} parameter related to each predictor \eqn{g}.
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative).
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000). If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: number of kept samples (default: 2000). If \code{nKeep} is too small, MCMC samples might not be represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{nCluster}: number of clusters used for the computation of reproducing kernels (default: 1). 
#'   When \code{nCluster} is greater than one, parallelization is used to apply \code{QUALYPSOSS} over multiple time steps or 
#'   grid points simultaneously.
#'   \item \strong{quantileCompress}: vector of probabilities (in [0,1]) for which we compute the quantiles from the posterior distributions
#'    \code{quantileCompress = c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)} by default.
#'    \code{uniqueFit}: logical, if \code{FALSE} (default), climate responses are fitted using Bayesian smoothing splines, otherwise,if \code{TRUE},
#'    a unique cubic smoothing spline is fitted for each run, using the function \link[stats]{smooth.spline}.
#'    \code{returnMCMC}: logical, if \code{FALSE} (default), the list \code{MCMC} is empty in the returned object.
#' }
#' @param lScen list of objects related to the scenario characteristics: item of the list obtained from \code{\link{QUALYPSOSS}}
#' @param RK Reproducing kernels: item of the list obtained from \code{\link{QUALYPSOSS}}
#' @param listCR Object for the extraction of the climate response: item of the list obtained from \code{\link{QUALYPSOSS}}
#'
#' @return list with the following fields:
#'
#' \itemize{
#'   \item \strong{MEAN}: list containing the mean estimate of different quantities: \code{ResidualVariability} (residual variability),
#'   \code{InternalVariability} (internal variability), \code{lambda} (smoothing parameters), \code{grandMean} (grand mean for all time
#'   steps), \code{effect} (list with one item per discrete predictor \code{i}, containing matrices \code{nT} x \code{nEffi}, 
#'   where \code{nEffi} is the number of possible values for the discrete predictor \code{i}).
#'   \item \strong{QUANT}: list containing quantiles of different estimated quantities, listed in \strong{MEAN}.
#'   \item \strong{DECOMPVAR}: list with the contribution of all components to the total uncertainty, provided in \code{TotalVar} for 
#'   all time steps. In addition, for each discrete predictor, \code{ContribEffect} provides the relative contribution of possible
#'   discrete value (e.g. the contribution of one RCM to the uncertainty due to RCMs).
#' }
#' 
#' @export
#'
#' @author Guillaume Evin
QUALYPSOSSlight = function(ClimateProjections,scenAvail,predCont,predContUnique,iCpredCont,iCpredContUnique,listOption,lScen,RK,listCR){
  # check options
  lOpt = QUALYPSOSS.check.option(listOption)
  
  # Dimensions
  nT = nrow(ClimateProjections)
  nS = ncol(ClimateProjections)
  nK = ncol(scenAvail)
  nP = length(predContUnique)
  
  #=================== EXTRACT CLIMATE CHANGE RESPONSE ===================
  climResponse = extract.climate.response(ClimateProjections=ClimateProjections,
                                          predCont=predCont,
                                          predContUnique=predContUnique,
                                          nMCMC = lOpt$nMCMC,
                                          lam=lOpt$lambdaClimateResponse,
                                          uniqueFit = lOpt$uniqueFit,
                                          listCR=listCR)
  phi = climResponse$phi
  eta = climResponse$eta
  
  
  
  #=================== COMPUTE CHANGE VARIABLES ===================
  # ANOVA is applied from the reference year until the end (assume no nas)
  trimEta = iCpredCont:nT
  nEta = length(trimEta)
  
  # compute absolute or relative climate change responses for all
  # the nMCMC fitted climate responses
  typeChangeVariable = lOpt$typeChangeVariable
  phiStar = array(dim=c(lOpt$nMCMC,nS,nP))
  etaStar = array(dim=c(lOpt$nMCMC,nS,nEta))
  
  # loop over the runs
  for(i in 1:nS){
    if(typeChangeVariable=='abs'){
      phiStar[,i,] = phi[,i,] - phi[,i,iCpredContUnique]
      etaStar[,i,] = eta[,i,trimEta]
    }else if(typeChangeVariable=='rel'){
      phiStar[,i,] = phi[,i,]/phi[,i,iCpredContUnique]-1
      etaStar[,i,] = eta[,i,trimEta]/phi[,i,iCpredContUnique]
    }
  }
  
  # if relative changes are computed on zeros, projections are meaningless
  if(any(is.nan(phiStar))) return(NULL)
  
  
  #=================== ANOVA ===================
  
  # matrix of predictors, such that columns of yMCMC corresponds
  # rows of lScen$Xfull
  yMCMC = matrix(nrow=lOpt$nMCMC,ncol=lScen$nFull)
  
  # collapse scenario characteristics
  vXFull <- apply(lScen$XFull[,1:nK], 1, paste, collapse='.')
  vscenAvail <- apply(scenAvail, 1, paste, collapse='.')
  
  for(i in 1:nS){
    zz = which(vXFull==vscenAvail[i])
    phiStarI = phiStar[,i,]
    yMCMC[,zz] = phiStarI
  }
  
  
  # Apply smoothing-spline ANOVA model
  QUALYPSOSS.ANOVA.OUT = QUALYPSOSS.ANOVA(lOpt, yMCMC, RK)
  g = QUALYPSOSS.ANOVA.OUT$g
  lambda = QUALYPSOSS.ANOVA.OUT$lambda
  
  
  
  #===================  PROCESS RESULTS   =================== 
  # quantiles from the posterior distributions
  qq = lOpt$quantileCompress
  nqq = length(qq)
  
  # continuous predictors
  XPred = lScen$XFull
  predContUnique = lScen$predContUnique
  nP = length(predContUnique)
  
  # ¨discrete predictors
  predDiscreteUnique = lScen$predDiscreteUnique
  nD = unlist(lapply(predDiscreteUnique,length))
  
  #### initialize lists
  MEAN = QUANT = DECOMPVAR = list()
  
  #### Residual variability
  DECOMPVAR$ResidualVariability = mean(QUALYPSOSS.ANOVA.OUT$sigma2)
  MEAN$ResidualVariability = mean(QUALYPSOSS.ANOVA.OUT$sigma2)
  QUANT$ResidualVariability = quantile(QUALYPSOSS.ANOVA.OUT$sigma2,probs = qq)
  
  #### Internal variability
  InternalVariability = apply(apply(etaStar[lOpt$vecKeepMCMC,,],c(1,2),var,na.rm=TRUE),1,mean)
  DECOMPVAR$InternalVariability = mean(InternalVariability)
  MEAN$InternalVariability = mean(InternalVariability)
  QUANT$InternalVariability = quantile(InternalVariability,probs = qq)
  
  #### lambda
  MEAN$lambda = colMeans(lambda)
  QUANT$lambda = apply(lambda,2,quantile,probs = qq)
  
  
  #### Main change for the continuous predictor
  MEAN$grandMean = vector(length=nP)
  QUANT$grandMean = matrix(nrow=nP,ncol=nqq)
  MEAN$effect = QUANT$effect = MEAN$effectWithMean = MEAN$effectWithMean = list()
  # DECOMPVAR$effect is the variance related to each effect
  DECOMPVAR$effect = matrix(nrow=nP,ncol=nK)
  # DECOMPVAR$ContribEffect is the relative (%) contribution of each type of effect to the total variance
  DECOMPVAR$ContribEffect = list()
  
  # index continuous predictor in RK and g
  iPredCont = nK + 1
  for(i in 1:nP){
    ii = which(XPred$PredCont==predContUnique[i])
    MEAN$grandMean[i] = mean(g[,ii,iPredCont])
    QUANT$grandMean[i,] = quantile(g[,ii,iPredCont],probs=qq)
    
    for(iPred in 1:nK) DECOMPVAR$effect[i,iPred] = mean(g[,ii,iPred]^2)
  }
  DECOMPVAR$TotalVar = rowSums(DECOMPVAR$effect) + DECOMPVAR$ResidualVariability + DECOMPVAR$InternalVariability
  
  #### changes for the effects
  for(iPred in 1:nK){
    nEff = nD[iPred]
    DECOMPVAR$ContribEffect[[iPred]] = matrix(nrow=nP,ncol=nEff)
    MEAN$effect[[iPred]] = matrix(nrow=nP,ncol=nEff)
    QUANT$effect[[iPred]] = array(dim=c(nP,nEff,nqq))
    MEAN$effectWithMean[[iPred]] = matrix(nrow=nP,ncol=nEff)
    QUANT$effectWithMean[[iPred]] = array(dim=c(nP,nEff,nqq))
    
    for(iCont in 1:nP){
      for(iEff in 1:nEff){
        ii = which(XPred$PredCont==predContUnique[iCont]&XPred[,iPred]==predDiscreteUnique[[iPred]][iEff])
        MEAN$effect[[iPred]][iCont,iEff] = mean(g[,ii,iPred])
        QUANT$effect[[iPred]][iCont,iEff,] = quantile(g[,ii,iPred],probs=qq)
        MEAN$effectWithMean[[iPred]][iCont,iEff] = mean(g[,ii,iPred]+g[,ii,iPredCont])
        QUANT$effectWithMean[[iPred]][iCont,iEff,] = quantile(g[,ii,iPred]+g[,ii,iPredCont],probs=qq)
        DECOMPVAR$ContribEffect[[iPred]][iCont,iEff] = mean(g[,ii,iPred]^2)/(DECOMPVAR$effect[iCont,iPred]*nEff)
      }
    }
  }
  
  
  # Climate Response
  phiQ = apply(phi[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress)
  etaQ = apply(eta[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress,na.rm=TRUE)
  phiStarQ = apply(phiStar[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress)
  etaStarQ = apply(etaStar[lOpt$vecKeepMCMC,,],c(2,3),quantile,lOpt$quantileCompress,na.rm=TRUE)
  climateResponse=list(phiStar=phiStarQ,etaStar=etaStarQ,phi=phiQ,eta=etaQ)
  
  
  #===================  return results   =================== 
  # return list of results
  return(list(MEAN=MEAN,QUANT=QUANT,DECOMPVAR=DECOMPVAR))
}
