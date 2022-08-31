### Guillaume Evin
### 01/02/2020, Grenoble
###  INRAE / ETNA
### guillaume.evin@inrae.fr
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

#================================================
#' get.matrix.KMS
#' Return the square Kac-Murdoch-Szego matrix for a rho correlation and n lines/colums
#' 
#' @param n nummber of lines/columns of the square matrix
#' @param rho correlation parameter in [0,1]
#'
#' @return n x n Kac-Murdock-Szego matrix 
#' 
#' @references Kac, M., W. L. Murdock, and G. Szego. 1953. 'On the Eigen-Values of Certain Hermitian Forms'
#' Journal of Rational Mechanics and Analysis 2: 767-800.
#' 
#' @author Guillaume Evin  
get.matrix.KMS <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) - 
                    (1:n - 1))
  rho^exponent
}


#================================================
#' get.matrix.KMSinv 
#' return the inverse of the square Kac-Murdock-Szego matrix for a rho correlation and n lines/colums
#' 
#' @param n nummber of lines/columns of the square matrix
#' @param rho correlation parameter in (-1,1)
#'
#' @return n x n Kac-Murdock-Szego matrix 
#' 
#' @references Kac, M., W. L. Murdock, and G. Szego. 1953. 'On the Eigen-Values of Certain Hermitian Forms'
#' Journal of Rational Mechanics and Analysis 2: 767-800.
#' 
#' @author Guillaume Evin 
get.matrix.KMSinv <- function(n, rho) {
  Minv = diag(rep(1+rho^2),n)
  Minv[1,1] = 1
  Minv[n,n] = 1
  
  Minv[cbind(2:n,1:(n-1))] = -rho
  Minv[cbind(1:(n-1),2:n)] = -rho
  
  Minv/(1-rho^2)
}

get.matrix.KMS.test = function(){return(get.matrix.KMS(10,0.8)%*%get.matrix.KMSinv(10,0.8) - diag(10))}


#================================================
#' get.vec.weight.hetero 
#' returns the vector of weights for the computation of heteroscedastic errors
#' corresponding to one simulation chain
#' 
#' @param nP length of the continuous predictor for which we want to obtain the prediction (e.g. time)
#' we suppose that continuous predictor is regularly spaced (e.g. 1990,2000,2010,...)
#' @param type.weight.hetero "constant" (homoscedastic) or "linear" (heteroscedastic)
#'
#' @return vector of square roots of weights of the same length than \code{predContUnique}
#' 
#' @author Guillaume Evin 
get.vec.weight.hetero = function(nP,type.weight.hetero){
  # return weights V such that W = VCV
  if(type.weight.hetero=="constant"){
    v = rep(1,nP)
  }else if(type.weight.hetero=="linear"){
    #v = c(0.01,((1:(nP-1))/(nP-1)))
    v = (1:nP/nP)
  }
  return(v)
}


#================================================
#' get.matrix.hetero 
#' returns the matrix of weights for the computation of heteroscedastic errors
#' corresponding to the entire ensemble
#' 
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#'
#' @return V matrix n x n of weights where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @author Guillaume Evin 
get.matrix.hetero = function(weight.hetero,nMO){
  W0_delta = diag(weight.hetero)
  
  # kronecker product
  IMO = diag(nMO)
  return(W0_delta %x% IMO)
}


#================================================
#' get.matrix.hetero.inv 
#' returns the inverse of the matrix of weights for the computation of heteroscedastic errors
#' corresponding to the entire ensemble
#' 
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#'
#' @return inverse matrix n x n of weights where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @author Guillaume Evin 
get.matrix.hetero.inv = function(weight.hetero,nMO){
  W0inv_delta = diag(1/weight.hetero)
  
  # kronecker product
  IMO = diag(nMO)
  return(W0inv_delta %x% IMO)
}

#WW = get.matrix.hetero(weight.hetero,nMO)
#WWinv = get.matrix.hetero.inv(weight.hetero,nMO)
#boxplot(WW%*%WWinv-diag(nFull))



#================================================
#' get.matrix.AR1 
#' return the matrix of AR(1) correlations corresponding to the entire ensemble
#' 
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#' @param nMO number of possible simulation chains (missing and non-missing)
#'
#' @return C matrix n x n of AR(1) correlations where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @author Guillaume Evin 
get.matrix.AR1 = function(nP,rho,nMO){
  W0 = get.matrix.KMS(nP,rho)
  IMO = diag(nMO)
  return(W0 %x% IMO)
}


#================================================
#' get.matrix.AR1.inv 
#' return the inverse matrix of AR(1) correlations corresponding to the entire ensemble
#' 
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#' @param nMO number of possible simulation chains (missing and non-missing)
#'
#' @return inverse matrix n x n of AR(1) correlations where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @author Guillaume Evin 
get.matrix.AR1.inv = function(nP,rho,nMO){
  W0inv = get.matrix.KMSinv(nP,rho)
  IMO = diag(nMO)
  return(W0inv %x% IMO)
}

#WW = get.matrix.AR1(nP,rho,nMO)
#WWinv = get.matrix.AR1.inv(nP,rho,nMO)
#boxplot(WW%*%WWinv-diag(nFull))


#================================================
#' get.matrix.W 
#' return the matrix of W = V x C x V for the treatment of heteroscedastic and AR(1) errors
#' see Wang (2011) section 5.3 for further details
#' 
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#'
#' @return matrix n x n where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @references Wang, Y. 2011. 'Spline Smoothing with Heteroscedastic and/or Correlated Errors.'
#' Smoothing Splines. Chapman and Hall/CRC. https://doi.org/10.1201/b10954-11.
#' 
#' @author Guillaume Evin 
get.matrix.W = function(weight.hetero,nMO,nP,rho){
  
  # weights for the heteroscedasticity
  W_delta = diag(get.matrix.hetero(weight.hetero,nMO))
  
  # weigths for the autocorrelation
  W_rho = get.matrix.AR1(nP,rho,nMO)
  
  # return product of weights
  return(t(W_delta * W_rho) * W_delta)
}


#================================================
#' get.matrix.Winv 
#' return the inverse matrix of W = V x C x V for the treatment of heteroscedastic and AR(1) errors
#' see Wang (2011) section 5.3 for further details
#' 
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#'
#' @return inverse matrix n x n of weights where code{n} is the total number of predictions 
#' (all the predictions for all the possible simulation chains)
#' 
#' @references Wang, Y. 2011. 'Spline Smoothing with Heteroscedastic and/or Correlated Errors' 
#' Smoothing Splines. Chapman and Hall/CRC. https://doi.org/10.1201/b10954-11.
#' 
#' @author Guillaume Evin 
get.matrix.Winv = function(weight.hetero,nMO,nP,rho){
  # weights for the heteroscedasticity
  Winv_delta = diag(get.matrix.hetero.inv(weight.hetero,nMO))
  
  # weigths for the autocorrelation
  Winv_rho = get.matrix.AR1.inv(nP,rho,nMO)
  
  # return product of weights
  return(t(Winv_delta * Winv_rho) * Winv_delta)
}

#================================================
#' get.det.KMS 
#' return the determinant of the KMS matrix
#' 
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#'
#' @return determinant of the KMS matrix
#' 
#' @author Guillaume Evin 
get.det.KMS <- function(nP, rho) {
  return((1-rho^2)^(nP-1))
}

#================================================
#' get.det.AR1 
#' return the determinant of the matrix provided by  \code{\link{get.matrix.AR1}}
#' 
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#' @param nMO number of possible simulation chains (missing and non-missing)
#'
#' @return determinant of the AR1 matrix
#' 
#' @author Guillaume Evin 
get.det.AR1 = function(nP,rho,nMO){
  return(get.det.KMS(nP,rho)^nMO)
}

#det.AR1.v1 = get.det.AR1(nP,rho,nMO)
#mat.AR1 = get.matrix.AR1(nP,rho,nMO)
#det.AR1.num = det(mat.AR1)


#================================================
#' get.logdet.W 
#' Return the logarithm of the determinant of the matrix W
#' 
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#'
#' @return logarithm of the determinant of the matrix W
#' 
#' @author Guillaume Evin 
get.logdet.W = function(weight.hetero,nMO,nP,rho){
  logdet.V0 = 2*sum(log(weight.hetero))
  logdet.C0 = (nP-1)*log(1-rho^2)
  return(nMO*(logdet.V0+logdet.C0))
}

# weight.hetero = get.vec.weight.hetero(nP,"linear")
# get.logdet.W(weight.hetero,nMO,nP,rho.hat0)
# W = get.matrix.W(weight.hetero,nMO,nP,rho.hat0)
# log(det(W))

#================================================
#' get.target.density.rho 
#' Return the log-density of the full conditional distribution for the parameter rho
#' 
#' @param nFull nP x nMO
#' @param deltaRV variance of the residual terms for the max value of the continuous predictor
#' @param distSS sum of square distances between the climate change responses and the ANOVA model
#' @param weight.hetero output of \code{\link{get.vec.weight.hetero}}
#' @param nMO number of possible simulation chains (missing and non-missing)
#' @param nP number of continuous predictors (e.g. future times)
#' @param rho AR(1) correlation parameter in (-1,1)
#'
#' @return log-density of the full conditional distribution
#' 
#' @author Guillaume Evin 
get.target.logdensity.rho = function(nFull,deltaRV,distSS,weight.hetero,nMO,nP,rho){
  # target log-density
  target.logdensity = - nFull*log(deltaRV)/2 - get.logdet.W(weight.hetero,nMO,nP,rho)/2 - (1/deltaRV)*distSS/2
  
  return(target.logdensity)
}


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
  
  # spar
  if('spar' %in% names(listOption)){
    spar = listOption[['spar']]
    if(!(is.numeric(spar)&spar>0)) stop('spar must be numeric and positive')
  }else{
    listOption[['spar']] = 1
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
    if(!(typeChangeVariable%in%c('abs','rel','none'))) stop("typeChangeVariable must be equal to 'abs', 'rel' or 'none'")
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
  
  # iid / AR1
  if('type.temporal.dep' %in% names(listOption)){
    type.temporal.dep = listOption[['type.temporal.dep']]
    if(!(type.temporal.dep%in%c("iid","AR1"))) stop('wrong value for type.temporal.dep')
  }else{
    listOption[['type.temporal.dep']] = "AR1"
  }
  
  # homoscedasticity / heteroscedasticity
  if('type.hetero' %in% names(listOption)){
    type.hetero = listOption[['type.hetero']]
    if(!(type.hetero%in%c("linear","constant"))) stop('wrong value for type.hetero')
  }else{
    listOption[['type.hetero']] = "linear"
  }
  
  # Version
  listOption[['version']] = 'v1.1.0'
  
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
#' @param spar smoothing parameter \code{spar} in \code{\link[stats]{smooth.spline}}: must be greater than zero
#' @param listCR list of objects for the extraction of the climate response
#'
#' @return list with the following fields:
#' \itemize{
#'   \item \strong{phi.MCMC}: MCMC draws of climate response
#'   \item \strong{eta.MCMC}: MCMC draws of deviation from the climate response
#'   \item \strong{deltaIV.MCMC}: MCMC draws of deltaRV
#'   \item \strong{listCR}: list of objects for faster computation on grids
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin
extract.climate.response = function(ClimateProjections,predCont,predContUnique,nMCMC,lam,uniqueFit,spar=spar,listCR=NULL){
  # dimensions
  nT = nrow(ClimateProjections)
  nS = ncol(ClimateProjections)
  nP = length(predContUnique)
  
  # prepare outputs
  deltaIV.MCMC = matrix(nrow=nMCMC,ncol=nS)
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
      smooth.spline.out<-stats::smooth.spline(predContS,z,spar=spar)
      #  fitted available responses
      phi = smooth.spline.out$y
      # fitted responses at unnown points
      phiNP = predict(smooth.spline.out, predContUnique)$y
      # residuals for available responses
      eta = z - phi
      
      #===== Store results
      deltaIV.MCMC[,iS] = rep(var(eta),nMCMC)
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
        
        # return objects for faster computation on grids
        listCR[[iS]] = list(X= X, XXinvX = XXinvX, XXinv=XXinv, matPredNP=matPredNP, XNP=XNP,
                            Q=Q, Qt=Qt, D=D, r=r)
      }
      
      
      #===================== i.MCMC = 0 ===========================
      # start with a cubic spline
      smooth.spline.out<-stats::smooth.spline(predContS,z,spar=spar)
      
      #  fitted available responses
      phi = smooth.spline.out$y
      
      # deltaIV: invariance prior
      deltaIV = mean(((z - phi)^2))
      
      # corresponding smoothing spline effect g
      g = phi
      
      
      #================= i.MCMC = 1..n.MCMC =======================
      for(iMCMC in 1:nMCMC){
        #===== beta
        mu.beta = XXinvX%*%(z-g)
        sig.beta = XXinv*deltaIV
        beta = mvtnorm::rmvnorm(n=1, mean = mu.beta, sigma = sig.beta)
        XBeta = X%*%t(beta)
        
        #===== nu
        Cinv = 1/(1 + lam/D)
        mu = diag(Cinv)%*%(Qt%*%(z-XBeta))
        nu = mu + rnorm(n=r,sd=sqrt(deltaIV*Cinv))
        
        #===== g and phi
        g = Q%*%nu
        phi = XBeta + g
        
        #===== for new predictors
        phiNP = XNP%*%t(beta) + matPredNP%*%g
        
        #===== deltaIV
        diff = sum((z - phi)^2)
        nuStar = sum(nu^2/D)
        deltaIV = 1/rgamma(n=1, shape=(n+r)/2, rate=(diff/2+lam*nuStar/2))
        
        #===== Store MCMC draws
        deltaIV.MCMC[iMCMC,iS] = deltaIV
        phi.MCMC[iMCMC,iS,] = phiNP
        eta.MCMC[iMCMC,iS,nz] = z - phi
      }
    }
  }
  
  # return results
  return(list(phi.MCMC=phi.MCMC,eta.MCMC=eta.MCMC,deltaIV.MCMC=deltaIV.MCMC,listCR=listCR))
}


#==============================================================================
#  compute.change.variable
#'
#' Compute change variables
#' @param climResponse output from \code{\link{extract.climate.response}}
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param iCpredContUnique index in \code{1:nP} indicating the reference continuous predictor 
#' for the computation of change variables.
#' @param iCpredCont index in \code{1:nT} indicating the reference period (reference period) 
#' for the computation of change variables.
#'
#' @return list with the following fields:
#' \itemize{
#'   \item \strong{phiStar.MCMC}: MCMC draws of climate change response
#'   \item \strong{etaStar.MCMC}: MCMC draws of deviation from the climate change response
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin
compute.change.variable = function(climResponse,lOpt,lDim,iCpredContUnique,iCpredCont){
  # retrieve phi and eta
  phi.MCMC = climResponse$phi.MCMC
  eta.MCMC = climResponse$eta.MCMC
  
  # retrieve some objects
  nS = lDim$nS
  nP = lDim$nP
  nT = lDim$nT
  
  # ANOVA is applied from the reference year until the end (assume no nas)
  trimEta = iCpredCont:nT
  nEta = length(trimEta)
  
  # compute absolute or relative climate change responses for all
  # the nMCMC fitted climate responses
  typeChangeVariable = lOpt$typeChangeVariable
  phiStar.MCMC = array(dim=c(lOpt$nMCMC,nS,nP))
  etaStar.MCMC = array(dim=c(lOpt$nMCMC,nS,nEta))
  
  # loop over the runs
  for(i in 1:nS){
    if(typeChangeVariable=='abs'){
      phiStar.MCMC[,i,] = phi.MCMC[,i,] - phi.MCMC[,i,iCpredContUnique]
      etaStar.MCMC[,i,] = eta.MCMC[,i,trimEta]
    }else if(typeChangeVariable=='rel'){
      phiStar.MCMC[,i,] = phi.MCMC[,i,]/phi.MCMC[,i,iCpredContUnique]-1
      etaStar.MCMC[,i,] = eta.MCMC[,i,trimEta]/phi.MCMC[,i,iCpredContUnique]
    }else if(typeChangeVariable=='none'){
      phiStar.MCMC[,i,] = phi.MCMC[,i,]
      etaStar.MCMC[,i,] = eta.MCMC[,i,trimEta]
    }
  }
  
  # if relative changes are computed on zeros, projections are meaningless
  if(any(is.nan(phiStar.MCMC))) return(NULL)
  
  # Climate Response
  change.variable=list(phiStar.MCMC=phiStar.MCMC,etaStar.MCMC=etaStar.MCMC)
  return(change.variable)
}


#=========================================================================
#' get.yMCMC
#'
#' Get matrix \code{nMCMC} x \code{nFull} of climate responses where \code{nMCMC}
#' is the number of MCMC draws and \code{nFull} is the number of possible combinations of 
#' predictors (discrete AND continuous),
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param lScen list of scenario characteristics, output from \code{\link{QUALYPSOSS.process.scenario}}
#' @param change.variable output from \code{\link{compute.change.variable}} containing MCMC draws 
#' of climate change response
#'
#' @return strong{yMCMC}: matrix \code{nMCMC} x \code{nFull} of climate responses
#' 
#' @export
#' 
#' @author Guillaume Evin  
get.yMCMC = function(lOpt,lDim,lScen,change.variable){
  
  # matrix of predictors, such that columns of yMCMC corresponds
  # rows of lDim$nfull
  yMCMC = matrix(nrow=lOpt$nMCMC,ncol=lDim$nFull)
  
  # collapse scenario characteristics
  vXFull <- apply(lScen$XFull[,1:lDim$nK,drop=FALSE], 1, paste, collapse='.')
  vscenAvail <- apply(lScen$scenAvail, 1, paste, collapse='.')
  vXFull.unique = unique(vXFull)
  nMO = length(vXFull.unique)
  
  for(i in 1:lDim$nS){
    zz = which(vXFull==vscenAvail[i])
    phiStarI = change.variable$phiStar.MCMC[,i,]
    yMCMC[,zz] = phiStarI
  }
  
  return(yMCMC)
}


#=========================================================================
#' QUALYPSOSS.get.RK
#'
#' Get reproducing kernel for each discrete predictor
#' @param X matrix of predictors
#' @param nK number of discrete predictors
#'
#' @return strong{RK}: list containing the reproducing kernels, obtained using spectral decomposition
#' 
#' @export
#' 
#' @author Guillaume Evin  
QUALYPSOSS.get.RK = function(X,nK){
  nE = nK+1
  
  # main reproducing kernels
  SIGMA.CONTINUOUS = reproducing.kernel(x=X$PredCont, type="continuous")
  SIGMA.PRED = list()
  for(i in 1:nK) SIGMA.PRED[[i]] = reproducing.kernel(x=X[,i], type="discrete")
  
  # predictors
  SIGMA.LIST = list()
  for(i in 1:nK) SIGMA.LIST[[i]] = SIGMA.PRED[[i]]*SIGMA.CONTINUOUS
  SIGMA.LIST[[nE]] = SIGMA.CONTINUOUS
  
  # spectral decomposition
  RK = list()
  for(iE in 1:nE){
    RK[[iE]] = get.spectral.decomp(SIGMA.LIST[[iE]])
    RK[[iE]]$Qt = t(RK[[iE]]$Q)
    RK[[iE]]$SIGMA = SIGMA.LIST[[iE]]
  }
  
  # return
  return(RK)
}


#=========================================================================
#' QUALYPSOSS.ANOVA.step1
#'
#' SSANOVA decomposition of the ensemble of climate change responses using a Bayesian approach.
#' The different fields of the returned list contain \code{n} samples from the posterior distributions
#' of the different inferred quantities. In this first step, the residual errors are assumed iid
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param yMCMC array \code{nMCMC} x \code{nFull} of climate change responses
#' @param RK large object containing the reproducing kernels, returned by \code{\link{QUALYPSOSS.get.RK}}
#'
#' @return list containing diverse information aboutwith the following fields:
#' \itemize{
#'   \item \strong{g.MCMC}: Smooth effects \code{g}: array \code{n} x \code{nFull} x \code{K} where
#'   \code{nFull} is the number of possible combinations of predictors (discrete AND continuous),
#'   \item \strong{nu.MCMC}: Smooth effects \code{nu}, a list with matrices of eigen vectors
#'   \item \strong{lambda.MCMC}: Smoothing parameters: matrix \code{n} x \code{K},  
#'   \item \strong{deltaRV.MCMC}: Residual variance: vector of length \code{n},
#'   \item \strong{g.hat}: Smooth effects estimates: matrix \code{nFull} x \code{K} where
#'   \item \strong{nu.hat}: Smooth effects estimates: a list with estimates of eigen vectors,
#'   \item \strong{lambda.hat}: Smoothing parameters estimates: vector of length \code{K},  
#'   \item \strong{deltaRV.hat}: Residual variance estimate.
#'   \item \strong{logLK}: vector of log-likelihood values of the draws
#'   \item \strong{logPost}: vector of log-posterior values of the draws
#'   \item \strong{Schwarz}: Schwarz criteria
#'   \item \strong{BIC}: BIC criteria
#'   
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin  
QUALYPSOSS.ANOVA.step1 = function(lOpt, lDim, yMCMC, RK){
  #====== dimensions
  nMCMC = lOpt$nMCMC
  MCMC.list = list()
  nE = lDim$nE
  nFull = lDim$nFull
  
  
  #====== dimension of the reproducing kernels
  r = sapply(RK, "[[", "r")
  
  
  #====== missing responses
  isMiss = is.na(yMCMC[1,])
  nMiss = sum(isMiss)
  
  
  #====== initialise matrix and arrays
  lambda.MCMC = matrix(nrow=nMCMC,ncol=nE)
  g.MCMC = array(dim=c(nMCMC,nFull,nE))
  deltaRV.MCMC = vector(length=nMCMC)
  nu = nu.Cinv = nu.mu = nu.MCMC = list()
  nuStar = vector(length=nE)
  g = matrix(nrow=nFull,ncol=nE)
  for(iE in 1:nE) nu.MCMC[[iE]] = matrix(nrow=nMCMC,ncol=RK[[iE]]$r)
  logLK = logPost = vector(length=nMCMC)
  
  # Hyperparameters for the predictors
  b = rep(lOpt$lambdaHyperParANOVA,nE)
  
  
  #====== i.MCMC = 0
  # deltaRV: empirical variance
  deltaRV = 1
  
  # lambda: gamma prior: Gamma(1/2,2bl)
  lam = rgamma(n = nE, shape = 1/2, scale = 2*b)
  
  # nu: normal prior: N(0,sigma^2/lambda_l D_l)
  for(iE in 1:nE) nu[[iE]] = rnorm(n=RK[[iE]]$r,sd=sqrt((deltaRV/lam[iE])*RK[[iE]]$D))
  
  # g: retrieve g
  for(iE in 1:nE) g[,iE] = RK[[iE]]$Q%*%nu[[iE]]
  
  # sum of main effects
  gSum = apply(g,1,sum)
  
  # climate change response
  y = yMCMC[1,]
  y[isMiss] = gSum[isMiss] + rnorm(n=nMiss, sd = sqrt(deltaRV))
  
  
  #====== i.MCMC = 1..n.MCMC
  pb <- txtProgressBar(min = 1, max = nMCMC, style = 3)
  for(iMC in 1:nMCMC){
    setTxtProgressBar(pb, iMC)
    y[!isMiss] = yMCMC[iMC,!isMiss]
    
    # nuStar = nu'%*%D^{-1}%*%nu
    for(iE in 1:nE) nuStar[iE] = sum(nu[[iE]]^2/RK[[iE]]$D[[iE]])
    
    # deltaRV
    SS = sum((y - gSum)^2) # sum of squares
    deltaRV = 1/rgamma(n=1, shape=(nFull+sum(r))/2, rate=(SS/2+sum(lam*nuStar)/2))
    
    # lambda 
    lam = rgamma(n=nE, shape=(r+1)/2, scale=(2*deltaRV*b)/(b*nuStar+deltaRV))
    
    # nu: Eq. 24
    for(iE in 1:nE){
      nu.Cinv = 1/(1 + lam[iE]/RK[[iE]]$D)
      diff = as.matrix(y-rowSums(g[,-iE,drop=FALSE]))
      nu.mu = nu.Cinv*(RK[[iE]]$Qt%*%diff)
      nu[[iE]] = nu.mu + rnorm(n=r[iE])*sqrt(deltaRV*nu.Cinv)                       
    }
    
    # g
    for(iE in 1:nE) g[,iE] = RK[[iE]]$Q%*%nu[[iE]]
    
    # sum of main effects
    gSum = rowSums(g)
    
    # completed obs: normal conditional posterior distribution
    y[isMiss] = gSum[isMiss] + rnorm(n=nMiss, sd = sqrt(deltaRV))
    
    # store values
    for(iE in 1:nE) nu.MCMC[[iE]][iMC,] = nu[[iE]]
    lambda.MCMC[iMC,] = lam
    g.MCMC[iMC,,] = g
    deltaRV.MCMC[iMC] = deltaRV
    
    
    ### Bayes factor ###
    
    # log-likelihood
    SS = sum((y - gSum)^2) # sum of squares
    logLK[iMC] = -SS/(2*deltaRV) - log(2*pi)*(nFull/2) - log(deltaRV)*(nFull/2)
    
    # log-prior (smooth part)
    logPriorSmoother_t1 = vector(length=nE)
    for(iE in 1:nE) logPriorSmoother_t1[iE] = sum(nu[[iE]]^2/RK[[iE]]$D[[iE]])*lam[iE]
    logPriorSmoother_t2 = (1/2)*(log(lam))
    for(iE in 1:nE) logPriorSmoother_t2[iE] = logPriorSmoother_t2[iE] - (1/2)*sum(log(RK[[iE]]$D[[iE]]))
    logPriorSmoother = -sum(logPriorSmoother_t1)/(2*deltaRV) -(sum(r)/2)*log(2*pi) -(sum(r)/2)*log(deltaRV) + 
      sum(logPriorSmoother_t2)
    
    # log-prior delta
    logPriorDelta = -log(deltaRV)
    
    # log-prior lambda
    logPriorLambda = -(1/2)*sum(log(lam)) - sum(lam/(2*b)) - nE*lgamma(1/2) -(1/2)*sum(log(2*b))
    
    # log-posterior
    logPost[iMC] = logLK[iMC] + logPriorSmoother + logPriorDelta + logPriorLambda
  }
  close(pb)
  
  # remove burn-in period
  zMCMC = lOpt$vecKeepMCMC
  g.sel = g.MCMC[zMCMC,,]
  lambda.sel = lambda.MCMC[zMCMC,]
  deltaRV.sel = deltaRV.MCMC[zMCMC]
  nu.sel = list()
  for(iE in 1:nE) nu.sel[[iE]] = nu.MCMC[[iE]][zMCMC,]
  logLK.sel = logLK[zMCMC]
  logPost.sel = logPost[zMCMC]
  
  # BIC & Schwarz
  Schwarz = max(logLK.sel) - (sum(r)+nE+1)*log(nFull)/2
  BIC = -2*Schwarz 
  
  # estimates
  g.hat = apply(g.sel,c(2,3),mean)
  lambda.hat = apply(lambda.sel,2,mean)
  deltaRV.hat = mean(deltaRV.sel)
  nu.hat = list()
  for(iE in 1:nE) nu.hat[[iE]] = apply(nu.sel[[iE]],2,mean)
  
  # check BIC
  std.err = (y - rowSums(g.hat))/sqrt(deltaRV.hat)
  SS = sum(std.err^2)
  check.BIC = list(y=y,gSum=rowSums(g.hat),std.err=std.err,SS=SS,
                   maxLK=max(logLK.sel),nTheta=sum(r)+nE+1,
                   deltaRV=deltaRV.hat,nFull=nFull)
  
  # return
  return(list(g.MCMC=g.sel,g.hat=g.hat,
              nu.MCMC=nu.sel,nu.hat=nu.hat,
              lambda.MCMC=lambda.sel,lambda.hat=lambda.hat,
              deltaRV.MCMC=deltaRV.sel,deltaRV.hat=deltaRV.hat,
              logPost=logPost.sel,logLK=logLK.sel,
              Schwarz=Schwarz,BIC=BIC,check.BIC=check.BIC))
}



#=========================================================================
#' QUALYPSOSS.ANOVA.step2
#'
#' SSANOVA decomposition of the ensemble of climate change responses using a Bayesian approach.
#' In this second step, we infer deltaRV (variance of the residual errors) and phi (autocorrelation lag-1)
#' considering hetero-autocorrelated residual errors, conditionally to smooth effects inferred in \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param yMCMC array \code{nMCMC} x \code{nFull} of climate change responses
#' @param gSum.step1 sum of smooth effect estimates provided by \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param deltaRV.step1 residual variance estimate provided by \code{\link{QUALYPSOSS.ANOVA.step1}}
#'
#' @return list containing diverse information aboutwith the following fields:
#' \itemize{
#'   \item \strong{rho.MCMC}: autocorrelation parameter of the AR(1) process: vector of length \code{n}
#'   \item \strong{deltaRV.MCMC}: Residual variance: vector of length \code{n},
#'   \item \strong{rho.hat}: autocorrelation parameter estimate of the AR(1) process, 
#'   \item \strong{deltaRV.hat}: Residual variance estimate.
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin  
QUALYPSOSS.ANOVA.step2 = function(lOpt, lDim, yMCMC, gSum.step1, deltaRV.step1){
  # retrieve dimensions
  nP = lDim$nP
  nFull = lDim$nFull
  nMCMC = lOpt$nMCMC
  nMO = nFull/nP
  
  # missing responses
  isMiss = is.na(yMCMC[1,])
  nMiss = sum(isMiss)
  
  # climate change response
  y = yMCMC[1,]
  y[isMiss] = gSum.step1[isMiss] + rnorm(n=nMiss, sd = sqrt(deltaRV.step1))
  
  # weigths specifying the heteroscedasticity
  weight.hetero = get.vec.weight.hetero(nP,lOpt$type.hetero)
  
  # initialize MCMC chain for rho
  if(lOpt$type.temporal.dep=="AR1"){
    rho.current = 0.5
  }else{
    # no temporal dependence
    rho.current = 0
  }
  
  # parameters for the Metropolis-Hastings algorithm
  rho.proposal.eps = 0.1 # for the proposal distribution
  cntProp = 0 # count iterations
  optimal.rate = 0.3 # optimal rate of acceptance
  
  # build W
  W = get.matrix.W(weight.hetero,nMO,nP,rho.current)
  Winv = get.matrix.Winv(weight.hetero,nMO,nP,rho.current)
  
  # initialize MCMC samples
  delta.step2.MCMC = rho.step2.MCMC = rho.step2.acc = vector(length=nMCMC)
  
  pb <- txtProgressBar(min = 1, max = nMCMC, style = 3)
  for(iMC in 1:nMCMC){
    cntProp = cntProp + 1
    setTxtProgressBar(pb, iMC)
    y[!isMiss] = yMCMC[iMC,!isMiss]
    
    # noise (hetero + ar1)
    epsRV = (y - gSum.step1)
    
    # sigma2
    distSS = t(epsRV)%*%Winv%*%epsRV
    deltaRV = 1/rgamma(n=1, shape=nFull/2, rate=distSS/2)
    
    #=== rho ===
    
    if(lOpt$type.temporal.dep=="AR1"){
      # every 100 MCMC samples during the burn-in period, we
      # adapt "rho.proposal.eps" which corresponds to the scale
      # of the proposal distribution IF the rate of acceptance
      # ("rate.acc") is too far from 0.3 (considered to be "optimal)
      if(cntProp==100&iMC<lOpt$nBurn){
        rate.acc = mean(rho.step2.acc[(iMC-99):iMC])
        if(abs(rate.acc-optimal.rate)>0.1){
          rho.proposal.eps = rho.proposal.eps*(rate.acc/optimal.rate)
        }
        cntProp = 1
      }
      
      # proposed rho ~ Unif[rho.curr - eps, rho.curr + eps]
      rho.proposed = runif(1)*2*rho.proposal.eps + rho.current - rho.proposal.eps
      
      # if negative values lt -0.99 or gt 0.99, we do not evaluate the proposal and reject this candidate
      if(rho.proposed>(-0.999) & rho.proposed<0.999){
        # log-density of the current value
        rho.current.target.dens = get.target.logdensity.rho(nFull,deltaRV,distSS,weight.hetero,nMO,nP,rho.current)
        
        # log-density of the proposed value: we need to compute W${-1} (Winv) and the sum of squared differences (distSS)
        Winv.proposed = get.matrix.Winv(weight.hetero,nMO,nP,rho.proposed)
        distSS.proposed = (t(epsRV)%*%Winv.proposed%*%epsRV)
        rho.proposed.target.dens = get.target.logdensity.rho(nFull,deltaRV,distSS.proposed,weight.hetero,nMO,nP,rho.proposed)
        
        # MH algorithm: compute transition probabily and accept proposal if r<proba.transition
        # (i.e. with probability proba.transition)
        transition.proba = exp(rho.proposed.target.dens - rho.current.target.dens)
        if(runif(1)<transition.proba){
          rho.current = rho.proposed
          rho.step2.acc[iMC] = 1
          Winv = get.matrix.Winv(weight.hetero,nMO,nP,rho.current)
          W = get.matrix.W(weight.hetero,nMO,nP,rho.current)
        }
      }
    }
    
    # completed obs: normal conditional posterior distribution
    if(any(isMiss)){
      W.noise = W[isMiss,isMiss]*deltaRV
      y[isMiss] = gSum.step1[isMiss] +  mvtnorm::rmvnorm(n = 1,mean=rep(0,nMiss),
                                                         sigma=W.noise)
    }
    
    # store values
    rho.step2.MCMC[iMC] = rho.current
    delta.step2.MCMC[iMC] = deltaRV
  }
  close(pb)
  
  # remove burn-in period
  zMCMC = lOpt$vecKeepMCMC
  rho.sel = rho.step2.MCMC[zMCMC]
  deltaRV.sel = delta.step2.MCMC[zMCMC]
  
  # estimates
  rho.hat = mean(rho.sel)
  deltaRV.hat = mean(deltaRV.sel)
  
  # return
  return(list(rho.MCMC=rho.sel,rho.hat=rho.hat,
              deltaRV.MCMC=deltaRV.sel,deltaRV.hat=deltaRV.hat))
}



#=========================================================================
get.mvtnorm = function(n,detW,Winv,deltaRV,x,mu){
  x.sh = x-mu
  d = -t(x.sh)%*%(Winv/deltaRV)%*%x.sh/2 - n*log(2*pi)/2 - (1/2)*detW - (1/2)*n*log(deltaRV)
  return(d)
}

#=========================================================================
#' QUALYPSOSS.ANOVA.step3
#'
#' SSANOVA decomposition of the ensemble of climate change responses using a Bayesian approach.
#' In this second step, we infer deltaRV (variance of the residual errors) and phi (autocorrelation lag-1)
#' considering hetero-autocorrelated residual errors, conditionally to smooth effects inferred in \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param yMCMC array \code{nMCMC} x \code{nFull} of climate change responses
#' @param RK large object containing the reproducing kernels, returned by \code{\link{QUALYPSOSS.get.RK}}
#' @param g.step1 smooth effect estimates provided by \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param lambda.step1 smooth parameter estimates provided by \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param rho.step2 lag-1 autocorrelation estimate provided by \code{\link{QUALYPSOSS.ANOVA.step2}}
#' @param deltaRV.step2 residual variance estimate provided by \code{\link{QUALYPSOSS.ANOVA.step2}}
#'
#' @return list containing diverse information aboutwith the following fields:
#' \itemize{
#'   \item \strong{g.MCMC}: Smooth effects \code{g}: array \code{n} x \code{nFull} x \code{K} where
#'   \code{nFull} is the number of possible combinations of predictors (discrete AND continuous),
#'   \item \strong{g.hat}: Smooth effects estimates: matrix \code{nFull} x \code{K} where
#'   \code{nFull} is the number of possible combinations of predictors (discrete AND continuous),
#'   \item \strong{Schwarz}: Schwarz criteria
#'   \item \strong{BIC}: BIC criteria
#' }
#' 
#' @export
#' 
#' @author Guillaume Evin  
QUALYPSOSS.ANOVA.step3 = function(lOpt, lDim, yMCMC, RK, g.step1, lambda.step1, rho.step2, deltaRV.step2){
  # retrieve dimensions
  nP = lDim$nP
  nFull = lDim$nFull
  nMCMC = lOpt$nMCMC
  nMO = nFull/nP
  nE = lDim$nE
  
  #====== dimension of the reproducing kernels
  r = sapply(RK, "[[", "r")
  
  # missing responses
  isMiss = is.na(yMCMC[1,])
  nMiss = sum(isMiss)
  
  # initialize quantities
  g.step3.MCMC = array(dim=c(nMCMC,nFull,nE))
  logLK = vector(length=nMCMC)
  
  # weigths specifying the heteroscedasticity
  weight.hetero = get.vec.weight.hetero(nP,lOpt$type.hetero)
  
  # W and W^{-1} with rho estimate
  W = get.matrix.W(weight.hetero,nMO,nP,rho.step2)
  W.noise = W[isMiss,isMiss]
  Winv = get.matrix.Winv(weight.hetero,nMO,nP,rho.step2)
  detW = get.logdet.W(weight.hetero,nMO,nP,rho.step2)
  
  # climate change response
  y = yMCMC[1,]
  if(any(isMiss)){
    y[isMiss] = rowSums(g.step1[isMiss,]) +  
    mvtnorm::rmvnorm(n = 1,mean=rep(0,nMiss),sigma=W.noise)*sqrt(deltaRV.step2)
  }
  
  # sampling strategy: sample direct full conditionals which are normal distributions N(mu3,S3)
  # N(mu3,S3) is the product of two multivariate normal distribution:
  # - likelihood N(mu1,S1) with mu1 = gtilde = y - sum_(e noteq k) g_e and S1 deltaRV W
  # - prior smooth effect N(mu2,S2) with mu2=0 and S2 (deltaRV/lambda)*RK
  #
  # In that case, the product can be obtained as (only one matrix inversion):
  # S3 = S1 (S1+S2)^-1 S2 with S1 = deltaRV, S2 (deltaRV/lambda) RK
  # mu3 = S2 (S1+S2)^-1 mu1 + S1 (S1+S2)^-1 mu2
  # which corresponds to
  # S3 = deltaRV*(W * C^-1 * (1/lambda)*RK) where C = W + (1/lambda)*RK
  # mu3 = (1/lambda)*RK*C^-1*gtilde
   
  gMatmu = g.sigma = list()
  for(iE in 1:nE){
    # C = W + (1/lambda)*RK
    Cinv = MASS::ginv(W + (1/lambda.step1[iE])*RK[[iE]]$SIGMA)
    
    # Covariance of the full conditionals
    # S3 = deltaRV*(W * C^-1 * (1/lambda)*RK) where C = W + (1/lambda)*RK
    S3 = deltaRV.step2*(W%*%Cinv%*%((1/lambda.step1[iE])*RK[[iE]]$SIGMA))
    # make symmetric if needed
    if(!isSymmetric(S3, tol = sqrt(.Machine$double.eps))){
      S3[lower.tri(S3)] = t(S3)[lower.tri(S3)]
    }
    # positive semidefinite
    ev <- eigen(S3, symmetric = TRUE)
    S3.sym = (ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0))))
    
    # store
    g.sigma[[iE]] = S3.sym
    
    # Matrix to compute the mean of the full conditional
    # (1/lambda)*RK*C^-1
    gMatmu[[iE]] = ((1/lambda.step1[iE])*RK[[iE]]$SIGMA)%*%Cinv
  }
  
  
  # resample with the MH algorithm
  #====== i.MCMC = 1..n.MCMC
  g.current = g.step1
  
  pb <- txtProgressBar(min = 1, max = nMCMC, style = 3)
  for(iMC in 1:nMCMC){
    setTxtProgressBar(pb, iMC)
    y[!isMiss] = yMCMC[iMC,!isMiss]
    
    for(iE in 1:nE){
      #===== Metropolis-Hastings algorithm to sample the product of two multivariate distributions
      # N(diffE,deltaRV*W) x N(0,(deltaRV/lambda * D))
      # diffE is the difference between the climate change response y and the sum of the main effects except the current one
      diffE = as.matrix(y-rowSums(g.current[,-iE,drop=FALSE]))
      
      # store values
      retval <- matrix(rnorm(nFull), nrow = 1, byrow = TRUE) %*% g.sigma[[iE]]
      g.current[,iE] = sweep(retval, 2, gMatmu[[iE]]%*%diffE, "+")
    }
    
    # store
    g.step3.MCMC[iMC,,] = g.current
    
    # sum of main effects
    gSum = rowSums(g.current)
    
    # completed obs: normal conditional posterior distribution
    if(any(isMiss)){
      y[isMiss] = gSum[isMiss] +  
      mvtnorm::rmvnorm(n = 1,mean=rep(0,nMiss),sigma = W.noise)*sqrt(deltaRV.step2)
    }
    
    # log-likelihood
    y.std = (y-gSum)/sqrt(deltaRV.step2)
    logLK[iMC] = get.mvtnorm(n=nFull,detW,Winv,1,y.std,rep(0,nFull))-nFull*log(sqrt(deltaRV.step2))
  }
  close(pb)
  
  
  # remove burn-in period
  zMCMC = lOpt$vecKeepMCMC
  g.step3.sel = g.step3.MCMC[zMCMC,,]
  g.step3.hat = apply(g.step3.sel,c(2,3),mean)
  logLK.sel = logLK[zMCMC]
  
  # BIC & Schwarz
  nTheta=sum(r)+nE+2
  Schwarz = max(logLK.sel) - nTheta*log(nFull)/2
  BIC = -2*Schwarz
  
  # check BIC
  Winv = get.matrix.Winv(weight.hetero,nMO,nP,rho.step2)
  SS = t(y.std)%*%Winv%*%y.std
  ww = diag(get.matrix.hetero(weight.hetero,nMO))
  check.BIC = list(y=y,gSum=gSum,err=(y - gSum),
                   std.err = (y-gSum)/(ww*sqrt(deltaRV.step2)),
                   SS=SS,maxLK=max(logLK.sel),nTheta=nTheta,
                   deltaRV=deltaRV.step2,nFull=nFull)
  
  return(list(g.MCMC=g.step3.sel,g.hat=g.step3.hat,Schwarz=Schwarz,BIC=BIC,check.BIC=check.BIC))
}




#==============================================================================
#' formatQUALYPSSoutput
#'
#' @param lOpt list of options, returned by \code{\link{QUALYPSOSS.check.option}}
#' @param lDim list of dimensions
#' @param lScen list of scenario characteristics, output from \code{\link{QUALYPSOSS.process.scenario}}
#' @param ANOVA.step1 list provided by \code{\link{QUALYPSOSS.ANOVA.step1}}
#' @param ANOVA.step2 list provided by \code{\link{QUALYPSOSS.ANOVA.step2}}
#' @param ANOVA.step3 list provided by \code{\link{QUALYPSOSS.ANOVA.step3}}
#' @param climResponse list containing phi, eta, provided by \code{\link{extract.climate.response}}
#' @param change.variable list containing phiStar, etaStar, provided by \code{\link{compute.change.variable}}
#' 
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{POINT}: list containing the mean estimate of different quantities: \code{RESIDUALVAR} 
#'   (residual variability), \code{INTERNALVAR} (internal variability), \code{GRANDMEAN} (grand mean for all time
#'   steps), \code{MAINEFFECT} (list with one item per discrete predictor \code{i}, containing matrices \code{nT} x 
#'   \code{nEffi}, where \code{nEffi} is the number of possible values for the discrete predictor \code{i}).
#'   \code{EFFECTVAR}, uncertainty related to the different main effect, \code{TOTVAR} Total variability,
#'   \code{DECOMPVAR}, decomposition of the total variability (percentages) for the different components,
#'   \code{CONTRIB_EACH_EFFECT}, contribution of each individual effects (percentages) to the corr. effect uncertainty.
#'   \item \strong{BAYES}: list containing quantiles of different estimated quantities, listed in \strong{POINT}.
#'   \item \strong{MCMC}: MCMC draws for the different quantities.
#' }
#'   
#' @export
#' 
#' @author Guillaume Evin
formatQUALYPSSoutput = function(lOpt, lDim, lScen, ANOVA.step1, ANOVA.step2, ANOVA.step3,
                                climResponse, change.variable){
  
  # quantiles from the posterior distributions
  qq = lOpt$quantileCompress
  nqq = length(qq)
  
  # dimensions
  nE = lDim$nE
  nK = lDim$nK
  nP = lDim$nP
  
  # continuous predictors
  XPred = lScen$XFull
  xP = lScen$predContUnique
  
  # discrete predictors
  predDiscreteUnique = lScen$predDiscreteUnique
  nD = unlist(lapply(predDiscreteUnique,length))
  
  
  #_________ initialize lists ______
  # - POINT: point estimates
  # - BAYES: Quantiles from the posterior distributions
  # - MCMC: MCMC simulation chains
  POINT = BAYES = MCMC = list()
  
  
  #_________Residual variability _________
  if(!is.null(ANOVA.step2)){
    # residual var for each future time
    ww = get.vec.weight.hetero(nP,lOpt$type.hetero)^2
    POINT$RESIDUALVAR = ww*ANOVA.step2$deltaRV.hat
    
    # quantiles of the residual var for each future time: matrix nQ x nP
    xq = quantile(ANOVA.step2$deltaRV.MCMC,probs = qq)
    BAYES$RESIDUALVAR = xq%o%ww
    
    # MCMC
    if(lOpt$returnMCMC){
      MCMC$RESIDUALVAR = ANOVA.step2$deltaRV.MCMC
    }
    
  }else{
    # residual var for each future time
    POINT$RESIDUALVAR = rep(ANOVA.step1$deltaRV.hat,nP)
    
    # quantiles of the residual var for each future time: matrix nQ x nP
    xq = quantile(ANOVA.step1$deltaRV.MCMC,probs = qq)
    BAYES$RESIDUALVAR = replicate(nP,xq) # nQ x nP
    
    # MCMC
    if(lOpt$returnMCMC){
      MCMC$RESIDUALVAR = ANOVA.step1$deltaRV.MCMC
    }
  }
  
  
  
  #_________ Internal variability _________
  iv.by.scen = apply(change.variable$etaStar.MCMC[lOpt$vecKeepMCMC,,],c(1,2),var,na.rm=TRUE)
  iv.MCMC = apply(iv.by.scen,1,mean)
  POINT$INTERNALVAR = rep(mean(iv.MCMC),nP)
  BAYES$INTERNALVAR = replicate(nP,quantile(iv.MCMC,probs = qq))
  if(lOpt$returnMCMC){
    MCMC$INTERNALVAR = iv.MCMC
  }
  
  
  #_________ rho _________
  if(lOpt$type.temporal.dep=="AR1"&(!is.null(ANOVA.step2))){
    POINT$RHO = ANOVA.step2$rho.hat
    BAYES$RHO = quantile(ANOVA.step2$rho.MCMC,probs = qq)
    if(lOpt$returnMCMC){
      MCMC$RHO = ANOVA.step2$rho.MCMC
    }
  }
  
  
  #_________ lambda _________
  POINT$LAMBDA = ANOVA.step1$lambda.hat
  BAYES$LAMBDA = apply(ANOVA.step1$lambda.MCMC,2,quantile,probs = qq)
  if(lOpt$returnMCMC){
    MCMC$LAMBDA = ANOVA.step1$lambda.MCMC
  }
  
  
  #_________ grand mean _________
  # initialise
  POINT$GRANDMEAN = vector(length=nP)
  BAYES$GRANDMEAN = matrix(nrow=nqq,ncol=nP)
  if(lOpt$returnMCMC){
    MCMC$GRANDMEAN = matrix(nrow=lOpt$nKeep,ncol=nP)
  }
  
  # retrieve grand mean for each value of the continuous predictor
  if(!is.null(ANOVA.step3)){
    for(iP in 1:nP){
      ii = which(XPred$PredCont==xP[iP])
      POINT$GRANDMEAN[iP] = mean(ANOVA.step3$g.hat[ii,nE])
      BAYES$GRANDMEAN[,iP] = quantile(ANOVA.step3$g.MCMC[,ii,nE],probs=qq)
      if(lOpt$returnMCMC){
        MCMC$GRANDMEAN[,iP] = apply(ANOVA.step3$g.MCMC[,ii,nE],1,mean)
      }
    }
  }else{
    for(iP in 1:nP){
      ii = which(XPred$PredCont==xP[iP])
      POINT$GRANDMEAN[iP] = mean(ANOVA.step1$g.hat[ii,nE])
      BAYES$GRANDMEAN[,iP] = quantile(ANOVA.step1$g.MCMC[,ii,nE],probs=qq)
      if(lOpt$returnMCMC){
        MCMC$GRANDMEAN[,iP] = apply(ANOVA.step1$g.MCMC[,ii,nE],1,mean)
      }
    }
  }
  
  
  #_________ main effects _________
  # initialise lists
  POINT$MAINEFFECT = list()
  BAYES$MAINEFFECT = list()
  if(lOpt$returnMCMC){
    MCMC$MAINEFFECT = list()
  }
  
  # loop over the different main effects (e.g. GCMs, RCMs)
  for(iK in 1:nK){
    # effects for this main effect
    xK = predDiscreteUnique[[iK]]
    
    # initialise objects for this main effect
    nEff = nD[iK]
    POINT$MAINEFFECT[[iK]] = matrix(nrow=nP,ncol=nEff)
    BAYES$MAINEFFECT[[iK]] = array(dim=c(nqq,nP,nEff))
    if(lOpt$returnMCMC){
      MCMC$MAINEFFECT[[iK]] = array(dim=c(lOpt$nKeep,nP,nEff))
    }
    
    # retrieve grand mean for each value of the continuous predictor
    if(!is.null(ANOVA.step3)){
      for(iP in 1:nP){
        for(iEff in 1:nEff){
          ii = which(XPred$PredCont==xP[iP]&XPred[,iK]==xK[iEff])
          POINT$MAINEFFECT[[iK]][iP,iEff] = mean(ANOVA.step3$g.hat[ii,iK])
          BAYES$MAINEFFECT[[iK]][,iP,iEff] = quantile(ANOVA.step3$g.MCMC[,ii,iK],probs=qq)
          if(lOpt$returnMCMC){
            MCMC$MAINEFFECT[[iK]][,iP,iEff] = apply(ANOVA.step3$g.MCMC[,ii,iK],1,mean)
          }
        }
      }
    }else{
      for(iP in 1:nP){
        for(iEff in 1:nEff){
          ii = which(XPred$PredCont==xP[iP]&XPred[,iK]==xK[iEff])
          POINT$MAINEFFECT[[iK]][iP,iEff] = mean(ANOVA.step1$g.hat[ii,iK])
          BAYES$MAINEFFECT[[iK]][,iP,iEff] = quantile(ANOVA.step1$g.MCMC[,ii,iK],probs=qq)
          if(lOpt$returnMCMC){
            MCMC$MAINEFFECT[[iK]][,iP,iEff] = apply(ANOVA.step1$g.MCMC[,ii,iK],1,mean)
          }
        }
      }
    }
  }
  
  
  #_________ uncertainty due to the main effect _________
  # initialize matrix
  POINT$EFFECTVAR = matrix(nrow=nP,ncol=nK)
  for(iK in 1:nK){
    for(iP in 1:nP){
      POINT$EFFECTVAR[iP,iK] = mean(POINT$MAINEFFECT[[iK]][iP,]^2)
    }
  }
  
  
  #_____ Total variability and its decomposition ________
  # collect variabilities and uncertainties in a matrix
  MATVAR = cbind(POINT$EFFECTVAR,POINT$RESIDUALVAR,POINT$INTERNALVAR)
  colnames(MATVAR) = c(colnames(lScen$scenAvail),"ResidualVar","InternalVar")
  
  # total variability
  POINT$TOTVAR = rowSums(MATVAR)
  
  # variance decomposition
  POINT$DECOMPVAR = MATVAR/replicate(ncol(MATVAR),POINT$TOTVAR)
  
  
  #_______ Contribution of each effect _____________
  POINT$CONTRIB_EACH_EFFECT = list()
  
  # loop over the different main effects (e.g. GCMs, RCMs)
  for(iK in 1:nK){
    nEff = nD[iK]
    POINT$CONTRIB_EACH_EFFECT[[iK]] = matrix(nrow=nP,ncol=nEff)
    
    for(iP in 1:nP){
      for(iEff in 1:nEff){
        POINT$CONTRIB_EACH_EFFECT[[iK]][iP,iEff] = mean(POINT$MAINEFFECT[[iK]][iP,iEff]^2)/(POINT$EFFECTVAR[iP,iK]*nEff)
      }
    }
  }
  
  
  #______________ Climate response ________________
  # initialise climate response outputs
  CLIMATE_RESPONSE = list()
  CLIMATE_RESPONSE$POINT = list()
  CLIMATE_RESPONSE$BAYES = list()
  if(lOpt$returnMCMC) CLIMATE_RESPONSE$MCMC = list()
  
  # phi
  phi = climResponse$phi.MCMC[lOpt$vecKeepMCMC,,]
  CLIMATE_RESPONSE$POINT$phi = apply(phi,c(2,3),mean)
  CLIMATE_RESPONSE$BAYES$phi = apply(phi,c(2,3),quantile,probs = qq)
  if(lOpt$returnMCMC) CLIMATE_RESPONSE$MCMC$phi = phi
  
  # eta
  eta = climResponse$eta.MCMC[lOpt$vecKeepMCMC,,]
  CLIMATE_RESPONSE$POINT$eta = apply(eta,c(2,3),mean)
  CLIMATE_RESPONSE$BAYES$eta = apply(eta,c(2,3),quantile,probs = qq)
  if(lOpt$returnMCMC) CLIMATE_RESPONSE$MCMC$eta = eta
  
  # phiStar
  phiStar = change.variable$phiStar.MCMC[lOpt$vecKeepMCMC,,]
  CLIMATE_RESPONSE$POINT$phiStar = apply(phiStar,c(2,3),mean)
  CLIMATE_RESPONSE$BAYES$phiStar = apply(phiStar,c(2,3),quantile,probs = qq)
  if(lOpt$returnMCMC) CLIMATE_RESPONSE$MCMC$phiStar = phiStar
  
  # etaStar
  etaStar = change.variable$etaStar.MCMC[lOpt$vecKeepMCMC,,]
  CLIMATE_RESPONSE$POINT$etaStar = apply(etaStar,c(2,3),mean)
  CLIMATE_RESPONSE$BAYES$etaStar = apply(etaStar,c(2,3),quantile,probs = qq)
  if(lOpt$returnMCMC) CLIMATE_RESPONSE$MCMC$etaStar = etaStar
  
  
  return(list(POINT=POINT,BAYES=BAYES,MCMC=MCMC,CLIMATE_RESPONSE=CLIMATE_RESPONSE))
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
#'   \item \strong{spar}: if \code{uniqueFit} is true, smoothing parameter passed to the function \link[stats]{smooth.spline}.
#'   \item \strong{lambdaClimateResponse}: smoothing parameter > 0 for the extraction of the climate response.
#'   \item \strong{lambdaHyperParANOVA}: hyperparameter \eqn{b} for the \eqn{\lambda} parameter related to each predictor \eqn{g}.
#'   \item \strong{typeChangeVariable}: type of change variable: "abs" (absolute, value by default) or "rel" (relative).
#'   \item \strong{nBurn}: number of burn-in samples (default: 1000). If \code{nBurn} is too small, the convergence of MCMC chains might not be obtained.
#'   \item \strong{nKeep}: number of kept samples (default: 2000). If \code{nKeep} is too small, MCMC samples might not be represent correctly the posterior
#'   distributions of inferred parameters.
#'   \item \strong{quantileCompress}: vector of probabilities (in [0,1]) for which we compute the quantiles from the posterior distributions
#'    \code{quantileCompress = c(0.005,0.025,0.05,0.5,0.95,0.975,0.995)} by default.
#'   \item \strong{uniqueFit}: logical, if \code{FALSE} (default), climate responses are fitted using Bayesian smoothing splines, otherwise,if \code{TRUE},
#'    a unique cubic smoothing spline is fitted for each run, using the function \link[stats]{smooth.spline}.
#'   \item \strong{returnMCMC}: logical, if \code{TRUE}, the list \code{MCMC} contains MCMC chains.
#'   \item \strong{returnOnlyCR}: logical, if \code{TRUE} (default), only Climate Responses are fitted and returned.
#'   \item \strong{type.temporal.dep}: "iid" for independent errors or "AR1" (default) for autocorrelated errors.
#'   \item \strong{type.hetero}: "constant" for homoscedastic errors or "linear" (default) for heteroscedastic errors.
#'   
#' }
#' @param RK Reproducing kernels: list
#' 
#' @return  list with the following fields:
#'
#' \itemize{
#'   \item \strong{POINT}: list containing the mean estimate of different quantities: \code{RESIDUALVAR} 
#'   (residual variability), \code{INTERNALVAR} (internal variability), \code{GRANDMEAN} (grand mean for all time
#'   steps), \code{MAINEFFECT} (list with one item per discrete predictor \code{i}, containing matrices \code{nT} x 
#'   \code{nEffi}, where \code{nEffi} is the number of possible values for the discrete predictor \code{i}).
#'   \code{EFFECTVAR}, uncertainty related to the different main effect, \code{TOTVAR} Total variability,
#'   \code{DECOMPVAR}, decomposition of the total variability (percentages) for the different components,
#'   \code{CONTRIB_EACH_EFFECT}, contribution of each individual effects (percentages) to the corr. effect uncertainty.
#'   \item \strong{BAYES}: list containing quantiles of different estimated quantities, listed in \strong{POINT}.
#'   \item \strong{MCMC}: list containing the MCMC chains (not returned by default).
#'   \item \strong{climateResponse}: list containing different objects related to the extraction of the climate response.
#'   phiStar (\eqn{\phi^*}) is an array \code{nQ} x \code{nS} x \code{nP} containing climate change responses, where \code{nQ} is the
#'   number of returned quantiles, \code{nS} is the number of scenarios and \code{nP} is the length of \code{predContUnique} (e.g. number
#'   of future years). 
#'   Similarly, etaStar (\eqn{\eta^*}) contains the deviation from the climate change response. 
#'   phi (\eqn{\phi}) contains the climate responses and eta (\eqn{\eta}) contains the deviations from the climate responses.
#'   \item \strong{listCR}: list containing objects created during the extraction of the climate responses
#'   \item \strong{ClimateProjections}: argument of the call to the function, for records.
#'   \item \strong{predCont}: (optional) argument of the call to the function, for records.
#'   \item \strong{predContUnique}:  (optional) argument of the call to the function, for records.
#'   \item \strong{predDiscreteUnique}: list of possible values taken by the discrete predictors given in \code{scenAvail}.
#'   \item \strong{listOption}: list of options
#'   \item \strong{listScenario}: list of scenario characteristics (obtained from \code{\link{QUALYPSOSS.process.scenario}})
#'   \item \strong{RK}: list containing the reproducing kernels
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
#' listOption = list(nBurn=20,nKeep=30,type.temporal.dep="iid",type.hetero="constant")
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
  
  # list of dimensions
  lDim = list(nT=nT, nP=nP, nK=nK, nE=nK+1, nS=nS, nFull=lScen$nFull)
  
  
  #=============================================================
  #                        CLIMATE RESPONSE
  #=============================================================
  
  #=================== EXTRACT CLIMATE RESPONSE ===================
  climResponse = extract.climate.response(ClimateProjections=ClimateProjections,
                                          predCont=predCont,
                                          predContUnique=predContUnique,
                                          nMCMC = lOpt$nMCMC,
                                          lam=lOpt$lambdaClimateResponse,
                                          uniqueFit = lOpt$uniqueFit,
                                          spar = lOpt$spar)
  
  
  
  #=================== COMPUTE CLIMATE CHANGE RESPONSE ===================
  # Compute change variables
  change.variable = compute.change.variable(climResponse,lOpt,lDim,iCpredContUnique,iCpredCont)
  
  # reformat climate responses
  yMCMC = get.yMCMC(lOpt,lDim,lScen,change.variable)
  
  
  #=============================================================
  #                           ANOVA
  #=============================================================
  
  # reproducing kernels (computationally intensive)
  returnRK = is.null(RK)
  if(is.null(RK)) RK = QUALYPSOSS.get.RK(X=lScen$XFull,nK=nK)
  
  # ANOVA: STEP 1: no hetero or autocorr
  ANOVA.step1 = QUALYPSOSS.ANOVA.step1(lOpt, lDim, yMCMC, RK)
  
  # if we consider autocorrelated residual errors, we continue the inference
  # -> steps 2 & 3
  if(lOpt$type.temporal.dep=="AR1"|lOpt$type.hetero!="constant"){
    # step 2: infer delta_RV and rho_RV conditionally to g.hat obtained in step1
    # (deltaRV just provided for initialisation)
    ANOVA.step2 = QUALYPSOSS.ANOVA.step2(lOpt = lOpt,lDim =  lDim,yMCMC = yMCMC,
                                         gSum.step1 = rowSums(ANOVA.step1$g.hat),
                                         deltaRV.step1 =  ANOVA.step1$deltaRV.hat)
    
    # step 3: update smooth effects g conditionally to delta_RV, rho_RV (step 2) and 
    # lambda (step 1). g.step1 & nu.step1 are just provided for initialisation
    ANOVA.step3 = QUALYPSOSS.ANOVA.step3(lOpt = lOpt,lDim =  lDim,yMCMC = yMCMC, RK=RK,
                                         g.step1 = ANOVA.step1$g.hat,
                                         lambda.step1 = ANOVA.step1$lambda.hat,
                                         rho.step2 = ANOVA.step2$rho.hat,
                                         deltaRV.step2 = ANOVA.step2$deltaRV.hat)
  }else{
    ANOVA.step2 = NULL
    ANOVA.step3 = NULL
  }
  
  
  
  #=============================================================
  #                       RETURN RESULTS
  #=============================================================
  
  # if returnRK is FALSE, we do not return RK (large object)
  if(!returnRK){
    RK = NULL
  }
  
  # format outputs
  output.MCMC = formatQUALYPSSoutput(lOpt, lDim, lScen, ANOVA.step1, ANOVA.step2, ANOVA.step3,
                                     climResponse,change.variable)
  
  # return list of results
  return(list(POINT=output.MCMC$POINT,BAYES=output.MCMC$BAYES,MCMC=output.MCMC$MCMC,
              CLIMATE_RESPONSE=output.MCMC$CLIMATE_RESPONSE,listCR = climResponse$listCR,
              ClimateProjections=ClimateProjections,predCont=predCont,
              predContUnique=predContUnique,predDiscreteUnique=lScen$predDiscreteUnique,
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
  phi = QUALYPSOSSOUT$CLIMATE_RESPONSE$POINT$phi
  
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
  phiStar = QUALYPSOSSOUT$CLIMATE_RESPONSE$POINT$phiStar
  
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
  med = QUALYPSOSSOUT$BAYES$GRANDMEAN[iMedian,][Iord]
  binf = QUALYPSOSSOUT$BAYES$GRANDMEAN[iBinf,][Iord]
  bsup = QUALYPSOSSOUT$BAYES$GRANDMEAN[iBsup,][Iord]
  
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
plotQUALYPSOSSeffect = function(QUALYPSOSSOUT,iEff,CIlevel=c(0.025,0.975),lim=NULL,
                                col=1:20,xlab="Continuous predictor",ylab="Effect",addLegend=TRUE,
                                ...){
  # Continuous predictor
  TTmix = QUALYPSOSSOUT$predContUnique
  Iord = order(TTmix)
  nTT = length(TTmix)
  TT = TTmix[Iord]
  
  # retrieve effects
  QUANTEff = QUALYPSOSSOUT$BAYES$MAINEFFECT[[iEff]]
  nEff = dim(QUANTEff)[3]
  
  # find index quantiles
  iMedian = which(QUALYPSOSSOUT$listOption$quantileCompress==0.5)
  iBinf = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[1])
  iBsup = which(QUALYPSOSSOUT$listOption$quantileCompress==CIlevel[2])
  
  # retrieve median and limits
  med = QUANTEff[iMedian,,]
  binf = QUANTEff[iBinf,,]
  bsup = QUANTEff[iBsup,,]
  
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
  VV = QUALYPSOSSOUT$POINT$DECOMPVAR
  nVV = ncol(VV)
  
  # figure
  col = col[1:nVV]
  cum=rep(0,nTT)
  plot(-1,-1,xlim=range(TT),ylim=c(0,1),xaxs="i",yaxs="i",las=1,
       xlab=xlab,ylab=ylab,...)
  for(i in 1:nVV){
    cumPrevious = cum
    cum = cum + VV[,i]
    polygon(c(TT,rev(TT)),c(cumPrevious,rev(cum)),col=rev(col)[i],lty=1)
  }
  abline(h=axTicks(side=2),col="black",lwd=0.3,lty=1)
  
  # legend
  if(addLegend){
    if(is.null(colnames(QUALYPSOSSOUT$listScenario$scenAvail))){
      namesEff = paste0("Eff",1:(nVV-2))
    }else{
      namesEff = colnames(QUALYPSOSSOUT$listScenario$scenAvail)
    }
    legend('topleft',bty='n',cex=1.1, fill=rev(col), legend=c(namesEff,'Res. Var.','Int. Variab.'))
  }
}