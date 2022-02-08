# this file will contain the the *Data functions [consider moving to repsective *ICT files]
#
#  .polyData
#  hyperData
#  invPolyData
#  expData
#  negExpData
#  logisticData
#
#  that will be used by
#
#  polyICT$makeData()
#  hyperICT$makeData()
#  invPolyICT$makeData()
#  expICT$makeData()
#  negExpICT$makeData()
#  logisticICT$makeData()
#
# These functions do not simulate data, rather they construct
# latent growth curve model data for one group within one phase
# from inputs of random effects and errrors.
#
# Mutligroup and/or multiphase data are achieved by concatenating data
# from multiple calls to *Data functions within *makeData methods.

# level 2 covariates -
# -- categorical: acheived via multiple groups
# -- continuous: implemented directly in lgmSim

# time varying covariates -
# -- continuous: implemented directly in lgmSim
# -- categorical: achieve by categarizing continuous covariates


#' .polyData - function to simulate polynomial growth data, used by the polyICT
#' class
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
#'
.polyData <- function(seed=123, n, nObs, mu, Sigma, self, dMsim, dM,
                     rFxVr, propErrVar, group=NULL)
{
  # get seeds
  # 1. random effects
  # 2. autocorrelated errors
  # 3. white noise errors
  seeds <- .makeSeeds(seed, 3)

  # random effects
  reVars <- rFxVr/sum(rFxVr)
  eta <- makeEta(n, mu, Sigma, seeds[1])
  seta <- scale(eta)
  attr(seta, "scaled:center") <- NULL
  attr(seta, "scaled:scale") <- NULL


  # rescale random effects by  by propErrVar and reVars, retaining means
  eta <- seta %*% diag( sqrt(propErrVar[1] * reVars) ) +
    matrix(mu, nrow(eta), ncol(eta), byrow=TRUE) # `mu` was `apply(eta, 2, mean)`

  # qc - should hold even in small samples
  #all(propErrVar[1] * reVars, apply(eta, 2, var))


  # expand random effects into matrices
  .L  <- list()
  for(i in 1:ncol(eta))
  {

    if(i==1) etaM <- matrix(eta[,i]) %*% rep(1, nObs)
    if(i>1)
    {
      # multiply by the time points
      etaM <- matrix(eta[,i]) %*% dMsim[,i]
      etaV <- scale(etaM)
      etaV[is.nan(etaV[,1]),1] <- 0
      # ensure the variances are exact
      etaM <- etaV*matrix(sd(eta[,i]), nrow(etaM), ncol(etaM), byrow=T) +
                   matrix(apply(etaM, 2, mean), nrow(etaM), ncol(etaM), byrow=T)
    }
    .L[[i]] <- etaM
  }

  # correct to target variances
  Y <- Reduce('+', .L)
  trueVar <- apply(do.call(rbind, lapply(.L, function(x) apply(x, 2, var))), 2, sum)
  Y <- scale(Y)*matrix(sqrt(trueVar),   nrow(Y), ncol(Y), byrow=T) +
    matrix(apply(Y,2,mean), nrow(Y), ncol(Y), byrow = T)

  # get total error variances by subtraction
  errVar <- 1-apply(Y, 2, var)

  # rescale err and merr variances to sum to 1
  errorVar <- propErrVar[2:3]
  errorVar <- errorVar/max(errorVar)

  # errors w/ scaling
   err <- scale(self$error$makeErrors(n, nObs, seeds[2]))
  merr <- scale(self$merror$makeErrors(n, nObs, seeds[3]))
  attr(err, "scaled:center")  <- NULL
  attr(err, "scaled:scale")   <- NULL
  attr(merr, "scaled:center") <- NULL
  attr(merr, "scaled:scale")  <- NULL
   err <- err*matrix(errorVar[1]*sqrt(errVar), nrow(err), ncol(err), byrow=T)
  merr <- merr*matrix(errorVar[2]*sqrt(errVar), nrow(merr), ncol(merr), byrow=T)

  # construct observed data
  Y <- Y + err + merr

  # create id and prepend group number to ensure unique ids
  id <- sort(rep(1:n, nObs))
  if(!is.null(group))
  {
    groupNumber <- as.numeric( gsub("[^[:digit:]]", "", group) )
    groupNumber <- (10*groupNumber)^max(nchar(id))
    id <- groupNumber + id
  }

  # construct a long dataset
  data <- data.frame(
    id    = id                                       ,
    y     = as.vector(t(Y))                          ,
    do.call(rbind, replicate(n, dM, simplify = FALSE))
  )

  # QC
  all(aggregate(data$y, list(data$Time), mean)$x ==
            unname(apply(Y, 2, mean)))
  all(aggregate(data$y, list(data$Time), var)$x ==
            unname(apply(Y, 2, var)))

  # make phase a factor (for later analysis and plotting)
  data$phase <- factor(data$phase)

  # add group identifiers, phase is already coming in via dM
  if(!is.null(group)) data$group <- group

  # return the data
  return(data)
}



#' makeEta - a wrapper for mvrnorm that currently does nothing else but set the
#' seed and call mvrnorm. Retained as we may need to extend its functionality
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
#' @param n See \code{\link{mvrnorm}}.
#' @param mu See \code{\link{mvrnorm}}.
#' @param Sigma See \code{\link{mvrnorm}}.
#' @param seed Default 1. A random seed for replication.
#'
#' @examples
#' Y  <- makeEta(1000, c(0,0), matrix(c(1,.3,.3,1), 2, 2))
#' pairs(Y)
makeEta <- function(n, mu, Sigma, seed=1)
{
  # simulate multivariate normal data
  set.seed(seed)
  Y <- mvrnorm(n, mu, Sigma)

  # return the data
  return(Y)
}


#' makeFam - a function to transform multivariate normal data to that with a
#' gamlss family distribution, following
#' "http://www.econometricsbysimulation.com/2014/02/easily-generate-correlated-variables.html"
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
#' @param Y Multivariate normal data such as that simulated by \code{\link{rnorm}},
#' \code{\link{mvrnorm}}, or \code{\link{arima.sim}}.
#' @param parms A named list of parameters for transforming normal random effects
#' to some other distribution. These will be passed to a quantile function from
#' \code{\link{gamlss.family}}, specifically mu, sigma, tau, and nu. If \code{family}
#' doesn't use one of these parameters it will be ignored.
#' @param family Default "qNO". A quoted \code{\link{gamlss.family}} quantile
#' distribution.
#'
#' @examples
#' parms  <- list(mu=.5, sigma=.35, nu=.15, tau=.25)
#' design = polyICT$new(n=1000)
#' par(mfrow=c(2,2))
#'
#' # beta distribution
#' YBEINF <- makeFam(design = design, parms = parms, family="qBEINF")
#'
#' # use the same parameters for normal, ignoring nu and tau
#' YNO    <- makeFam(design = design, parms = parms, family="qNO")
#'
#' # visualize
#' psych::pairs.panels(data.frame(YNO, YBEINF))
#'
#' # quantile functions not in gamlss.family also work
#' norm <- makeFam(design, parms=list(mean=0, sd=20, lower.tail=FALSE), "qnorm")
#' #pois <- makeFam(design, parms=list(lambda=.3, lower.tail=FALSE), "qpois") # not working
#'
#' # this example illustrates how makeFam works under the hood
#' \dontrun{
#' # requires library(psych), not required by PersonAlyticsPower
#' randFxCorMat <- matrix(.5, nrow=3, ncol=3) + diag(3)*.5
#' randFxCorMat[1,3] <- randFxCorMat[3,1] <- .25
#' randFxMean <- list(randFx=list(intercept=0, slope=.5, quad=.1),
#'                     fixdFx=list(phase=.5, phaseTime=.5, phaseTime2=.04))
#' randFxVar <- c(1,.2,.1)
#'
#'
#' design <- polyICT$new(n=10000, randFxCorMat=randFxCorMat, randFxMean=randFxMean,
#'                       randFxVar=randFxVar)
#'
#' # first simulate multivariate normal data
#' Y <- mvrnorm( design$n, rep(0, design$randFxOrder + 1),  design$randFxCorMat)
#' psych::pairs.panels(Y)
#' all.equal(cor(Y), randFxCorMat, tolerance=.03)
#'
#' # now get the propabilities
#' YpNorm <- data.frame(pnorm(Y))
#' psych::pairs.panels(YpNorm)
#' all.equal(unname(cor(YpNorm)), design$randFxCorMat, tolerance=.08)
#'
#' # use the quantile functions on the probabilities to get non-normal data
#' # with approximately the same correlation matrix
#'
#' # works well for symmetric distributions
#' YBEINF <-. doLapply(YpNorm, .fcn="qBEINF", mu=.5, sigma=.35)
#' psych::pairs.panels(YBEINF)
#' all.equal(unname(cor(YBEINF)), randFxCorMat, tolerance=.08)
#'
#' # does not work as well for skewed distributions
#' YLOGNO <- .doLapply(YpNorm, "qLOGNO", mu=3, sigma=1)
#' psych::pairs.panels(YLOGNO)
#' all.equal(unname(cor(YLOGNO)), randFxCorMat, tolerance=.08)
#'
#'
#' }
makeFam <- function(Y, parms, family="NO")
{
  # qc inputs
  .checkFam(family, parms)
  family <- paste('q', family, sep='')

  # get the propabilities
  YpNorm <- data.frame(pnorm(Y))

  # transform to the specified distribution
  Y <- .doLapply(YpNorm, family, mu=parms$mu, sigma=parms$sigma,
                 nu=parms$nu, tau=parms$tau, ...)

  # rescale the data - this may be problamtic if the user specifies a
  # discrete distribution for random effects, may add a warning here

  return(Y)
}



