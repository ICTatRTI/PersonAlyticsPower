# TODO: save data for use in SAS

# TODO: to move to the description file imports
library(MASS)
library(tidyr)

#' simICTdata - this is a user interface function for one dataset,
#' switching between different designs so the function can be updated
#' @author Stephen Tueller \email{stueller@@rti.org}

simICTdata <- function(design      = polyICT$new()           ,
                       randFxParms = list(mu=.5 , sigma=1)   ,
                       family      = "qNO"                   ,
                       randFxSeed  = 1                       ,
                       errorParms  = list(ar=c(.5), ma=c(0)) ,
                       errorFUN    = arima.sim               ,
                       errorSeed   = 2                       )
{
  # simulate random effects
  randFx <- mvrFam(design, randFxParms, family, randFxSeed)

  # simulate the errors
  errors <- ICTerror(design, errorParms, errorFUN, randFxSeed)


  # polynomial ICT
  if('polyICT' %in% class(design))
  {
    data <- simPolyICT(design, randFxSeed, errorSeed)
  }

  # TODO see ALDA pp 234-235 for implementing these, we may need to create a
  # separate function for each ala `simPolyICT`, or we may be able to do
  # some code injection via formula

  # hyperbolic
  if('hyperICT' %in% class(design))
  {

  }

  # inverse polynomial
  if('invPolyICT' %in% class(design))
  {

  }

  # exponential
  if('expICT' %in% class(design))
  {

  }

  # negative exponential
  if('negExpICT' %in% class(design))
  {

  }

  # logistic
  if('logisticICT' %in% class(design))
  {

  }

  return(data)

}


#' simPolyICT
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param parms A named list of parameters to be passed to a quantile function from
#' \code{\link{gamlss.family}}, specifically mu, sigma, tau, and nu. If \code{family}
#' doesn't use one of these parameters it will be ignored.
#' @param family Default "qNO". A quoted \code{\link{gamlss.family}} quantile
#' distribution.
#' @param seed Default 1. A random seed for replication.

simPolyICT <- function(design, randFxSeed, errorSeed, family)
{


}



#' mvrFam - a function to simulate multivariate data from a gamlss family
#' distribution
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' follow http://www.econometricsbysimulation.com/2014/02/easily-generate-correlated-variables.html
#'
#' @param design An \code{\link{polyICT}} object, or its not yet implemented siblings
#' \code{hyperICT}, \code{invPolyICT}, \code{expICT}, \code{negExpICT},
#' \code{logisticICT}.
#' @param parms A named list of parameters to be passed to a quantile function from
#' \code{\link{gamlss.family}}, specifically mu, sigma, tau, and nu. If \code{family}
#' doesn't use one of these parameters it will be ignored.
#' @param family Default "qNO". A quoted \code{\link{gamlss.family}} quantile
#' distribution.
#' @param seed Default 1. A random seed for replication.
#'
#' @examples
#' parms  <- list(mu=.5, sigma=.35, nu=.15, tau=.25)
#' design = polyICT$new(n=1000)
#' par(mfrow=c(2,2))
#'
#' # beta distribution
#' YBEINF <- mvrFam(design = design, parms = parms, family="qBEINF")
#'
#' # use the same parameters for normal, ignoring nu and tau
#' YNO    <- mvrFam(design = design, parms = parms, family="qNO")
#'
#' # visualize
#' psych::pairs.panels(data.frame(YNO, YBEINF))
#'
#' # quantile functions not in gamlss.family also work
#' norm <- mvrFam(design, parms=list(mean=0, sd=20, lower.tail=FALSE), "qnorm")
#' #pois <- mvrFam(design, parms=list(lambda=.3, lower.tail=FALSE), "qpois") # not working
#'
#' # this example illustrates how mvrFam works under the hodo
#' \dontrun{
#' corMat <- matrix(.5, nrow=3, ncol=3) + diag(3)*.5
#' corMat[1,3] <- corMat[3,1] <- .25
#' effectSizes <- list(randFx=list(intercept=0, slope=.5, quad=.1),
#'                     fixdFx=list(phase=.5, phaseTime=.5, phaseTime2=.04))
#' randFxVar <- c(1,.2,.1)
#'
#'
#' design <- polyICT$new(n=10000, corMat=corMat, effectSizes=effectSizes,
#'                       randFxVar=randFxVar)
#'
#' # first simulate multivariate normal data
#' Y <- mvrnorm( design$n, rep(0, design$polyOrder + 1),  design$corMat)
#' psych::pairs.panels(Y) # TODO make the plots not run
#' all.equal(cor(Y), corMat, tolerance=.03)
#'
#' # now get the propabilities
#' YpNorm <- data.frame(pnorm(Y))
#' psych::pairs.panels(YpNorm)
#' all.equal(unname(cor(YpNorm)), design$corMat, tolerance=.08)
#'
#' # use the quantile functions on the probabilities to get non-normal data
#' # with approximately the same correlation matrix
#'
#' # works well for symmetric distributions
#' YBEINF <- doLapply(YpNorm, .fcn="qBEINF", mu=.5, sigma=.35)
#' psych::pairs.panels(YBEINF)
#' all.equal(unname(cor(YBEINF)), corMat, tolerance=.08)
#'
#' # does not work as well for skewed distributions
#' YLOGNO <- doLapply(YpNorm, "qLOGNO", mu=3, sigma=1)
#' psych::pairs.panels(YLOGNO)
#' all.equal(unname(cor(YLOGNO)), corMat, tolerance=.08)
#'
#'
#' }
mvrFam <- function(design, parms, family="qNO", seed=1)
{
  # attach the parameters so they can be passed by name
  attach(parms)

  # first simulate multivariate normal data
  set.seed(seed)
  Y <- mvrnorm( design$n, rep(0, design$polyOrder + 1), design$corMat)

  if(! family %in% c("qNo", "qnorm"))
  {
    # now get the propabilities
    YpNorm <- data.frame(pnorm(Y))

    # transform to the specified distribution
    Y <- doLapply(YpNorm, family, mu=mu, sigma=sigma, nu=nu, tau=tau, ...)
  }

  # rescale the data - this may be problamtic if the user specifies a
  # discrete distribution for random effects, may add a warning here
  Y <- scale(Y, FALSE, TRUE) * sqrt(1-design$propErrVar)

  # return the data
  return(Y)
}



#' ICTerror
#' @author Stephen Tueller \email{stueller@@rti.org}
#'

ICTerror <- function(design     = polyICT$new()           ,
                     errorParms = list(ar=c(.5),ma=c(0))  ,
                     errorFUN   = arima.sim               ,
                     seed       = 1                       ,
                     ...                                  )
{
  # get the number of Observations
  nObservations <- length(c(unlist((design$phases))))

  # get seeds
  set.seed(seed)
  seeds <- as.list( ceiling(runif(design$n, 0, 9e6) ) )

  # sim errors, use sd = 1 for now, it can be rescaled in other functions
  FUN <- function(x, errorFUN, model, n=design$n, sd=1, ...)
  {
    set.seed(x)
    errorFUN(model = parms, n = n, sd = sd, ...)
  }
  errors <- lapply(seeds, FUN,
                   errorFUN = errorFUN      ,
                   model    = parms         ,
                   n        = nObservations ,
                   sd       = 1             ,
                   ...)
  errors <- data.frame(do.call(rbind, errors))

  # rescale to have unit variances within time points then multiply
  # by the proportion of error variance; need to transpose twice so
  # that the rescaling applies within person (which also gets the
  # between person scaling right, but not vice versa in initial tests)
  errorsr <- t(scale(t(errors), FALSE, TRUE)) * sqrt(design$propErrVar)

  doQC <- FALSE
  if(doQC)
  {
    # QC, only holds asymptotically
    # within persons:
    all.equal(mean(apply(errorsr, 1, var)), design$propErrVar, tolerance = .05)
    # between persons:
    all.equal(mean(apply(errorsr, 2, var)), design$propErrVar, tolerance = .05)
    tseries <- lapply(data.frame(t(errorsr)), ts)
    aFun <- function(x, order=c(1,0,0)) try(arima(x, order), silent = TRUE)
    sigma2s <- lapply(lapply(tseries, aFun), function(x) x$sigma2)
    # tends toward underestimation when nObservations are small, but very precise
    # when large
    hist(unlist(sigma2s))
  }

  # return
  return(errorsr)
}
















