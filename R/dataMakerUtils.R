
#' makeEta - a wrapper for mvrnorm that currenty does nothing else but set the seed
#' and call mvrnorm. Retained as we may need to extend its functionality
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
#' Y <- mvrnorm( design$n, rep(0, design$maxRandFx + 1),  design$randFxCorMat)
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
#' YBEINF <- doLapply(YpNorm, .fcn="qBEINF", mu=.5, sigma=.35)
#' psych::pairs.panels(YBEINF)
#' all.equal(unname(cor(YBEINF)), randFxCorMat, tolerance=.08)
#'
#' # does not work as well for skewed distributions
#' YLOGNO <- doLapply(YpNorm, "qLOGNO", mu=3, sigma=1)
#' psych::pairs.panels(YLOGNO)
#' all.equal(unname(cor(YLOGNO)), randFxCorMat, tolerance=.08)
#'
#'
#' }
makeFam <- function(Y, parms, family="NO")
{
  # qc inputs
  checkFam(family, parms)
  family <- paste('q', family, sep='')

  # get the propabilities
  YpNorm <- data.frame(pnorm(Y))

  # transform to the specified distribution
  Y <- doLapply(YpNorm, family, mu=parms$mu, sigma=parms$sigma,
                nu=parms$nu, tau=parms$tau, ...)

  # rescale the data - this may be problamtic if the user specifies a
  # discrete distribution for random effects, may add a warning here
  # Note: zero mean is specified hear b/c the mean off the random effects
  # is added later in the $makeData method of an designICT child class. This
  # may not be the best approach, especially with non-normal data, but I don't
  # know whether the method of running pnorm through q* quantile functions
  # preserves means at all
  # -- Note: rescaling not done here, instead we sim according to design$randFxCovMat
  # and choose the error variance commensurate to design$propErrVar

  #Y <- scale(Y, TRUE, TRUE) * sqrt(1-design$propErrVar)

  return(Y)
}
