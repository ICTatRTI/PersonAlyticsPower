if(1==2)
{

# RESUME HERE: this needs to be moved to the make data
# loops b/c it is phase/group specific; or, move it
# to checkPolyICT; also propErrVar should be allowed to be
# phase and group specific
# get the total variance and the error variance @ time = 1
#*# moving from polyICT to either polyICT$makeData or a dataMaker function
totalVar  <- sum(randFxCovMat)/(1-propErrVar)
errorVar  <- propErrVar * totalVar
variances <- list(totalVar = totalVar,
				 errorVar = errorVar,
				 randFxVar = totalVar - errorVar)

# get the expected varariances
# initial values
expectedVariances <- expectedVar(randFxCovMat, designMat, variances,
								randFxMean, nObservations, n)
unStdInputMat      <- expectedVariances$randFxMean

#*# find homes for these documentation items
#'
#' @param errorParms Parameters for simulating random errors. Currently only
#' supports \code{\link{arima.sim}}.
#'
#' @param errorFUN A function for simulating errors. Currently only
#' \code{\link{arima.sim}} is supported.
#'
#' @param errorFamily A quantile function fo simulating non-normal errors. Not
#' currently implemented.
#'
#' @param randFxFamily A quantile function such as those in
#' \code{\link{gamlss.family}}.
#' If \code{family} is not "qNO" or "qnorm", multivariate normal random effects
#' are first simulated, then transformed using the quantile function.
#'
#' @param randFxParms Ignored unless \code{family} is not "qNO" or "qnorm".
#' Parameters to be passed to a quantile function given in \code{family}. For
#' \code{\link{gamlss.family}} quantile functions, the parameters are mu,
#' sigma, nu, and tau.                    ,

}



