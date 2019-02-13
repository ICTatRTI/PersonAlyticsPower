# TODO: save data for use in SAS

#' simICTdata - this is a user interface function for one dataset,
#' switching between different designs so the function can be updated
#' @author Stephen Tueller \email{stueller@@rti.org}

simICTdata <- function(design      = polyICT$new()           ,
                       randFxParms = list(mu=.5 , sigma=1)   ,
                       family      = "qNO"                   ,
                       randFxSeed  = 11                      ,
                       errorParms  = list(ar=c(.5), ma=c(0)) ,
                       errorFUN    = arima.sim               ,
                       errorSeed   = 21                      ,
                       designCheck = FALSE                   )
{
  # parallelize here and detect design from class(design); better yet,
  # we make polyICT, etc. methods for the class

  # do a design check with large N
  # TODO: move to a function
  if(designCheck)
  {
    # message
    message("\n\nStarting a design check with n=1000 participants\n",
            "WARNING: this check is currently only done for the standard MLM")

    # reset n temporarily, saving user specified n for later
    userN <- design$n
    design$n <- 1000

    # simulate random effects
    randFx <- mvrFam(design, randFxParms, family, randFxSeed)

    # simulate the errors
    errors <- ICTerror(design, errorParms, errorFUN, randFxSeed)

    # construct the data
    dat <- design$makeData(randFx, errors)

    # compare expected to observed variances
    expObsVar <- cbind( design$expectedVar(),
                        aggregate(dat$y, by=list(dat$Time), var)$x)

    # get the correlation of expected and observed variances
    expObsCor <- round(cor(expObsVar)[1,2],4)
    cat("The correlation between the expected variance and the observed\n",
        "variances is ", expObsCor)

    #
    ggplot(dat[dat$id<=100,], aes(x=Time, y=y, group=id, col=phase)) +
      geom_line() + geom_smooth(se = FALSE, size=.5) +
      ggtitle('Raw data and smoothed average trajectories for first 100 participants')

    # TODO: generalize the equation to the implied model, hmmm, need to generate
    # that from the inputs, currently only works for slopes model
    ctrl <- lmeControl(opt="optim")
    mod0 <- lme(y~phase*Time, data=dat, random = ~ Time | id,
                control=ctrl, correlation = corARMA(p=1,q=0))

    cat("\n\nMODEL RESULTS\n")
    print( round(rbind(summary(mod0)$tTable[,1]),3 ) )

    cat("\nMODEL INPUTS\n")
    print( c(unlist(design$effectSizes)) )

    cat("\n\n\n")
  }


  # parralelization




  #return(data)

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
  # simulate multivariate normal data
  set.seed(seed)
  Y <- mvrnorm( design$n, rep(0, design$polyOrder + 1), design$covMat)

  if(! family %in% c("qNO", "qnorm"))
  {
    # now get the propabilities
    YpNorm <- data.frame(pnorm(Y))

    # transform to the specified distribution
    Y <- doLapply(YpNorm, family, mu=parms$mu, sigma=parms$sigma,
                  nu=parms$nu, tau=parms$tau, ...)
  }

  # rescale the data - this may be problamtic if the user specifies a
  # discrete distribution for random effects, may add a warning here
  # Note: zero mean is specified hear b/c the mean off the random effects
  # is added later in the $makeData method of an ICTdesign child class. This
  # may not be the best approach, especially with non-normal data, but I don't
  # know whether the method of running pnorm through q* quantile functions
  # preserves means at all
  # -- Note: rescaling not done here, instead we sim according to design$covMat
  # and choose the error variance commensurate to design$propErrVar

  #Y <- scale(Y, TRUE, TRUE) * sqrt(1-design$propErrVar)

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

  # get seeds
  set.seed(seed)
  seeds <- as.list( ceiling(runif(design$n, 0, 9e6) ) )

  # sim errors, use sd = 1 for now, it can be rescaled in other functions
  FUN <- function(x, errorFUN, model, n, sd, ...)
  {
    set.seed(x)
    errorFUN(model = model, n = n, sd = sd, ...)
  }
  errors <- lapply(seeds, FUN,
                   errorFUN = errorFUN             ,
                   model    = errorParms           ,
                   n        = design$nObservations ,
                   sd       = 1                    ,
                   ...)
  errors <- data.frame(do.call(rbind, errors))

  # rescale to have unit variances within time points then multiply
  # by the proportion of error variance; need to transpose twice so
  # that the rescaling applies within person (which also gets the
  # between person scaling right, but not vice versa in initial tests)
  errorsr <- t(scale(t(errors), FALSE, TRUE)) * sqrt(design$variances$errorVar)

  doQC <- FALSE
  if(doQC)
  {
    # QC, only holds asymptotically
    # within persons:
    all.equal(mean(apply(errorsr, 1, var)), errorVar, tolerance = .05)
    # between persons:
    all.equal(mean(apply(errorsr, 2, var)), errorVar, tolerance = .05)
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


















