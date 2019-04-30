

#' samplingDist - function to plot the sampling distributions of
#' \code{\link{PersonAlytic}} output produce by \code{\link{ICTpower}}.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
samplingDist <- function(paout)
{
  # parse statNames
  statNames <- names(paout)[ grepl('statName', names(paout)) ]
  statsAre  <- list()
  for(i in seq_along(statNames))
  {
    statsAre[[statNames[i]]] <- statNameParse(paout[[statNames[i]]])
  }

  # plot statValues
  if(statsAre[[1]] == 'Correlations Between Y and Time')
  {
    statValues <- names(paout)[ grepl('statValue', names(paout)) ]
    hist(paout[[statValues[1]]], main=statsAre[[1]], xlab='Correlation')
    m <- paste('Mean =', round(mean(paout[[statValues[1]]], na.rm=TRUE),2))
    s <- paste('SD   =', round(  sd(paout[[statValues[1]]], na.rm=TRUE),2))
    legend("topright", legend = c(m, s))
  }

  # TODO delete, this distorts when unequal variance
  #xlim <- range( paout[, statValues[2:length(statValues)]], na.rm=TRUE )
  ##TODO - generalized list to any number of phases
  #cols <- c(rgb(1,0,0,0.5), rgb(0,0,1,0.5))
  #
  #for(i in seq_along(statsAre)[-1])
  #{
  #  hist(paout[[statValues[i]]], col=cols[i-1], xlim=xlim,
  #       add=ifelse(i==2, FALSE, TRUE), main='', xlab='Mean of Y')
  #  if(i==2)
  #  {
  #    title("Histogram of data means by phase")
  #  }
  #}
  #box()
  #legend("topright", legend=statsAre[2:length(statsAre)], col=cols, lty=1,
  #       lwd=6)


  # ggplot approach
  long <- list()
  wc   <- statValues[2:length(statValues)]
  for(i in seq_along(wc))
  {
    long[[i]] <- data.frame(Phase = i-1, Ymeans = paout[,wc[i]])
  }
  long <- do.call(rbind, long)
  long$Phase <- factor(long$Phase)
  ggplot(long, aes(x=Ymeans, fill=Phase)) + geom_density(alpha=0.25)


  ggplot(paout, aes(x=phase1.Value, fill=factor(phase1.p.value<=.05))) +
    geom_density(alpha=0.25) +
    guides(fill=guide_legend(title="P <= .05")) +
    xlab("Estimated Phase Effect")

  ggplot(paout, aes(x=Time.Value, fill=factor(Time.p.value<=.05))) +
    geom_density(alpha=0.25) +
    guides(fill=guide_legend(title="P <= .05")) +
    xlab("Estimated Time Effect")
}


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
#' # this example illustrates how mvrFam works under the hood
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
makeFam <- function(Y, parms, family="qNO",)
{
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
  # -- Note: rescaling not done here, instead we sim according to design$randFxCovMat
  # and choose the error variance commensurate to design$propErrVar

  #Y <- scale(Y, TRUE, TRUE) * sqrt(1-design$propErrVar)

  return(Y)
}



#' makeErr - simulate errors
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
ICTerror <- function(err           = armaErr$new() ,
                     n             = 100           ,
                     nObservations = 40            ,
                     seed          = 1             ,
                     ...                  )
{

  # get seeds
  seeds <- as.list( makeSeeds(seed, nObservations) )

  # sim errors, use sd = 1 for now, it can be rescaled in other functions
  FUN <- function(x, errorFUN, model, n, sd, ...)
  {
    set.seed(x)
    errorFUN(model = model, n = n, sd = sd, ...)
  }
  errors <- lapply(seeds, FUN,
                   errorFUN = err$FUN   ,
                   model    = err$parms ,
                   n        = n         ,
                   sd       = 1         ,
                   ...)
  errors <- data.frame(do.call(rbind, errors))

  # TODO move this to a validation test or examples for this function.
  # rescale to have unit variances within time points then multiply
  # by the proportion of error variance; need to transpose twice so
  # that the rescaling applies within person (which also gets the
  # between person scaling right, but not vice versa in initial tests)
  #errorsr <- t(scale(t(errors), FALSE, TRUE)) * sqrt(design$variances$errorVar)
#
  #doQC <- FALSE
  #if(doQC)
  #{
  #  # QC, only holds asymptotically
  #  # within persons:
  #  all.equal(mean(apply(errorsr, 1, var)), errorVar, tolerance = .05)
  #  # between persons:
  #  all.equal(mean(apply(errorsr, 2, var)), errorVar, tolerance = .05)
  #  tseries <- lapply(data.frame(t(errorsr)), ts)
  #  aFun <- function(x, order=c(1,0,0)) try(arima(x, order), silent = TRUE)
  #  sigma2s <- lapply(lapply(tseries, aFun), function(x) x$sigma2)
  #  # tends toward underestimation when nObservations are small, but very precise
  #  # when large
  #  hist(unlist(sigma2s))
  #}

  # return
  return(t(errors))
}


#' statNameParse - used by samplingDist()
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
statNameParse <- function(statName)
{

  # correlations
  if( all(grepl('correlation', statName)) )
  {
    statIs <- 'Correlations Between Y and Time'
  }

  # means
  if( all(grepl('mean', statName)) )
  {
    phase  <- strsplit(statName[1], "\\=")[[1]][2]
    statIs <- paste('Means of Phase =', phase)
  }

  # error
  if( ! exists('statIs') )
  {
    stop("The type of descriptive statistic cannot be determined.")
  }

  # return
  return(statIs)
}
