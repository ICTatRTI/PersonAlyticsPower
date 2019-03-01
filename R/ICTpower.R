#' ICTpower - get simulated power for an ICT design using parametric bootstrap.
#'
#' @export
#' @import gamlss
#' @import nlme
#' @import foreach
#' @import snow
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param file The file name for saving power analysis results.
#' It may be a full path. The simulated data are to be saved by including the
#' file extension in a list. For example,
#' \code{file=c('n10_medSlopeEffect', 'csv')}
#' or \code{file=c('condition2', 'Rdata')}. Only `csv` and `RData` are supported.
#'
#' @param design An \code{\link{ICTdesign}} object such as
#' \code{\link{polyICT}}.
#'
#' @param B The number of simulated dataset (or parametric bootstrap replications).
#'
#' @param checkDesign Logical. Default is \code{FALSE}. If \code{TRUE}, an
#' initial sumulated dataset with n=1,000 participants will be simulated and the
#' model fit to the data. This will show whether the design inputs in \code{design}
#' are recovered under large sample. The data for the 1st 10 participants will
#' also be plotted. This allows users to check their inputs.
#'
#' @param alpha Numeric. Default is .05. The Type I error rate for computing
#' emprical power (the proportion of p-values across \code{B} data sets that are
#' less than or equal to \code{alpha}.
#'
#' @param randFxFamily A quantile function such as those in
#' \code{\link{gamlss.family}}.
#' If \code{family} is not "qNO" or "qnorm", multivariate normal random effects
#' are first simulated, then transformed using the quantile function.
#'
#' @param randFxParms Ignored unless \code{family} is not "qNO" or "qnorm".
#' Parameters to be passed to a quantile function given in \code{family}. For
#' \code{\link{gamlss.family}} quantile functions, the parameters are mu,
#' sigma, nu, and tau.
#'
#' @param randFxSeed A random seed for generating random effects.
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
#' @param errorseed A random seed for generating errors.
#'
#' @param cores The number of cores used in parralelized simulation. The default
#' is to use one few cores than are detected on the computer. Do not exceed the
#' maximum available cores or unexpect results may occur or R may crash.
#'
#' @param savePowerReport Should the power report be saved using the \code{file}
#' name? If \code{TRUE}, a `txt` file is saved in the working directory. This
#' option is set to \code{FALSE} in the function \code{\link{ICTpowerSim}} which
#' instead saves the reports from several conditions in a `csv` file.
#'
#' @param ... Further arguments to be passed to \code{\link{PersonAlytics}}.
#'
#' @examples
#'
#' # TODO change to dontrun
#' \dontrun{
#' ICTpower(c('testICTpower10', 'csv'),
#'   polyICT$new(n=10),
#'   randFxSeed = 54, errorSeed=89496,
#'   checkDesign='only')
#' ICTpower(c('testICTpower10', 'csv'),
#'   polyICT$new(n=200),
#'   randFxSeed = 23, errorSeed=923,
#'   checkDesign='only')
#' ICTpower(c('testICTpower20_2', 'csv'), polyICT$new(randFxVar=c(5,2)),
#'   B=10, checkDesign='only', randFxSeed=20342, errorSeed=4202)
#' }


ICTpower <- function(file                                        ,
                     design          = polyICT$new()             ,
                     B               = 100                       ,
                     checkDesign     = 'no'                      ,
                     alpha           = .05                       ,
                     randFxFamily    = "qNO"                     ,
                     randFxParms     = list(mu=.5 , sigma=1)     ,
                     randFxSeed      = 11                        ,
                     errorParms      = list(ar=c(.5), ma=c(0))   ,
                     errorFUN        = arima.sim                 ,
                     errorFamily     = NULL                      ,
                     errorSeed       = 21                        ,
                     cores           = parallel::detectCores()-1 ,
                     savePowerReport = TRUE                      ,
                     ...
)
{
  #argList <-  as.list(match.call(expand.dots = TRUE)[-1]) # only explicitly passed
  argList <- mget(names(formals()),sys.frame(sys.nframe()))
  print(argList)

  # check file name
  if(!is.null(file)) file <- checkFile(file)

  #...........................................................................
  # do a design check with large N
  #...........................................................................
  if(! checkDesign %in% c('no', 'yes', 'only'))
  {
    stop("`checkDesign` must be one of 'no', 'yes', or 'only',\n",
         "where 'only' only checks the design but does not simulate data.\n",
         "Instead, `checkDesign` was given as ", checkDesign, '.\n\n')
  }
  if(checkDesign %in% c('yes', 'only'))
  {
    designCheck(design      = design       ,
                file        = file$file    ,
                family      = randFxFamily ,
                randFxParms = randFxParms  ,
                randFxSeed  = randFxSeed   ,
                errorParms  = errorParms   ,
                errorFUN    = errorFUN     ,
                errorFamily = errorFamily  ,
                errorSeed   = errorSeed
                )
  }

  if(checkDesign != 'only')
  {
    #...........................................................................
    # parralelization
    #...........................................................................

    # generate seeds
    set.seed(randFxSeed)
    randFxSeeds <- ceiling(runif(B, 0, 9e6) )
    set.seed(errorSeed)
    errorSeeds <- ceiling(runif(B, 0, 9e6) )

    # message and timing
    message("\nStarting simulation of B=", B, " data sets.\n",
            ifelse(!is.null(file$sfile),
                   paste("Data will be saved in the file:\n\n", file$sfile,
                         "\n\nwhere the outcome for each replicate will be",
                         " labeled y1, ..., y", B,
                         "\n\n", sep=""),
                   '\n\n')
    )
    DIM <- 1:B
    start <- Sys.time()

    # parralelization set up
    pkgs     <- c("gamlss", "nlme", "foreach")
    capture.output( pb <- txtProgressBar(max = length(DIM), style = 3),
                    file='NUL')
	  progress <- function(n) setTxtProgressBar(pb, n)
    opts     <- list(progress = progress)
    cl       <- snow::makeCluster(cores, type="SOCK", outfile="")
    snow::clusterExport(cl, c())
    doSNOW::registerDoSNOW(cl)

    Data <- foreach( b=DIM, .packages = pkgs, .options.snow = opts) %dopar%
    {

      # simulate random effects
      randFx <- mvrFam(design, randFxParms, randFxFamily, randFxSeeds[b])

      # simulate the errors
      errors <- ICTerror(design, errorParms, errorFUN, errorSeeds[b])

      # construct the data
      dat <- design$makeData(randFx, errors)

      # rename for merging
      names(dat)[2] <- paste('y', b, sep='')

      # return
      return(dat)
    }
    cat("\n\n") # to separate progress bar from warnings and messages
    # stop the cluster
    parallel::stopCluster(cl)

    #...........................................................................
    # data clean up and save
    #...........................................................................
    # merge data in Data
    mergeby <- names(Data[[1]])
    mergeby <- mergeby[mergeby != 'y1']
    Data    <- Reduce(function(df1, df2) merge(df1, df2, by=mergeby, all = TRUE),
                      Data)

    # if requested, save the data
    if(!is.null(file$sfile))
    {
      if(file$isRData)
      {
        save(Data, file=file$sfile)
      }
      if(file$iscsv)
      {
        write.csv(Data, file=file$sfile, row.names = FALSE)
      }
    }

    #...........................................................................
    # process inputs in `...` that may be passed to `PersonAlytic`
    #...........................................................................
    # get the ARMA order from
    if(!exists('correlation'))
    {
      ar <- length(errorParms$ar[errorParms$ar!=0])
      ma <- length(errorParms$ar[errorParms$ma!=0])
      correlation <- paste('corARMA(p=', ar, ',q=', ma, ')', sep='')
      detectAR <- FALSE
    }
    if(!exists('detectTO')) detectTO <- FALSE

    #...........................................................................
    # analyze the data using PersonAlytics, treating y1,...,yB
    # as separate outcomes
    #...........................................................................
    paout <- PersonAlytic(output      = file$file                        ,
                          data        = Data                             ,
                          ids         = 'id'                             ,
                          dvs         = as.list(paste('y', 1:B, sep='')) ,
                          time        = 'Time'                           ,
                          phase       = 'phase'                          ,
                          time_power  = design$polyOrder                 ,
                          correlation = correlation                      ,
                          detectAR    = detectAR                         ,
                          detectTO    = detectTO                         ,
                          ...
    )

    #...........................................................................
    # power analysis specific summary of results
    #...........................................................................
    powerL <- powerReport(paout, alpha, file=file$file,
                          saveReport=savePowerReport)


    #...........................................................................
    # distributions of the estimates
    #...........................................................................
    plotDists <- FALSE
    if(plotDists)
    {
      samplingDist(paout)
    }

    #...........................................................................
    return(powerL)

  }
}


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



#' mvrFam - a function to simulate multivariate data from a gamlss family
#' distribution, following
#' "http://www.econometricsbysimulation.com/2014/02/easily-generate-correlated-variables.html"
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
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
#' # this example illustrates how mvrFam works under the hood
#' \dontrun{
#' # requires library(psych), not required by PersonAlyticsPower
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
#' psych::pairs.panels(Y)
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


















