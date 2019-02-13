

#' cor2cov - a version of this function was in the deprecated package `stremo`
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

cor2cov <- function(corMat    = matrix(c(1,.2,.2,1), 2, 2) ,
                    variances = c(1, .1)                   )
{
  sds <- sqrt(variances)
  b <- sds %*% t(sds)
  Sigma <- b * corMat
  # if QC is wanted use
  # all.equal(cov2cor(Sigma), corMat)
  return(Sigma)
}

#' makePhase - a function to create the phase variable
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

makePhase <- function(nObsPerCondition = c(10,20,10) ,
                       conditions = c("A", "B", "A")  )
{
  phases <- list()
  for(i in seq_along(nObsPerCondition))
  {
    phases[[i]] <- rep(conditions[i], nObsPerCondition[i])
  }
  return( phases )
}


#' checkCorMat
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
checkCorMat <- function(corMat, cor=TRUE)
{
  mNm <- 'corMat'
  if(!cor) mNm <- 'varMat'
  if(!is.matrix(corMat)) stop('`corMat` must be a numeric matrix.')
  if( is.matrix(corMat))
  {
    isNum <- is.numeric(corMat)
    isSym <- isSymmetric(corMat)
    isCor <- all(diag(corMat)==1)
    isSlv <- class(try(solve(corMat), silent = TRUE)) %in% 'matrix'
    if(!isNum) stop('`', mNm, '` is not numeric')
    if(!isSym) stop('`', mNm, '` is not symmetric')
    if(cor & !isCor) stop('`corMat` does not have ones on the diagonal')
    if(!isSlv) stop('`', mNm, '` cannot be inverted')
  }
}

#' catMat - function to make a numeric matrix compatible with `cat`
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
catMat <- function(x=matrix(1:4, 2, 2))
{
  xList <- data.frame(t(x))
  xList <- paste(c('', rep(' ', length(xList)-1)), xList, '\n')
  xCat  <- paste(xList, collapse='')
  xCat
}


#' getICTdesign - create the fixed effects design matrix for n=1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal

getICTdesign <- function(phases      = makePhase() ,
                         polyOrder   = 2           ,
                         design      = 'polyICT'
)
{
  # get the number of observations
  nObservations <- length(c(unlist((phases))))

  if(design == 'polyICT')
  {
    # generate the times
    times      <- list()
    times[['Time']] <- seq(0, nObservations-1, 1)
    if(polyOrder>1)
    {
      for(i in 2:polyOrder)
      {
        times[[paste('Time', i, sep='')]] <- times[['Time']]^i
      }
    }
    time    <- data.frame( do.call(cbind, times) )
    #varTime <- apply(time, 2, var)

    # clean up phases
    phase <- as.numeric(factor( c(unlist(phases)) ) ) - 1

    # generate interactions
    phaseTime <- list()
    for(i in seq_along(time))
    {
      phaseTime[[paste('phase', names(time)[i], sep='')]] <- phase * time[,i]
    }
    phaseTime    <- data.frame( do.call(cbind, phaseTime) )
    #varPhaseTime <- apply(phaseTime, 2, var)

    designMat <- cbind(phase, time, phaseTime)

    return(designMat)
  }

}


#' doIt - helper function for doLapply for ignoring extra arguments, written
#' to allow calls to \code{\link{gamlss.family}} distributions but ignore unused
#' parameters
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param X a vector of probabilities to be based to a quantile function
#' like \code{\link{qN0}} in \code{\link{gamlss}} or \code{\link{qnorm}}
#' @param .fcn see \code{\link{doCall}} in the `R.utils` package
#'
#' @keywords internal
doIt <- function(X, .fcn="qNO", ...)
{
  doCall(.fcn, p=X, ...)
}

#' doLapply - lapply with \code{\link{doCall}} to ignore extra arguments to
#' \code{\link{gamlss.family}} distributions
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param Yp X a vector of probabilities to be based to a quantile function
#' @param .fcn A quantile function like \code{\link{qN0}} in \code{\link{gamlss}}
#' or \code{\link{qnorm}}. See also \code{\link{doCall}} in the `R.utils` package
#' @param ... other options passed to \code{\link{gamlss.family}} distribution
#' functions, most commonly mu, sigma, nu, and tau.
#'
#' @keywords internal
doLapply <- function(Yp, .fcn="qNO", ...)
{
  temp <- lapply(Yp, FUN=doIt, .fcn, ...)
  do.call(cbind, temp)
}


#' checkFile - check whether a file is writable and a csv or Rdata file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
checkFile <- function(file)
{
  file    <- tolower(file)
  isRData <- grepl('.rdata', file)
  iscsv   <- grepl('.csv', file)
  hasPath <- grepl(glob2rx("*/*"), file)

  if(!hasPath) file <- paste(getwd(), file, sep='/')

  if(!isRData & !iscsv)
  {
    stop("The file name file=\n\n'", file, "'\n\n is not a *.csv or an *.RData file.")
  }

  fileExists <- file.exists(file)

  if(fileExists)
  {
    message("The file\n\n", file, "\n\nexists and will be overwritten.")
  }
  if(!fileExists)
  {
    # test if the file can be created
    test=1
    if(iscsv)   write.csv(test, file)
    if(isRData) save(test, file=file)
  }
  invisible( list(file=file, isRData=isRData, iscsv=iscsv))
}


#' designCheck - check whether a file is writable and a csv or Rdata file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
designCheck <- function(design, family, randFxParms, randFxSeed,
                        errorParms, errorFUN, errorFamily, errorSeed)
{
  # message
  message("\n\nStarting a design check with n=1000 participants\n",
          "WARNING: this check is currently only done for the standard MLM")

  # reset n
  design$n <- 1000

  # simulate random effects
  randFx <- mvrFam(design, randFxParms, family, randFxSeed)

  # simulate the errors
  errors <- ICTerror(design, errorParms, errorFUN, errorSeed)

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
