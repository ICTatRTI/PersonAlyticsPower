

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
#' @export

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
  sfile <- isRData <- iscsv <- NULL
  if(length(file)>1)
  {
    iscsv   <- file[2] %in% c('csv', 'CSV', 'Csv')
    isRData <- file[2] %in% c('RData', 'rdata', 'Rdata', 'RDATA')

    sfile <- paste(file[1],file[2], sep='.')

    hasPath <- grepl(glob2rx("*/*"), file[1])
    if(!hasPath) sfile <- paste(getwd(), sfile, sep='/')

    if(! iscsv & ! isRData)
    {
      stop("The file name file=\n\n'", file, "'\n\n is not a *.csv or an *.RData file.")
    }

    fileExists <- file.exists(sfile)

    if(fileExists)
    {
      message("The file\n\n", sfile, "\n\nexists and will be overwritten.")
    }
    if(!fileExists)
    {
      # test if the file can be created
      test=1
      if(iscsv  ) write.csv(test, sfile)
      if(isRData) save(test, file=sfile)
    }
  }
  invisible( list(file=file[1], isRData=isRData, iscsv=iscsv, sfile=sfile))
}


#' designCheck - check whether a file is writable and a csv or Rdata file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
designCheck <- function(design, file, family, randFxParms, randFxSeed,
                        errorParms, errorFUN, errorFamily, errorSeed)
{
  # message
  message("\n\nStarting a design check with n=1000 participants\n",
          "WARNING: this check is currently only done for the standard MLM")

  # save and reset n, due to inheritance design will get overwritten, fix
  # n below
  originaln <- design$n
  design$n  <- 1000

  # simulate random effects
  randFx <- mvrFam(design, randFxParms, family, randFxSeed)

  # simulate the errors
  errors <- ICTerror(design, errorParms, errorFUN, errorSeed)

  # construct the data
  dat <- design$makeData(randFx, errors)

  # get expected means and variances
  expectedVar <- design$expectedVariances

  # compare expected to observed variances
  expObsVar   <- cbind( expectedVar$Variances,
                        aggregate(dat$y, by=list(dat$Time), var)$x)
  # get the correlation of expected and observed variances
  expObsVcor <- round(cor(expObsVar)[1,2],4)
  cat("The correlation between the expected variances and the observed\n",
      "variances are ", expObsVcor, "\n\n")

  # compare expected to observed means
  expObsMean <- cbind( expectedVar$Means,
                       aggregate(dat$y, by=list(dat$Time), mean)$x)
  expObsMcor <- round(cor(expObsMean)[1,2],4)
  cat("The correlation between the expected means and the observed\n",
      "means are ", expObsMcor, "\n\n")

  # compare expected to observed total
  TotalMean <- mean(dat$y)
  TotalVar  <- var(dat$y)
  cat("The observed total mean is ", TotalMean,
      "\nThe expected total mean is ", expectedVar$TotalMean,
      "\nThe observed total variance is ", TotalVar,
      "\nThe expected total variance is ", expectedVar$TotalVar, "\n\n")

  # a plot
  g <- ggplot(dat[dat$id<=10,], aes(x=Time, y=y, group=id, col=phase)) +
    geom_line() + geom_smooth(se = FALSE, size=.5) +
    ggtitle('Raw data and smoothed average trajectories for first 10 participants')
  suppressMessages( print( g ) )

  # TODO: generalize the equation to the implied model, hmmm, need to generate
  # that from the inputs, currently only works for slopes model
  ctrl <- lmeControl(opt="optim")
  mod0 <- lme(y~phase*Time, data=dat, random = ~ Time | id,
              control=ctrl, correlation = corARMA(p=1,q=0))

  save(mod0, file=paste(file, 'designCheck.RData', sep='_'))
  print( mod0 )

  cat("\n\nMODEL RESULTS\n")
  print( round(rbind(summary(mod0)$tTable[,1]),3 ) )

  cat("\nMODEL INPUTS\n")
  print( c(unlist(design$effectSizes)) )

  cat("\n\n\n")

  design$n <- originaln
}


#' powerReport - print power results to screen and to a file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
powerReport <- function(paout, alpha, file, saveReport=TRUE)
{
  whichP <- names(paout)[ grepl('p.value', names(paout)) ]
  powerL <- list()
  for(i in whichP)
  {
    powerL[[i]] <- mean(paout[[i]] <= alpha)
  }

  # print the report to the screen
  names(powerL) <- gsub('.p.value', '', names(powerL))
  names(powerL) <- gsub("\\s", " ", format(names(powerL),
                                           width=max(nchar(names(powerL)))) )
  powerOutput <- paste(names(powerL), '\t', round(unlist(powerL),2), '\n' )
  cat( powerOutput )

  # save the report
  if(saveReport)
  {
   powerfile <- paste(file, 'PowerReport.txt', sep='_')
    cat( powerOutput, file = powerfile)
  }

  # return results
  return( unlist(powerL) )
}
