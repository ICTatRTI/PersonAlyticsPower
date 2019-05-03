# TODO this file is a mix of internal and exported functions, split into a
# file containing internals and a file with exports

#' cor2cov - a version of this function was in the deprecated package `stremo`
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @param randFxCorMat Numeric matrix. An invertible correlation matrix.
#' @param variances Numeric vector. The variances.

cor2cov <- function(randFxCorMat    = matrix(c(1,.2,.2,1), 2, 2) ,
                    variances = c(1, .1)                   )
{
  sds <- sqrt(variances)
  b <- sds %*% t(sds)
  Sigma <- b * randFxCorMat
  # if QC is wanted use
  # all.equal(cov2cor(Sigma), randFxCorMat)
  return(Sigma)
}


#' checkCorMat
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'

checkCorMat <- function(randFxCorMat, cor=TRUE)
{
  mNm <- 'randFxCorMat'
  if(!cor) mNm <- 'varMat'
  if(!is.matrix(randFxCorMat)) stop('`randFxCorMat` must be a numeric matrix.')
  if( is.matrix(randFxCorMat))
  {
    isNum <- is.numeric(randFxCorMat)
    isSym <- isSymmetric(randFxCorMat)
    isCor <- all(diag(randFxCorMat)==1)
    isSlv <- class(try(solve(randFxCorMat), silent = TRUE)) %in% 'matrix'
    if(!isNum) stop('`', mNm, '` is not numeric')
    if(!isSym) stop('`', mNm, '` is not symmetric')
    if(cor & !isCor) stop('`randFxCorMat` does not have ones on the diagonal')
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
# TODO consider making this a method for ICTdesign that is passed by inheritence
# to polyICT, etc.
getICTdesign <- function(phases      = makePhase() ,
                         phaseNames  = NULL        ,
                         maxRandFx   = 2           ,
                         design      = 'polyICT'
)
{
  # make names if null
  if(is.null(phaseNames)) phaseNames <- paste('phase', 1:length(phases))

  # get the number of observations
  nObservations <- length(c(unlist((phases))))

  if(design == 'polyICT')
  {
    # generate the times
    times      <- list()
    times[['Time']] <- seq(0, nObservations-1, 1)
    if(maxRandFx>1)
    {
      for(i in 2:maxRandFx)
      {
        times[[paste('Time', i, sep='')]] <- times[['Time']]^i
      }
    }
    time    <- data.frame( do.call(cbind, times) )
    #varTime <- apply(time, 2, var)

    # clean up phases
    phase <- as.numeric(factor( c(unlist(phases)) ) ) - 1
    phase <- factor(phase, labels=phaseNames)

    designMat <- cbind(phase, time)

    return(designMat)
  }

}

#' totalVar - get the total variance for simulated data
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
# TODO consider making this a method for ICTdesign that is passed by inheritence
# to polyICT, etc.
totalVar <- function(randFxCovMat, propErrVar, designMat, randFxMean,
                     nObservations, n)
{

}

#' doIt - helper function for doLapply for ignoring extra arguments, written
#' to allow calls to \code{\link{gamlss.family}} distributions but ignore unused
#' parameters
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param X a vector of probabilities to be passed to a quantile function
#' like \code{\link{qN0}} in \code{\link{gamlss}} or \code{\link{qnorm}}.
#'
#' @param .fcn see \code{\link{R.utils::doCall}} in the `R.utils` package.
#'
#' @keywords internal

doIt <- function(X, .fcn="qNO", ...)
{
  R.utils::doCall(.fcn, p=X, ...)
}

#' doLapply - lapply with \code{\link{R.utils::doCall}} to ignore extra arguments to
#' \code{\link{gamlss.family}} distributions
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param Yp X a vector of probabilities to be passed to a quantile function
#'
#' @param .fcn A quantile function like \code{\link{qN0}} in \code{\link{gamlss}}
#' or \code{\link{qnorm}}. See also \code{\link{R.utils::doCall}} in the `R.utils` package
#'
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


#' checkFam - check whether \code{fam} is a gamlss.family object and
#' whether its parameters in \code{parms} are legitimate
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
checkFam <- function(fam, parms)
{
  # check family
  if( ! is(gamlss.dist::gamlss.family(fam)) == "gamlss.family")
  {
    stop("The value of `fam`=", fam, " is not a `gamlss.family` distribution.")
  }

  # clean up parms for lazy people -- may be superceded by err active bindings
  if(is.null(names(parms))) names(parms) <- c('mu', 'sigma', 'nu', 'tau')
  if(length(parms) > 4 |
     any( ! names(parms) %in%  c('mu', 'sigma', 'nu', 'tau') ))
  {
    stop("`parms` should be a named list of length 1 to 4 with possible names\n",
         "'mu', 'sigma', 'nu', or 'tau'.")
  }
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
    powerL[[i]] <- mean(paout[[i]] <= alpha, na.rm = TRUE)
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

#' seeds - make a vector of random seeds
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
makeSeeds <- function(seed, S)
{
  set.seed(seed)
  ceiling(runif(S, 0, 9e6))
}


#' randFxCorMatPop - correlation matrix populator, see \code{\link{checkPolyICT}}
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
randFxCorMatPop <- function(phaseNames, groupNames, randFxCorMat, wn=c('p', 'g', 'n'))
{
  randFxCorMatL <- list()
  for(p in seq_along(phaseNames))
  {
    for(g in seq_along(groupNames))
    {
      if(wn=='p') returnMat <- randFxCorMat[[p]]
      if(wn=='g') returnMat <- randFxCorMat[[g]]
      if(wn=='n') returnMat <- randFxCorMat
      randFxCorMatL[[phaseNames[p]]][[groupNames[g]]] <- returnMat
    }
  }
  return(randFxCorMatL)
}

#' randFxVarPop - random effects variance populator, see \code{\link{checkPolyICT}}
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
randFxVarPop <- function(phaseNames, groupNames, randFxVar, wn=c('p', 'g', 'n'))
{
  randFxVarL <- list()
  for(p in seq_along(phaseNames))
  {
    for(g in seq_along(groupNames))
    {
      if(wn=='p') returnMat <- randFxVar[[p]]
      if(wn=='g') returnMat <- randFxVar[[g]]
      if(wn=='n') returnMat <- randFxVar
      randFxVarL[[phaseNames[p]]][[groupNames[g]]] <- returnMat
    }
  }
  return(randFxVarL)
}

# TODO this isn't a fully functional function
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


