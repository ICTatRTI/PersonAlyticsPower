

#' cor2cov - a version of this function was in the deprecated package `stremo`
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @param corMat Numeric matrix. An invertible correlation matrix.
#' @param variances Numeric vector. The variances.

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
#'
#' @param nObsPerPhase Numeric vector. The length is the number of phases
#' and each entry is the number of observations per phase.
#'
#' @param phaseNames Numeric or character vector. The names of the phases,
#' e.g., phase 0 and phase 1, or phases "A", "B", and "C". Phase names must
#' be unique, i.e., even if your design is ABA, you should use 1, 2, 3, or
#' "A", "B", "C".

makePhase <- function(nObsPerPhase = c(10,20,10) ,
                       phaseNames = c("phase1", "phase2", "phase3")  )
{
  if(any(duplicated(phaseNames)))
  {
    stop("`phaseNames` must be unique, see ?makePhase")
  }

  phases <- list()
  for(i in seq_along(nObsPerPhase))
  {
    phases[[i]] <- rep(phaseNames[i], nObsPerPhase[i])
  }
  return( phases )
}


#' studySetup - function to help users set up a study design matrix
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @param phases List. The length of phases is the number of phases. Each
#' item in the list is a vector repeating the phase name or number as many
#' times as there are time points within that phase. See \code{\link{makePhase}}.
#'
#' @param nGroups Integer. The number of groups.
#'
#' @param effectSizes See \link{\code{polyICT2}}. If \code{effectSizes} is
#' provided, \code{nGroups} will be set to \code{length(effectSizes)}.
#'
#' @return
#'
#' A data.frame with a `time` and `phase` variable, plus the following pairs of
#' variables each group:
#'   \itemize{
#'     \item lower_Group#. The lower bound at each time point.
#'     \item upper_Group#. The upper bound at each time point.
#'   }
#' where the `#` is replaced by the group number. This naming convention is required
#' for other functions in the package.
#'
#' Lower/upper bounds give the truncated range of the parent distribution
#' from which observed values are sampled. ICTs often look at subpopulations whose
#' scores lie on one end of the distribution, and a study intervetion aims to move
#' them toward scores that are more normal.
#'
#' For example, in a three phase study with 5 time points in each phase. The first
#' phase is a static baseline, the second phase occurs after an intervention which
#' raises scores, and the third phase shows linear decay back to the original
#' bounds after the removal of the intervention:
#'
#' \tabular{rrrr}{
#' time \tab phase \tab lower_Group1 \tab upper_Group1 \cr
#'    1 \tab     A \tab           -2 \tab            0 \cr
#'    2 \tab     A \tab           -2 \tab            0 \cr
#'    3 \tab     A \tab           -2 \tab            0 \cr
#'    4 \tab     A \tab           -2 \tab            0 \cr
#'    5 \tab     A \tab           -2 \tab            0 \cr
#'    6 \tab     B \tab           -1 \tab            1 \cr
#'    7 \tab     B \tab           -1 \tab            1 \cr
#'    8 \tab     B \tab           -1 \tab            1 \cr
#'    9 \tab     B \tab           -1 \tab            1 \cr
#'   10 \tab     B \tab           -1 \tab            1 \cr
#'   11 \tab     A \tab           -1 \tab            1 \cr
#'   12 \tab     A \tab        -1.25 \tab         0.75 \cr
#'   13 \tab     A \tab        -1.50 \tab         0.50 \cr
#'   14 \tab     A \tab        -1.75 \tab         0.25 \cr
#'   15 \tab     A \tab           -2 \tab            0
#' }
#'
#' In this example, the modal value will be halfway between the lower and upper
#' bound. The mean and median will be determined by the distribution selected
#' by the user for data sampling and the truncated range of that distribution
#' defined by the lower and upper bounds.
#'
#' This matrix can be constructed by the user, the studySetup function simply
#' creates a shell for the user to fill in.
#'
#' @examples
#'
#' # view an example study design matrix that could have been initiated
#' # with `studySetup`
#' g2exampleDesignMatrix
#'
#' # visualize this design using truncated standard normal
#' ICTviz(DIST = 'NO', parms=list(mu=0, sigma=1),
#'        designMatrix = g2exampleDesignMatrix)
#'
#' \dontrun{
#' # two phase example with 5 time points each
#' designMatrix <- studySetup(makePhase(c(5,5), c("A","B")), nGroup=1)
#'
#' # view the empty matrix
#' designMatrix
#'
#' # visualize an observed data distribution density plot
#' ICTviz()
#'
#' # now the user needs to manually fill in the design matrix, values for
#' # the lower/upper bounds are taken from the x-axis of the
#' designMatrix <- edit(designMatrix)
#' }

studySetup <- function(phases = makePhase(), nGroups = 1, effectSizes = NULL)
{
  # TODO: add a check that phases and effectSizes conform, see checkPolyICT2()

  # if effectSizes is provided, override nGroups
  if(!is.null(effectSizes))
  {
    if(nGroups != length(effectSizes))
    {
      message('`effectSizess` was provided, `nGroups` is being changed to ', nGroups)
    }
    nGroups <- length(effectSizes)
    gNames <- names(effectSizes)
  }
  if( is.null(effectSizes))
  {
    gNames <- paste('Group', 1:nGroups, sep='')
  }

  nPhases <- length(phases)
  phases <- unlist(phases)
  designMatrix <- data.frame(time  = 0:(length(phases)-1),
                             phase = phases )
  for(g in 1:nGroups)
  {
    designMatrix[[paste('means_', gNames[g], sep='')]] <- as.numeric(NA)
    designMatrix[[paste('lower_', gNames[g], sep='')]] <- as.numeric(NA)
    designMatrix[[paste('upper_', gNames[g], sep='')]] <- as.numeric(NA)
  }

  # if effectSizes are provided, populate the mean_ columns of designMatrix
  if(!is.null(effectSizes))
  {
    phaseNames <- unique(phases)
    for(g in 1:nGroups)
    {
      wc <- paste('means_', gNames[g], sep='')
      for(p in phaseNames)
      {
        temp <- designMatrix[[wc]][designMatrix$phase==p]
        temp <- makeMeans(effectSizes[[gNames[g]]][[p]],  length(temp))
        designMatrix[[wc]][designMatrix$phase==p] <- temp
      }
    }
  }

  return(designMatrix)
}

#' makeMeans - a function to create mean values given within-phase
#' growth model parameters (i, s, q, etc.; see \code{effectSizes} in
#' \link{\code{polyICT2}}).
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
makeMeans <- function(phaseParms=c(100, .5, -.25), nObs=10)
{
  # set up growth model factor loadings
  times     <- 0:(nObs-1)
  times     <- times %*% t( 0:(length(phaseParms)-1) )
  times[,1] <- 1

  # compute the means at each time point
  means <- times * matrix(phaseParms, nrow(times), length(phaseParms), byrow = TRUE)
  means <- apply(means, 1, sum)
  return(means)
}


#' ICTviz - a function to create the phase variable
#' and select mean values
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export

ICTviz <- function(DIST = 'NO', parms=list(mu=0, sigma=1, nu=2, tau=2),
                   designMatrix=NULL)
{
  # clean up parms for lazy people
  if(is.null(names(parms))) names(parms) <- c('mu', 'sigma', 'nu', 'tau')

  # set up items for plot title
  faminfo <- gamlss.dist::gamlss.family(DIST)
  family  <- paste(faminfo$family[2], 'family')
  wparms  <- names(parms) %in% names(faminfo$parameters)
  pparms  <- paste(names(parms[wparms]), parms[wparms], sep="=")
  pparms  <- paste(pparms, collapse = "; ")

  # find the theoretical limits of the distribution given the parameters
  qdist <- paste('q', DIST, sep='')
  xmin  <- doCall(qdist, p=1/9e9)
  xmax  <- doCall(qdist, p=9e9/(9e9+1))
  x     <- seq(xmin, xmax, length=1000)

  # get the density
  ddist  <- paste('d', DIST, sep='')
  dY <- doCall(ddist, x=x,
               mu=parms$mu,
               sigma=parms$sigma,
               nu=parms$nu, tau=parms$tau)
  dY <- data.frame(x=x, y=dY)
  g1 <- ggplot(dY, aes(x=x, y=y)) + geom_line() + ylab('Density') +
    xlab('Outcome') + ggtitle(paste(family, pparms, sep=" with: "))

  if( is.null(designMatrix)) print( g1 )

  if(!is.null(designMatrix))
  {
    # get # groups
    nGroups <- sum(grepl('Group', names(designMatrix))/2)

    # wide to long
    varying <- names(designMatrix)[ grep("Group", names(designMatrix)) ]
    times   <- expand.grid(c("lower", "upper"), paste("Group", 1:nGroups, sep=""))
    times   <- paste(times$Var2, times$Var1, sep = "_")
    designMatrixL <- reshape(designMatrix, varying, v.names = "Y", direction = 'long',
                             timevar = "Group", times = times)
    group_ul <- do.call(rbind, strsplit(designMatrixL$Group, "_"))
    designMatrixL$Group_ul <- designMatrixL$Group
    designMatrixL$ul <- group_ul[,2]
    designMatrixL$Group <- group_ul[,1]

    # find start points of phases
    # TODO make this its own function if used elsewhere
    phases <- designMatrixL$phase
    nTimes <- max(designMatrixL$time)
    wphases <- which(phases[1:(nTimes-1)] != phases[2:nTimes])
    rects <- designMatrixL[wphases,]
    rects <- rbind(rects, designMatrixL[nTimes,])
    rects$xstart <- c(min(designMatrixL$time), rects$time[1:(nrow(rects)-1)])

    # set up phase colors
    cols <- RColorBrewer::brewer.pal(length(table(rects$phase))+1, 'Accent')
    rects$cols <- NA; j <- 1
    for(i in unique(rects$phase))
    {
      rects$cols[rects$phase==i] <- cols[j]; j <- j+1
    }

    # design plot
    g2 <- ggplot(designMatrixL, aes(x=time, y=Y, col=Group, group=Group_ul)) +
      geom_line(aes(linetype=ul), size=2) +
      geom_rect(data=rects, aes(xmin=xstart, xmax=time, ymin=-Inf, ymax=Inf, fill=phase),
                alpha=0.4, color = NA, inherit.aes = F) +
      scale_fill_manual(values = rects$cols, name='phase') +
      ylim(range(x)) + ylab('') +
      scale_linetype_manual(values = c("solid", "twodash"),
                            name = "bounds",
                            guide = guide_legend(reverse = TRUE) ) +
      theme(legend.key.width=unit(2,"line"))

    g1 <- g1 + coord_flip() + scale_y_reverse()

    gridExtra::grid.arrange(g1, g2, nrow=1, ncol=2, widths = c(1,4))
  }

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
#' @param X a vector of probabilities to be passed to a quantile function
#' like \code{\link{qN0}} in \code{\link{gamlss}} or \code{\link{qnorm}}.
#'
#' @param .fcn see \code{\link{doCall}} in the `R.utils` package.
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
#' @param Yp X a vector of probabilities to be passed to a quantile function
#'
#' @param .fcn A quantile function like \code{\link{qN0}} in \code{\link{gamlss}}
#' or \code{\link{qnorm}}. See also \code{\link{doCall}} in the `R.utils` package
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


#' designCheck - check whether a file is writable and a csv or Rdata file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
designCheck <- function(design, file, family, randFxParms, randFxSeed,
                        errorParms, errorFUN, errorFamily, errorSeed)
{
  # message
  message("\n\nStarting a design check with n=5000 participants\n",
          "WARNING: this check is currently only done for the standard MLM")

  # save and reset n, due to inheritance design will get overwritten, fix
  # n below
  originaln <- design$n
  design$n  <- 5000

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
  # that from the inputs, currently only works for slopes model with AR(1)
  pa   <- Palytic$new(data=dat, ids='id', dv='y', time='Time', phase='phase')
  mod0 <- pa$lme()

  save(mod0, file=paste(file[1], 'designCheck.RData', sep='_'))
  print( mod0 )

  #TODO this only works for linear models
  #TODO needs better matching, force name matching in polyICT
  Estimates <- round(rbind(summary(mod0)$tTable[,1]),3 )
  Inputs    <- round(c(unlist(design$unStdEffects))[c(1,3,2,4)],3)
  cat('\n\nCheck the effect size estimates against inputs:\n')
  print( data.frame(Inputs=Inputs, Estimates=t(Estimates)) )

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
