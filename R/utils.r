
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

#' makeDesignMat - create the fixed effects design matrix for n=1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
makeDesignMat <- function(phases      = makePhase() ,
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
    statsAre[[statNames[i]]] <- .statNameParse(paout[[statNames[i]]])
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

