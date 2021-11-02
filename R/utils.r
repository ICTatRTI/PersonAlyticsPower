# utility functions ####

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

#' cutY - categorize x by a vector of proportions that sum to 1
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
cutY <- function(y, yCut)
{
  breaks <- c(-Inf, quantile(y, probs = cumsum(yCut)) )
  u <- cut(y, breaks = breaks)
  u <- as.numeric(factor(u, labels = 1:length(yCut))) - 1
  return(u)
}

# TODO this isn't a fully functional function, fixing it is a low priority as
# of 20190510
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

# NOTE: this is unfinished. It was to be used in the `inputMat` function within
# `.active()`, but the code currently there may be sufficient
#' .update - function to update multiple objects when any
#' one object is updated
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
.update <- function(..., anICT)
{
  if(length(ls())!=2)
  {
    stop("\nMore than two items were passed to `.update()`. Only 2 are allowed.\n\n")
  }
  if(exists("inputmat"))
  {
    if(!identical(anICT$inputMat, inputMat))
    {
      # update
    }
  }
  if(exists("groups"))
  {
    if(!identical(anICT$groups, groups))
    {

    }
  }
  if(exists("phases"))
  {
    if(!identical(anICT$phases, phases))
    {

    }
  }
  if(exits("designMat"))
  {
    if(!identical(anICT$designMat, designMat))
    {

    }
  }
}



















