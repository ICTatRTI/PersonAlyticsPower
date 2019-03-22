################################################################################
# functions for setting up ICT designs
################################################################################

#' @name ICTsetup
#' @aliases makePhase
#' @aliases polyICTsetup
#' @aliases ICTviz
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @title functions for setting up ICT designs.
#'
#' @note \code{makePhase} is used to generate a `phase` variable.
#'
#' \code{polyICTsetup} is used to create a polynomial ICT design.
#'
#' \code{ICTviz} visualizes ICT sampling and resulting designs.
#'
#' @param nObsPerPhase Numeric vector. The length is the number of phases
#' and each entry is the number of observations per phase.
#'
#' @param phaseNames Numeric or character vector. The names of the phases,
#' e.g., phase 0 and phase 1, or phases "A", "B", and "C". Phase names must
#' be unique, i.e., even if your design is ABA, you should use 1, 2, 3, or
#' "A", "B", "C".
#'
#' @param designMatrix Matrix. See \link{\code{polyICTsetup}}.
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
#' @param DIST Character. A \link{\code{gamlss.family}} distribution.
#'
#' @param parms Named list. The parameters for \code{DIST}. \code{parms} should
#' be from length 1 to 4 with possible names including 'mu', 'sigma', 'nu', and
#' 'tau'.
#'
#' @param BLrange Numeric vector. The range of possible values at baseline. This
#' range should be smaller than, and within, the theoretical range of the
#' distribution specified by \code{DIST}.

#' @rdname ICTsetup
#' @export
#' @return A list of phase names replicated by the number of observations per
#' phase
#'
#' @examples
#' makePhase()
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

#' @rdname ICTsetup
#'
#' @export
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
#' This matrix can be constructed by the user, the polyICTsetup function simply
#' creates a shell for the user to fill in.
#'
#' @examples
#'
#' # view an example study design matrix that could have been initiated
#' # with `polyICTsetup`
#' g2exampleDesignMatrix
#'
#' # visualize this design using truncated standard normal
#' ICTviz(DIST = 'NO', parms=list(mu=0, sigma=1),
#'        designMatrix = g2exampleDesignMatrix)
#'
#' \dontrun{
#' # two phase example with 5 time points each
#' designMatrix <- polyICTsetup(makePhase(c(5,5), c("A","B")), nGroup=1)
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

polyICTsetup <- function(phases = makePhase(), nGroups = 1, effectSizes = NULL,
                       DIST = 'NO', parms = list(mu=0, sigma=1),
                       BLrange = c(-3, -2))
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

  designMatrix <- makeBounds(designMatrix, DIST, parms, BLrange)

  return(designMatrix)
}




#' @rdname ICTsetup
#'
#' @export
#'
#' @return If \code{DIST} and \code{parms} are given (and optionally,
#' \code{BLrange}), a density plot for the \code{DIST} is plotted. If
#' \code{BLrange} is given, lines indicating the truncated part of the distribution
#' from which baseline values will be sampled are shown.
#'
#' Additionally,
#' if \code{designMatrix} is provided, a visualization of the study desgin will
#' be plotted next to the density plot.
#'
#' @examples
#' # standard normal distribution with baseline values sampled from -3 to -2
#' ICTviz(DIST = 'NO', parms=list(mu=0, sigma=1), BLrange = c(-3, -2))
#'
#' ICTviz(DIST = 'NO', parms=list(mu=0, sigma=1),
#'        designMatrix = g2exampleDesignMatrix, BLrange = c(-2, -0))

ICTviz <- function(DIST = 'NO', parms=list(mu=0, sigma=1, nu=2, tau=2),
                   designMatrix=NULL, BLrange = NULL)
{
  # check inputs
  if( ! is(gamlss.dist::gamlss.family(DIST)) == "gamlss.family")
  {
    stop("The value of `DIST`=", DIST, " is not a `gamlss.family` distribution.")
  }

  # clean up parms for lazy people
  if(is.null(names(parms))) names(parms) <- c('mu', 'sigma', 'nu', 'tau')
  if(length(parms) > 4 |
     any( ! names(parms) %in%  c('mu', 'sigma', 'nu', 'tau') ))
  {
    stop("`parms` should be a named list of length 1 to 4 with possible names\n",
         "'mu', 'sigma', 'nu', or 'tau'.")
  }

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
  if(!is.null(BLrange))
  {
    if(length(BLrange) != 2) stop("`BLrange` must be a vector of length 2")
    if(! BLrange[2] > BLrange[1])
    {
      stop("The first value of `BLrange` must be smaller than the second value.\n",
           "You provided `BLrange` = c(", BLrange[1], ", ", BLrange[2], ")" )
    }
    if(BLrange[1]<xmin | BLrange[2]>xmax)
    {
      stop("`BLrange` = c(", BLrange[1], ", ", BLrange[2], ") is NOT within the ",
           "theoretical range of the ", family,
           "\ndistribution which is c(", round(xmin,2), ", ", round(xmax,2), ").")
    }
  }

  # get the density
  ddist  <- paste('d', DIST, sep='')
  dY <- doCall(ddist, x=x,
               mu=parms$mu,
               sigma=parms$sigma,
               nu=parms$nu, tau=parms$tau)
  dY <- data.frame(x=x, y=dY)
  g1 <- ggplot(dY, aes(x=x, y=y)) + geom_line() + ylab('Density') +
    xlab('Outcome')

  if(!is.null(BLrange))
  {
    g1 <- g1 +
      geom_vline(xintercept = BLrange[1], col='grey', linetype="dashed", size=2) +
      geom_vline(xintercept = BLrange[2], col='grey', linetype="dashed", size=2)
  }

  # if no desgin matrix is provided, print the plot
  if( is.null(designMatrix))
  {
    print( g1 + ggtitle(paste(family, pparms, sep=' with: ')) )
  }
  if(!is.null(BLrange) & is.null(designMatrix))
  {
    g11 <- g1 + ggtitle(paste(family, pparms, sep=' with: '),
                        "Dashed grey lines are the sampling bounds at baseline")
    print(g11)
  }

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

    top <- paste(
      paste(family, pparms, sep=' with: '),
      ifelse(!is.null(BLrange),
             "\nDashed grey lines on the left panel are the sampling bounds at baseline",
             "")
    )

    gridExtra::grid.arrange(g1, g2, nrow=1, ncol=2, widths = c(2,4), top=top)
  }

}

################################################################################
# internal utility functions below this line
################################################################################


#' makeBounds - function to get the bounds from the baseline bounds and
#' vectors of effect sizes
#'
#' Note that the bounds are approximate on the assumptiong that fixed width
#' bounds have sufficiently similar SD
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

makeBounds <- function(designMatrix, DIST, parms, BLrange)
{

  # define the truncated distribution
  trDIST <- gamlss.tr::gen.trun(par = BLrange,
                                family = NO, # TODO should be DIST, not taking it directly
                                type = "both")

  # find the trucated distributions mean and variance via approximate asymptotics
  y <- rNOtr(100000)  # TODO we need to refer to rNOtr dynamically
  segM  <- mean(y)
  segSD <- sd(y)

  # update the designMatrix
  designMatrixNew <- designMatrix
  meanCols <- names(designMatrix)[ grep('means_', names(designMatrix)) ]
  lowCols  <- names(designMatrix)[ grep('lower_', names(designMatrix)) ]
  uppCols  <- names(designMatrix)[ grep('upper_', names(designMatrix)) ]
  for(i in seq_along(meanCols))
  {
    designMatrixNew[[meanCols[i]]] <- designMatrix[[meanCols[i]]]*segSD + segM
    # and now we need to solve for the new segment bounds...here's an approximation
    # for now
    designMatrixNew[[lowCols[i]]] <- BLrange[1] + designMatrix[[meanCols[i]]]*segSD
    designMatrixNew[[uppCols[i]]] <- BLrange[2] + designMatrix[[meanCols[i]]]*segSD
  }
  return(designMatrixNew)
}




#' makeMeans - a function to create mean values given within-phase
#' growth model parameters (i, s, q, etc.; see \code{effectSizes} in
#' \link{\code{polyICT2}}).
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

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