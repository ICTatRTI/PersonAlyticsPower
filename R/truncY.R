

# try last answer from
# https://stats.stackexchange.com/questions/113230/generate-random-numbers-following-a-distribution-within-an-interval
# noting that this is very similar in structure to the approach taken in makeFam() for
# transforming data

# note, we use R.utils::doCall from R.utils to allow for ignore arguments (e.g., nu, tau)
truncY0 <- function(N, parms=list(mu=0, sigma=1, nu=2, tau=2),
                   DIST='NO', a = -Inf, b = Inf, seed=123)
{
  if(a > b) stop('`b` must be > `a`')

  # TODO this check is also in ICTviz, consider moving common checks to one function
  if( ! is(gamlss.dist::gamlss.family(DIST)) == "gamlss.family")
  {
    stop("The value of `DIST`=", DIST, " is not a `gamlss.family` distribution.")
  }

  pdist <- paste('p', DIST, sep='')
  qdist <- paste('q', DIST, sep='')

  set.seed(seed)
  p <- runif(N,
             R.utils::doCall(pdist, q=a,
                    mu=parms$mu,
                    sigma=parms$sigma,
                    nu=parms$nu, tau=parms$tau),
             R.utils::doCall(pdist, q=b,
                    mu=parms$mu,
                    sigma=parms$sigma,
                    nu=parms$nu, tau=parms$tau))

  R.utils::doCall(qdist, p=p,
         mu=parms$mu,
         sigma=parms$sigma,
         nu=parms$nu, tau=parms$tau)

}

# TODO this needs to be a method for a class, we're moving the same
# parameters around too much
# UL <- designMatrix[,4:5]
truncY <- function(N, UL, DIST='NO', parms=list(mu=0, sigma=1, nu=2, tau=2),
                   seed = 123)
{
  # this approach can't handle varying...so we'll have to loop or lapply?
  trCall <- paste("gamlss.tr::gen.trun(par = UL, family = ", DIST,
                  ", type = 'both', varying = TRUE)")
  eval(parse(text=trCall))
}


if(1==2)
{
  #dY <- data.frame(y=dY)
  #library(ggplot2)
  #ggplot(dY, aes(x=y)) + geom_density()
  # produces error
  #library(ggmap)
  #gglocator()

  rm3 <- function(x) round(mean(x), 3)
  pdf('truncYexample_Normal.pdf', height = 11, width = 8.5)
  par(mfrow=c(4,2))
  y0 <- truncY(10000, a=-Inf, b=Inf, seed=1);  hist(y0, main=paste('Mean=', rm3(y0))); qqnorm(y0); qqline(y0, col=2)
  y1 <- truncY(10000, a=-4,   b=-2,  seed=2);  hist(y1, main=paste('Mean=', rm3(y1))); qqnorm(y1); qqline(y1, col=2)
  y2 <- truncY(10000, a=-3,   b=-1,  seed=3);  hist(y2, main=paste('Mean=', rm3(y2))); qqnorm(y2); qqline(y2, col=2)
  y3 <- truncY(10000, a=-2,   b= 0,  seed=4);  hist(y3, main=paste('Mean=', rm3(y3))); qqnorm(y3); qqline(y3, col=2)
  y4 <- truncY(10000, a=-1,   b= 1,  seed=5);  hist(y4, main=paste('Mean=', rm3(y4))); qqnorm(y4); qqline(y4, col=2)
  y5 <- truncY(10000, a= 0,   b= 2,  seed=6);  hist(y5, main=paste('Mean=', rm3(y5))); qqnorm(y5); qqline(y5, col=2)
  y6 <- truncY(10000, a= 1,   b= 3,  seed=7);  hist(y6, main=paste('Mean=', rm3(y6))); qqnorm(y6); qqline(y6, col=2)
  y7 <- truncY(10000, a= 2,   b= 4,  seed=8);  hist(y7, main=paste('Mean=', rm3(y7))); qqnorm(y7); qqline(y7, col=2)
  dev.off()

  par(mfrow=c(1,2))

  pdf('truncYexample_Skew.pdf', height = 11, width = 8.5)
  par(mfrow=c(4,2))
  #y0 <- truncY(100, pdist='pSEP1',           , seed=1);  hist(y0); qqnorm(y0); qqline(y0, col=2)
  hist(rSEP1(1000, mu=0, sigma=1, nu=0, tau=.2))
  y1 <- truncY(1000, pdist='pSEP1', a=-4, b=-2, seed=2);  hist(y1); qqnorm(y1); qqline(y1, col=2)
  y2 <- truncY(1000, pdist='pSEP1', a=-3, b=-1, seed=3);  hist(y2); qqnorm(y2); qqline(y2, col=2)
  y3 <- truncY(1000, pdist='pSEP1', a=-2, b= 0, seed=4);  hist(y3); qqnorm(y3); qqline(y3, col=2)
  y4 <- truncY(1000, pdist='pSEP1', a=-1, b= 1, seed=5);  hist(y4); qqnorm(y4); qqline(y4, col=2)
  y5 <- truncY(1000, pdist='pSEP1', a= 0, b= 2, seed=6);  hist(y5); qqnorm(y5); qqline(y5, col=2)
  y6 <- truncY(1000, pdist='pSEP1', a= 1, b= 3, seed=7);  hist(y6); qqnorm(y6); qqline(y6, col=2)
  y7 <- truncY(1000, pdist='pSEP1', a= 2, b= 4, seed=8);  hist(y7); qqnorm(y7); qqline(y7, col=2)
  dev.off()

}

  ## get user inputs
  ## TODO: we need error handling and restarts here!!
  #
  #cat('\n')
  #nGroups <- readinteger("Enter the number of groups")
  #groupSizes <- list()
  #for(g in 1:nGroups)
  #{
  #  groupName <- readcharacter("Enter the name of group", g)
  #  groupSizes[[groupName]] <- readinteger("Enter the # of participants in group",
  #                                     groupName)
  #}
  #
  #nPhases <- readinteger("Enter the number of phases: ")
  #phaseLengths <- list()
  #for(p in 1:nPhases)
  #{
  #  phaseName <- readcharacter("Enter the name of phase", p)
  #  phaseLengths[[ phaseName ]] <- readinteger( "Enter the # of timepoints in phase",
  #                                              phaseName    )
  #}
  #
  #locator(type = 'p', col = 'red')




# 20190729 Material to be evaluated for deprecation ####


# #' @param designMatrix Matrix. See \code{\link{polyICTsetup}}.
# #'
# #' @param phases List. The length of phases is the number of phases. Each
# #' item in the list is a vector repeating the phase name or number as many
# #' times as there are time points within that phase. See \code{\link{makePhase}}.
# #'
# #' @param nGroups Integer. The number of groups.
# #'
# #' @param randFxMean See \code{polyICT2}. If \code{randFxMean} is
# #' provided, \code{nGroups} will be set to \code{length(randFxMean)}.
# #'
# #' @param fam Character. A \code{\link{gamlss.family}} distribution.
# #'
# #' @param parms Named list. The parameters for \code{fam}. \code{parms} should
# #' be from length 1 to 4 with possible names including 'mu', 'sigma', 'nu', and
# #' 'tau'.
# #'
# #' @param BLrange Numeric vector. The range of possible values at baseline. This
# #' range should be smaller than, and within, the theoretical range of the
# #' distribution specified by \code{fam}.

# #' @rdname ICTsetup
# #'
# #' @export
# #'
# #' @return
# #'
# #' For \code{polyICTsetup},
# #' a data.frame with a `time` and `phase` variable, plus the following pairs of
# #' variables each group:
# #'   \itemize{
# #'     \item lower_Group#. The lower bound at each time point.
# #'     \item upper_Group#. The upper bound at each time point.
# #'   }
# #' where the `#` is replaced by the group number. This naming convention is required
# #' for other functions in the package.
# #'
# #' Lower/upper bounds give the truncated range of the parent distribution
# #' from which observed values are sampled. ICTs often look at subpopulations whose
# #' scores lie on one end of the distribution, and a study intervetion aims to move
# #' them toward scores that are more normal.
# #'
# #' For example, in a three phase study with 5 time points in each phase. The first
# #' phase is a static baseline, the second phase occurs after an intervention which
# #' raises scores, and the third phase shows linear decay back to the original
# #' bounds after the removal of the intervention:
# #'
# #' \tabular{rrrr}{
# #' time \tab phase \tab lower_Group1 \tab upper_Group1 \cr
# #'    1 \tab     A \tab           -2 \tab            0 \cr
# #'    2 \tab     A \tab           -2 \tab            0 \cr
# #'    3 \tab     A \tab           -2 \tab            0 \cr
# #'    4 \tab     A \tab           -2 \tab            0 \cr
# #'    5 \tab     A \tab           -2 \tab            0 \cr
# #'    6 \tab     B \tab           -1 \tab            1 \cr
# #'    7 \tab     B \tab           -1 \tab            1 \cr
# #'    8 \tab     B \tab           -1 \tab            1 \cr
# #'    9 \tab     B \tab           -1 \tab            1 \cr
# #'   10 \tab     B \tab           -1 \tab            1 \cr
# #'   11 \tab     A \tab           -1 \tab            1 \cr
# #'   12 \tab     A \tab        -1.25 \tab         0.75 \cr
# #'   13 \tab     A \tab        -1.50 \tab         0.50 \cr
# #'   14 \tab     A \tab        -1.75 \tab         0.25 \cr
# #'   15 \tab     A \tab           -2 \tab            0
# #' }
# #'
# #' In this example, the modal value will be halfway between the lower and upper
# #' bound. The mean and median will be determined by the distribution selected
# #' by the user for data sampling and the truncated range of that distribution
# #' defined by the lower and upper bounds.
# #'
# #' This matrix can be constructed by the user, the polyICTsetup function simply
# #' creates a shell for the user to fill in.
# #'
# #' @examples
# #'
# #' # view an example study design matrix that could have been initiated
# #' # with `polyICTsetup`
# #' g2exampleDesignMatrix
# #'
# #' # visualize this design using truncated standard normal
# #' #ICTviz(fam = 'NO', parms=list(mu=0, sigma=1),
# #' #       designMatrix = g2exampleDesignMatrix)
# #'
# #' \dontrun{
# #' # two phase example with 5 time points each
# #' designMatrix <- polyICTsetup(makePhase(c(5,5), c("A","B")), nGroup=1)
# #'
# #' # view the empty matrix
# #' designMatrix
# #'
# #' # visualize an observed data distribution density plot
# #' ICTviz()
# #'
# #' # now the user needs to manually fill in the design matrix, values for
# #' # the lower/upper bounds are taken from the x-axis of the
# #' designMatrix <- edit(designMatrix)
# #' }
# #'
# #'  randFxMean = list(
# #'    group1 = list(
# #'      phase1 = c(i=0.0, s=0.0, q=0.0),
# #'      phase2 = c(i=0.2, s=0.0, q=0.0),
# #'      phase3 = c(i=0.0, s=0.0, q=0.0)
# #'    ),
# #'    group2 = list(
# #'      phase1 = c(i=0.0, s=0.0, q=0.0),
# #'      phase2 = c(i=0.0, s=0.0, q=0.0),
# #'      phase3 = c(i=0.0, s=0.0, q=0.0)
# #'    )
# #'    )
# #' designMatrix <- polyICTsetup(phases = makePhase(), randFxMean=randFxMean,
# #'   fam = 'NO')
# #' # TODO really need a portable method of moving items around, probably need to
# #' # make this an R6 class with methods.
# #' ICTviz(fam="NO", parms = list(mu=0, sigma=1), designMatrix=designMatrix,
# #' BLrange=c(-3,-2))
#
# polyICTsetup <- function(phases = makePhase(), nGroups = 1, randFxMean = NULL,
#                          fam = 'NO', parms = list(mu=0, sigma=1),
#                          BLrange = c(-3, -2))
# {
#   # TODO: add a check that phases and randFxMean conform, see checkPolyICT2()
#
#   # if randFxMean is provided, override nGroups
#   if(!is.null(randFxMean))
#   {
#     if(nGroups != length(randFxMean))
#     {
#       message('`randFxMeans` was provided, `nGroups` is being changed to ', nGroups)
#     }
#     nGroups <- length(randFxMean)
#     gNames  <- names(randFxMean)
#   }
#   if( is.null(randFxMean))
#   {
#     gNames <- paste('Group', 1:nGroups, sep='')
#   }
#
#   nPhases <- length(phases)
#   phases  <- unlist(phases)
#   designMatrix <- data.frame(time  = 0:(length(phases)-1),
#                              phase = phases )
#   for(g in 1:nGroups)
#   {
#     designMatrix[[paste('means_', gNames[g], sep='')]] <- as.numeric(NA)
#     designMatrix[[paste('lower_', gNames[g], sep='')]] <- as.numeric(NA)
#     designMatrix[[paste('upper_', gNames[g], sep='')]] <- as.numeric(NA)
#   }
#
#   # if randFxMean are provided, populate the mean_ columns of designMatrix
#   if(!is.null(randFxMean))
#   {
#     phaseNames <- unique(phases)
#     for(g in 1:nGroups)
#     {
#       wc <- paste('means_', gNames[g], sep='')
#       for(p in phaseNames)
#       {
#         temp <- designMatrix[[wc]][designMatrix$phase==p]
#         temp <- makeMeans(randFxMean[[gNames[g]]][[p]],  length(temp))
#         designMatrix[[wc]][designMatrix$phase==p] <- temp
#       }
#     }
#   }
#
#   designMatrix <- makeBounds(designMatrix, fam, parms, BLrange)
#
#   return(designMatrix)
# }
#
#
# #' @rdname ICTsetup
# #'
# #' @export
# #'
# #' @return For \code{ICTviz}, if \code{fam} and \code{parms} are given
# #' (and optionally,
# #' \code{BLrange}), a density plot for the \code{fam} is plotted. If
# #' \code{BLrange} is given, lines indicating the truncated part of the distribution
# #' from which baseline values will be sampled are shown.
# #'
# #' Additionally,
# #' if \code{designMatrix} is provided, a visualization of the study desgin will
# #' be plotted next to the density plot.
# #'
# #' @examples
# #' # standard normal distribution with baseline values sampled from -3 to -2
# #' ICTviz(fam = 'NO', parms=list(mu=0, sigma=1), BLrange = c(-3, -2))
# #'
# #' ICTviz(fam = 'NO', parms=list(mu=0, sigma=1),
# #'        designMatrix = g2exampleDesignMatrix, BLrange = c(-2, -0))
#
# ICTviz <- function(fam = 'NO', parms=list(mu=0, sigma=1, nu=2, tau=2),
#                    designMatrix=NULL, BLrange = NULL)
# {
#   .checkFam(fam, parms)
#
#   # set up items for plot title
#   faminfo <- gamlss.dist::gamlss.family(fam)
#   family  <- paste(faminfo$family[2], 'family')
#   wparms  <- names(parms) %in% names(faminfo$parameters)
#   pparms  <- paste(names(parms[wparms]), parms[wparms], sep="=")
#   pparms  <- paste(pparms, collapse = "; ")
#
#   # find the theoretical limits of the distribution given the parameters
#   qdist <- paste('q', fam, sep='')
#   xmin  <- R.utils::doCall(qdist, p=1/9e9)
#   xmax  <- R.utils::doCall(qdist, p=9e9/(9e9+1))
#   x     <- seq(xmin, xmax, length=1000)
#   if(!is.null(BLrange))
#   {
#     if(length(BLrange) != 2) stop("`BLrange` must be a vector of length 2")
#     if(! BLrange[2] > BLrange[1])
#     {
#       stop("The first value of `BLrange` must be smaller than the second value.\n",
#            "You provided `BLrange` = c(", BLrange[1], ", ", BLrange[2], ")" )
#     }
#     if(BLrange[1]<xmin | BLrange[2]>xmax)
#     {
#       stop("`BLrange` = c(", BLrange[1], ", ", BLrange[2], ") is NOT within the ",
#            "theoretical range of the ", family,
#            "\ndistribution which is c(", round(xmin,2), ", ", round(xmax,2), ").")
#     }
#   }
#
#   # get the density
#   ddist  <- paste('d', fam, sep='')
#   dY <- R.utils::doCall(ddist, x=x,
#                         mu=parms$mu,
#                         sigma=parms$sigma,
#                         nu=parms$nu, tau=parms$tau)
#   dY <- data.frame(x=x, y=dY)
#   g1 <- ggplot(dY, aes(x=x, y=y)) + geom_line() + ylab('Density') +
#     xlab('Outcome')
#
#   if(!is.null(BLrange))
#   {
#     g1 <- g1 +
#       geom_vline(xintercept = BLrange[1], col='grey', linetype="dashed", size=2) +
#       geom_vline(xintercept = BLrange[2], col='grey', linetype="dashed", size=2)
#   }
#
#   # if no desgin matrix is provided, print the plot
#   if( is.null(designMatrix))
#   {
#     print( g1 + ggtitle(paste(family, pparms, sep=' with: ')) )
#   }
#   if(!is.null(BLrange) & is.null(designMatrix))
#   {
#     g11 <- g1 + ggtitle(paste(family, pparms, sep=' with: '),
#                         "Dashed grey lines are the sampling bounds at baseline")
#     print(g11)
#   }
#
#   if(!is.null(designMatrix))
#   {
#     # get # groups
#     nGroups <- sum(grepl('group', tolower(names(designMatrix))))/2
#
#     # wide to long
#     # TODO: 'group' 'Group' - this cannot take arbitray group names!!
#     varying <- names(designMatrix)[ grep("er_group", tolower(names(designMatrix))) ]
#     times   <- expand.grid(c("lower", "upper"), paste("group", 1:nGroups, sep=""))
#     times   <- paste(times$Var2, times$Var1, sep = "_")
#     designMatrixL <- reshape(designMatrix, varying, v.names = "Y", direction = 'long',
#                              timevar = "group", times = times)
#     group_ul <- do.call(rbind, strsplit(designMatrixL$group, "_"))
#     designMatrixL$Group_ul <- designMatrixL$group
#     designMatrixL$ul <- group_ul[,2]
#     designMatrixL$Group <- group_ul[,1]
#
#     # find start points of phases
#     # TODO make this its own function if used elsewhere
#     phases <- designMatrixL$phase
#     nTimes <- max(designMatrixL$time)
#     wphases <- which(phases[1:(nTimes-1)] != phases[2:nTimes])
#     rects <- designMatrixL[wphases,]
#     rects <- rbind(rects, designMatrixL[nTimes,])
#     rects$xstart <- c(min(designMatrixL$time), rects$time[1:(nrow(rects)-1)])
#
#     # set up phase colors
#     cols <- RColorBrewer::brewer.pal(length(table(rects$phase))+1, 'Accent')
#     rects$cols <- NA; j <- 1
#     for(i in unique(rects$phase))
#     {
#       rects$cols[rects$phase==i] <- cols[j]; j <- j+1
#     }
#
#     # design plot
#     # TODO: names will be affected here too, cannot be arbitrary
#     g2 <- ggplot(designMatrixL, aes(x=time, y=Y, col=Group, group=Group_ul)) +
#       geom_line(aes(linetype=ul), size=2) +
#       geom_rect(data=rects, aes(xmin=xstart, xmax=time, ymin=-Inf, ymax=Inf, fill=phase),
#                 alpha=0.4, color = NA, inherit.aes = F) +
#       scale_fill_manual(values = rects$cols, name='phase') +
#       ylim(range(x)) + ylab('') +
#       scale_linetype_manual(values = c("solid", "twodash"),
#                             name = "bounds",
#                             guide = guide_legend(reverse = TRUE) ) +
#       theme(legend.key.width=unit(2,"line"))
#
#     g1 <- g1 + coord_flip() + scale_y_reverse()
#
#     top <- paste(
#       paste(family, pparms, sep=' with: '),
#       ifelse(!is.null(BLrange),
#              "\nDashed grey lines on the left panel are the sampling bounds at baseline",
#              "")
#     )
#
#     gridExtra::marrangeGrob(g1, g2, nrow=1, ncol=2, widths = c(2,4), top=top)
#   }
#
# }
#
# # internal utility functions below this line ###################################
#
#
# #' makeBounds - function to get the bounds from the baseline bounds and
# #' vectors of effect sizes
# #'
# #' Note that the bounds are approximate on the assumptiong that fixed width
# #' bounds have sufficiently similar SD
# #'
# #' @author Stephen Tueller \email{stueller@@rti.org}
# #'
# #' @keywords internal
# makeBounds <- function(designMatrix, fam, parms, BLrange)
# {
#   # check that fam is a gamlss.family distribution, esp. since eval(parse())
#   # used below
#   if( ! is(gamlss.dist::gamlss.family(fam)) == "gamlss.family")
#   {
#     stop("The value of `fam`=", fam, " is not a `gamlss.family` distribution.")
#   }
#
#   ### this code is preserved to prevent future developers from trying and failing
#   #   to use these approaches
#   # the following returns the function value, not the function
#   #famg <- get(fam)()
#   #eval(call(fam))
#   #famg <- eval(expression(fam))
#   #eval(quote(fam))
#   ### this fails b/c if fam='NO', famn still calls itself 'NO' internally which
#   #   fails to create famntr etc., instead producing NOtr etc.
#   # to make the names useable later, copy the distribution, density,
#   # distribution function, quantile function, and random generation functions
#   # to the internal distribution famn
#   #famn  <- eval(parse(text=fam))
#   #dfamn <- eval(parse(text=paste("d", fam, sep='')))
#   #pfamn <- eval(parse(text=paste("p", fam, sep='')))
#   #qfamn <- eval(parse(text=paste("q", fam, sep='')))
#   #rfamn <- eval(parse(text=paste("r", fam, sep='')))
#   #gamlss.tr::gen.trun(par = BLrange, family = "famn", type = 'both')
#
#   # eval-parse is frowned upon, but gamlss use of 'inherits' in 'as.gamlss.family'
#   # as used in 'gen.trun()' passes 'fam' instead of its value. To prevent
#   # unsafe code, there is a check that fam is a gamlss.family object. I tried
#   # the options presented here, but none worked
#   # https://stackoverflow.com/questions/1743698/evaluate-expression-given-as-a-string
#   trCall <- paste("gamlss.tr::gen.trun(par = BLrange, family = ", fam,
#                   ", type = 'both')")
#   eval(parse(text=trCall))
#   # now create R.utils::doCall -able names (just rfam for now)
#   # TODO - this should be passed so it doesn't need to be recreated!
#   # maybe not...we have several values for a,b across the study...
#   # see `varying` in ?gen.trun
#   rfamtr <- paste("r", fam, "tr", sep="")
#
#   # find the trucated distributions mean and variance via approximate asymptotics
#   set.seed(111) # TODO should we make this user controllable??
#   y     <- R.utils::doCall(rfamtr, n=100000)   # TODO we need to refer to rNOtr dynamically
#   segM  <- mean(y)
#   segSD <- sd(y)
#
#   # update the designMatrix
#   designMatrixNew <- designMatrix
#   meanCols <- names(designMatrix)[ grep('means_', names(designMatrix)) ]
#   lowCols  <- names(designMatrix)[ grep('lower_', names(designMatrix)) ]
#   uppCols  <- names(designMatrix)[ grep('upper_', names(designMatrix)) ]
#   for(i in seq_along(meanCols))
#   {
#     designMatrixNew[[meanCols[i]]] <- designMatrix[[meanCols[i]]]*segSD + segM
#     # and now we need to solve for the new segment bounds...here's an approximation
#     # for now
#     designMatrixNew[[lowCols[i]]] <- BLrange[1] + designMatrix[[meanCols[i]]]*segSD
#     designMatrixNew[[uppCols[i]]] <- BLrange[2] + designMatrix[[meanCols[i]]]*segSD
#   }
#   return(designMatrixNew)
#
#   # test simulation/this will be move to its own function/method
#   truncY
# }
#
# #' makeMeans - a function to create mean values given within-phase
# #' growth model parameters (i, s, q, etc.; see \code{randFxMean} in
# #' \code{polyICT2}).
# #'
# #' @author Stephen Tueller \email{stueller@@rti.org}
# #'
# #' @keywords internal
# makeMeans <- function(phaseParms=c(100, .5, -.25), nObs=10)
# {
#   # set up growth model factor loadings
#   times     <- 0:(nObs-1)
#   times     <- times %*% t( 0:(length(phaseParms)-1) )
#   times[,1] <- 1
#
#   # compute the means at each time point
#   means <- times * matrix(phaseParms, nrow(times), length(phaseParms), byrow = TRUE)
#   means <- apply(means, 1, sum)
#   return(means)
# }
