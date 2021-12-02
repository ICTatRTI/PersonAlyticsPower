#-------------------------------------------------------------------------------
# functions for setting up polyICT2 ####
#-------------------------------------------------------------------------------


#' Functions to help in setting up an ICT design
#'
#' \code{setupDist} and \code{fastDist} are used to check a range of parameters
#' for a non-normal residual distribution. \code{setupDist} requires the user to
#' specify a complete design (e.g., \code{\link{polyICT}} and residual
#' correlation struction (e.g., \code{\link{armaErr}}), but is very slow.
#' \code{fastDist} specifies a generic normal distribution (standing in for
#' the random effects) and an error distribution (with no residual autocorrelation)
#' which is less accurate but is much faster and can be used to ballpark inputs
#' for \code{setupDist}. See the examples.
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param design See \code{\link{designICT}} and \code{\link{polyICT}}.
#'
#' @param err,fam,famParms See \code{\link{armaErr}}. For \code{famParms},
#' a range of values can be specified using a numeric vector and
#' the resulting distributions from all possible combinations will be plotted
#' with the resulting
#'
#' @param propErrVar See \code{\link{polyICT}}.
#'
#' @param file A character string used to name output files. Extensions (e.g.,
#' pdf, csv) will be added later.
#'
#' @return For \code{fastDist} and \code{setupDist},
#' a pdf (with the prefix specified in \code{file}),
#' with density plots for normal data (standing in for the random effects
#' distribution), the error data (using the distribution in \code{fam} and its
#' parameters extracted from the ranges in \code{famParms}), and their
#' \code{propErrVar} weighted combined distribution. The plots are sorted
#' left to right according to the values of \code{propErrVar}, and rows according
#' to all possible combinations of the range of values in \code{famParms}. The
#' resulting skewness and kurtosis are printed in the upper right corner of
#' each plot.
#'
#' There is also a csv file (with the prefix specified in \code{file}) is
#' created tabulating all the setting and the resulting skewness and kurtosis
#' in the 'observed' distribution (i.e., the distribution you wish to simulate).
#'
#' @export
#' @examples
#'
#' \dontrun{
#' # Say we want to simulate data that has a skewness of 1.5. We'll first look
#' # at a large range of values passed to the Weibull (Type 3) distribution
#' # using `fastDist`.
#' fastDist(fam      = "WEI3"              ,
#'          famParms = list(mu    = 1:20,
#'                          sigma = c(.25,.5,.75,1)  ) ,
#'          propErrVar = c(.90, .75, .50 ) ,
#'          file = 'fastDist'
#'          )
#' }
#'
#' # After looking at the output, we decide to try a condition that had a
#' # skewness of 2.53 because we know we need a larger skewness from `fastDist`
#' # than our target of 1.5
#'
#' example(polyICT) # generate `myPolyICT` from ?polyICT
#' err <- armaErr$new(model = list(ar=c(.5), ma=c(.5)), fam = 'WEI3',
#'                     famParms = list(mu=1, sigma=.5))
#' myPolyICT$propErrVar <- c(randFx = .5, res    = .49, mserr  = .01)
#' myPolyICT$error <- err
#' datstat <- myPolyICT$designCheck(return='datstat', justData=TRUE, seed = 456)
#' datstat$descriptives

fastDist <- function(fam, famParms, propErrVar, file = 'fastDist.pdf')
{
  # n
  n <- 1e5

  # sort propErrVar
  propErrVar <- sort(propErrVar)
  if(! all(propErrVar > 0 & propErrVar < 1))
  {
    stop("\nAll values in `propErrVar` must be > 0 and < 1")
  }

  # set up a grid
  grids <- c(list(propErrVar=propErrVar), famParms)
  grids <- expand.grid(grids)

  # parralelization set up
  pkgs     <- c("gamlss", "nlme", "foreach", "PersonAlytics")
  capture.output( pb <- txtProgressBar(max = nrow(grids), style = 3),
                  file='NUL')
  progress <- function(n) setTxtProgressBar(pb, n)
  opts     <- list(progress = progress)
  cl       <- snow::makeCluster(spec=parallel::detectCores()-1,
                                type="SOCK", outfile="")
  snow::clusterExport(cl, c())
  doSNOW::registerDoSNOW(cl)

  # loop
  gg <- foreach( i=1:nrow(grids),  .packages = pkgs, .options.snow = opts) %dopar%
  {
    # parameters
    args <- list(n=n)
    for(j in 1:length(famParms))
    {
      parm <- names(famParms[j])
      args[[parm]] <- grids[[parm]][i]
    }

    # normal data
    d1 <- rNO(n)
    d1 <- data.frame(y = scale(d1), distribution = "Normal", alpha=.25)

    # error data
    d2 <- do.call(paste("r", fam, sep=""), args=args)
    d2 <- data.frame(y = scale(d2), distribution = as.gamlss.family(fam)$family[2],
                     alpha = .25)

    # observed data
    pN <- sqrt(1-grids$propErrVar[i])
    pE <- sqrt(  grids$propErrVar[i])
    d3 <- data.frame(y = pN*d1$y + pE*d2$y, distribution = "Observed", alpha=1)

    # plotting frame
    d <- rbind(d1,d2,d3)

    #
    descriptives <- dstats(d$y, d$distribution)
    dnms <- c("NO", as.gamlss.family(fam)$family[1], "Obs")
    ds   <- format(dnms, width=max(nchar(dnms)))
    ds   <- paste( paste(ds, ": skew=",
                     format(round(descriptives[1:3,4],2), nsmall=2),
                     ", kurt=",
                     format(round(descriptives[1:3,5],2), nsmall=2),
                     sep=''),
                  collapse = '\n')

    # plot for deciding on options
    g <-
    ggplot(d, aes(x=y, col=distribution, fill=distribution)) +
      geom_density(size=1.5, alpha=.2) +
      annotate("text", x=Inf, y=Inf, label = ds, vjust=1, hjust=1) +
      ggtitle(paste(paste(names(grids),
                          format(round(unlist(grids[i,]),2), nsmal = 2),
                          sep=' = '), collapse='; '))

    return( list(g=g, descriptives=descriptives) )
  }
  cat("\n\n") # to separate progress bar from warnings and messages
  # stop the cluster
  suppressWarnings( parallel::stopCluster(cl) )

  # extract plots from stats
  descriptives <- lapply(gg, function(x) x$descriptives[3,4:5])
  grids        <- data.frame(grids, do.call('rbind', descriptives))
  plots        <- lapply(gg, function(x) x$g)

  # set up pdf
  npev <- length(propErrVar)
  npages  <- ceiling( nrow(grids)/(3*npev) )
  grids$page <- sort(rep(1:npages, 9))[1:nrow(grids)]
  pdf(file = paste(file, 'pdf', sep='.'), width = 20, height=12)
  for(i in unique(grids$page))
  {
    do.call("grid.arrange", c(plots[which(grids$page==i)], ncol=npev))
  }
  dev.off()

  # save csv
  write.csv(grids, file=paste(file, 'csv', sep='.'))
}


#' @rdname fastDist
#' @export
#'
#' @examples
#'
#' # rerun the `fastDist` example with a subset of parameters to WEI3
#'
#' \dontrun{
#' # regenerate `myPolyICT` from ?polyICT
#' example(polyICT)
#'
#' # set up an error object, values for `famParms` are not needed and will be
#' # ignored by `setupDist`
#' skew <- armaErr$new(model = list(ar=c(.5), ma=c(.5)), fam = 'WEI3')
#'
#' # use a subset of distributions
#' setupDist(design = myPolyICT                    ,
#'           err    = skew                         ,
#'           famParms = list(mu    = c(1,3,5)   ,
#'                           sigma = c(.5, .75)    ,
#'           propErrVar = c(.90, .75, .50 )        ,
#'           file = 'setupDist'
#'           )
#' }
#'
setupDist <- function(design, err, famParms, propErrVar, file = 'setupDist')
{
  # set up a grid
  grids <- c(list(propErrVar=propErrVar), famParms)
  grids <- expand.grid(grids)

  # parralelization set up
  pkgs     <- c("gamlss", "nlme", "foreach", "PersonAlytics")
  capture.output( pb <- txtProgressBar(max = nrow(grids), style = 3),
                  file='NUL')
  progress <- function(n) setTxtProgressBar(pb, n)
  opts     <- list(progress = progress)
  cl       <- snow::makeCluster(spec=parallel::detectCores()-1,
                                type="SOCK", outfile="")
  snow::clusterExport(cl, c())
  doSNOW::registerDoSNOW(cl)

  # loop
  gg <- foreach( i=1:nrow(grids),  .packages = pkgs, .options.snow = opts) %dopar%
  {
    # clone the design
    d0 <- design$clone(deep=TRUE)

    # update err
    for(j in c('mu', 'sigma', 'nu', 'tau'))
    {
      if( j %in% names(grids) )
      {
        err$famParms[[j]] <- grids[[j]][i]
      }
    }
    #err$checkModel()

    # update design
    #d0$propErrVar <- c(randFx=.01, res=.98, mserr=.01)
    suppressMessages(
      d0$propErrVar <- c(randFx = 1 - (grids$propErrVar[[i]] + .01),
                         res    = grids$propErrVar[[i]],
                         mserr  = .01)
    )
    #d0$inputMat

    # add err to design
    d0$error <- err

    # run the design check
    datstat <- d0$designCheck(return='datstat', justData=TRUE, seed = i)
    #hist(datstat$data$y)
    #dstats(datstat$data$y, datstat$data$phase, print=TRUE)

    # plot for deciding on options
    g <- ggplot(datstat$data, aes(x=y, fill=phase)) + geom_density(alpha=0.25)
    ginfo <- ggplot_build(g)
    ypos <- max(ginfo$data[[1]]$density, na.rm=TRUE)
    xpos <- max(datstat$data$y, na.rm=TRUE)
    stats <- datstat$descriptives
    statnms <- colnames(stats)
    statnms <- gsub("\\s", " ", format(statnms, width=max(nchar(statnms))))
    stats <- format(round(stats, 2), nsmall = 2)
    stats <- paste(paste(statnms, stats, sep=' = '), collapse = '\n')

    g <- g + ggtitle(paste(paste(names(grids),
                                 format(round(unlist(grids[i,]),2), nsmal = 2),
                                 sep=' = '), collapse='; '))

    g <- g + annotate("text", x=Inf, y=Inf, label = stats, vjust=1, hjust=1)

    return( list(g=g, descriptives=descriptives) )
  }
  cat("\n\n") # to separate progress bar from warnings and messages
  # stop the cluster
  suppressWarnings( parallel::stopCluster(cl) )

  # extract plots from stats
  descriptives <- lapply(gg, function(x) x$descriptives[3,4:5])
  grids        <- data.frame(grids, do.call('rbind', descriptives))
  plots        <- lapply(gg, function(x) x$g)

  # set up pdf
  npev <- length(propErrVar)
  npages  <- ceiling( nrow(grids)/(3*npev) )
  grids$page <- sort(rep(1:npages, 9))[1:nrow(grids)]
  pdf(file = paste(file, '.pdf', sep='.'), width = 20, height=12)
  for(i in unique(grids$page))
  {
    do.call("grid.arrange", c(plots[which(grids$page==i)], ncol=npev))
  }
  dev.off()
}


#' makePhase - function to make phases
#'
#' @param nObsPerPhase Numeric vector. The length is the number of phases
#' and each entry is the number of observations per phase.
#'
#' @param phaseNames Numeric or character vector. The names of the phases,
#' e.g., \code{c("phase0", "phase1")}. Phase
#' names must each be unique, start with a letter, and contain only letters,
#' numbers, underscores, or decimals. The defualt is \code{NULL}, in which case
#' the \code{P} phases are labeled \code{c("phase1", "phase2", ..., "phaseP")}
#'
#' @export
#' @return For \code{makePhase}, a list of phase names replicated by the number
#' of observations per phase.
#'
#' @examples
#'
#' makePhase()
makePhase <- function(nObsPerPhase = c(10,20,10) ,
                      phaseNames = NULL  )
{
  if(any(duplicated(phaseNames)))
  {
    stop("`phaseNames` must be unique, see ?makePhase")
  }

  if(is.null(phaseNames))
  {
    phaseNames <- paste("phase", 1:length(nObsPerPhase), sep="")
  }

  phases <- list()
  for(p in seq_along(nObsPerPhase))
  {
    phases[[phaseNames[p]]] <- rep(phaseNames[p], nObsPerPhase[p])
  }
  return( phases )
}


