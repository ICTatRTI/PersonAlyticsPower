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
#' @param seed A random seed for replicating the power analysis.
#'
#' @param cores The number of cores used in parralelized simulation. The default
#' is to use one fewer cores than are detected on the computer. Do not exceed the
#' maximum available cores or unexpected results may occur and R may crash.
#'
#' @param savePowerReport Should the power report be saved using the \code{file}
#' name? If \code{TRUE}, a `txt` file is saved in the working directory. This
#' option is set to \code{FALSE} in the function \code{\link{ICTpowerSim}} which
#' instead saves the reports from several conditions in a `csv` file.
#'
#' @param ... Further arguments to be passed to \code{\link{PersonAlytic}} for
#' analysis. All options in \code{PersonAlytic} can be passed except for
#' \code{output}, \code{data}, \code{ids}, \code{dvs}, \code{time}, and
#' \code{phase} as these are taken from the \code{design}.
#'
#' @examples
#'
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
                     seed            = 123                       ,
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

  #
  # do a design check with large N
  #
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
    #
    # parralelization
    #

    # generate seeds
    seeds <- makeSeeds(seed, B)

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
      # construct the data
      dat <- design$makeData(seed)

      # rename for merging
      names(dat)[2] <- paste('y', b, sep='')

      # return
      return(dat)
    }
    cat("\n\n") # to separate progress bar from warnings and messages
    # stop the cluster
    parallel::stopCluster(cl)

    #
    # data clean up and save
    #
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

    #
    # process inputs in `...` that may be passed to `PersonAlytic`
    #
    # get the ARMA order from `correlation` if it is passed by ...
    if(!exists('correlation'))
    {
      ar          <- length(errorParms$ar[errorParms$ar!=0]) # TODO needs to come from design
      ma          <- length(errorParms$ar[errorParms$ma!=0])
      correlation <- paste('corARMA(p=', ar, ',q=', ma, ')', sep='')
      detectAR    <- FALSE
    }
    if(!exists('detectTO'))   detectTO   <- FALSE
    if(!exists('time_power')) time_power <- design$maxRandFx

    #
    # analyze the data using PersonAlytics, treating y1,...,yB
    # as separate outcomes
    #
    paout <- PersonAlytic(output      = file$file                        ,
                          data        = Data                             ,
                          ids         = 'id'                             ,
                          dvs         = as.list(paste('y', 1:B, sep='')) ,
                          time        = 'Time'                           ,
                          phase       = 'phase'                          ,
                          time_power  = time_power                       ,
                          correlation = correlation                      ,
                          detectAR    = detectAR                         ,
                          detectTO    = detectTO                         ,
                          ...
    )

    #
    # power analysis specific summary of results
    #
    powerL <- powerReport(paout, alpha, file=file$file,
                          saveReport=savePowerReport)


    #
    # distributions of the estimates
    #
    plotDists <- FALSE
    if(plotDists)
    {
      samplingDist(paout)
    }

    #
    return(powerL)

  }
}


