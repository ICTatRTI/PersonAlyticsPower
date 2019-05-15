#' ICTpower - get simulated power for an ICT design using parametric or
#' non-parametric bootstrap.
#'
#' @export
#' @import gamlss
#' @import nlme
#' @import foreach
#' @import snow
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param outFile The file name for saving power analysis results.
#' It may be a full path. The simulated data are to be saved by including the
#' file extension in a list. For example,
#' \code{file=c('n10_medSlopeEffect', 'csv')}
#' or \code{file=c('condition2', 'Rdata')}. Only `csv` and `RData` are supported.
#'
#' @param design An \code{\link{designICT}} object such as
#' \code{\link{polyICT}}. This must be provided for a parametric bootstrap
#' esstimate of power.
#'
#' @param B The number of simulated dataset (or parametric bootstrap replications).
#'
#' @param dataFile Character. A file name for exntant data. Use this option to
#' do a non-parametric bootstrap (e.g., for finite sample power).
#' If \code{dataFile} is provide, \code{design}
#' will be ignored. The \code{dataFile} must be a *.csv file with the variables
#' 'id' and 'Time', optionally 'phase' and/or 'group', and dependent variables
#' labeled 'y'. Alternatively, \code{dataFile} may be a *.Rdata file named 'Data'
#' with the same columns as described for the *.csv file.
#'
#' @param sampleSizes Numeric vector of the same length as the number of groups
#' in the data. The group specific sample sizes to be drawn \code{B} times from
#' the \code{dataFile} data.
#'
#' @param alpha Numeric. Default is .05. The Type I error rate for computing
#' emprical power (i.e., the proportion of p-values across \code{B} data sets
#' that are less than or equal to \code{alpha}).
#'
#' @param seed A random seed for replicating the power analysis.
#'
#' @param cores The number of cores used in parralelized simulation (for a
#' parametric bootstrap) and for fitting models to the data for both bootstrap
#' types (see \code{\link{PersonAlytic}}). The default
#' is to use one fewer cores than are detected on the computer. Do not exceed
#' the maximum available cores or unexpected results may occur and R may crash.
#'
#' @param savePowerReport Should the power report be saved using the
#' \code{outFile} name? If \code{TRUE}, a `txt` file is saved in the working
#' directory. This option is set to \code{FALSE} in the function
#' \code{\link{ICTpowerSim}} which instead saves the reports from several
#' conditions in a `csv` file.
#'
#' @param standardize Named list of length 3. Default is
#' \code{list(dv = TRUE, ivs = FALSE, byids = TRUE)}. See
#' \code{\link{PersonAlytic}} or \code{\link{Palytic}}.
#'
#' @param ... Further arguments to be passed to \code{\link{PersonAlytic}} for
#' analysis. All options in \code{PersonAlytic} can be passed except for
#' \code{output}, \code{data}, \code{ids}, \code{dvs}, \code{time}, and
#' \code{phase} as these are taken from the \code{design}.
#'
#' @examples
#'
#' \dontrun{
#'
#' example(polyICT)
#' myPolyICT$inputMat
#'
#' # parametric bootstrap examples
#'
#' # lazy cloning - this can lead to errors as `myPolyICT` is updated with each
#' # call to `myPolyICT$update`
#'
#' #testICTpower10 <- ICTpower(c('testICTpower10', 'csv'),
#' #  myPolyICT, B=3,
#' #  seed = 54)
#' #testICTpower20 <- ICTpower(c('testICTpower20', 'csv'),
#' #  myPolyICT$update(groups=c(group1=20, group2=20)), B=3,
#' #  seed = 23)
#' #testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#' #  myPolyICT$update(groups=c(group1=20, group2=20),
#' #  phases=makePhase(c(20,60,20))),
#' #  B=3, seed = 23)
#'
#' # safe cloning
#' myPolyICT2 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT2$update(groups=c(group1=20, group2=20))
#' testICTpower20 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT2, B=3, seed = 23)
#'
#' myPolyICT3 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT3$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT3, B=3, seed = 23)
#'
#' myPolyICT4 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT4$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT4, B=3, seed = 23)
#'
#'
#' # non-parametric bootstrap examples
#'
#' # create a population with 500 participants per group
#' myPolyICTnonPar <- myPolyICT$clone(deep=TRUE)
#' myPolyICTnonPar$update(groups=c(group1=500, group2=500))
#' Data <- myPolyICTnonPar$makeData()
#' save(Data, file = "Data.RData")
#'
#' # non parametric bootstrap samlpes of 25 participants each group
#' ICTpower(outFile    = c("npbsTest", "csv"),
#'          B          = 100                 ,
#'          dataFile   = "Data.RData"        ,
#'          sampleSizes = c(25,25)           )
#'
#' # clean up
#' txts <- dir(getwd(), glob2rx("*.txt"))
#' csvs <- dir(getwd(), glob2rx("*.csv"))
#' file.remove("Data.RData", txts, csvs)
#'  }


ICTpower <- function(outFile         = NULL                      ,
                     design          = NULL                      ,
                     B               = 100                       ,
                     dataFile        = NULL                      ,
                     sampleSizes     = NULL                      ,
                     alpha           = .05                       ,
                     seed            = 123                       ,
                     cores           = parallel::detectCores()-1 ,
                     savePowerReport = TRUE                      ,
                     standardize     = list(dv    = TRUE ,
                                            ivs   = FALSE,
                                            byids = TRUE )      ,
                     ...
)
{
  #
  if(is.null(design) & is.null(dataFile))
  {
    stop('Please provide `design` or `dataFile`, both are `NULL`.')
  }
  if(!is.null(design) & !is.null(dataFile))
  {
    stop('Please provide only on of `design` or `dataFile` but not both.')
  }

  #argList <-  as.list(match.call(expand.dots = TRUE)[-1]) # only explicitly passed
  argList <- mget(names(formals()),sys.frame(sys.nframe()))
  #print(argList) # this is cluttering up the console

  # check file name
  if(!is.null(outFile)) outFile <- .checkFile(outFile)

  # generate seeds
  seeds <- .makeSeeds(seed, B)

  # parametric bootstrap
  if(is.null(dataFile))
  {
    #
    # parralelization
    #

    # message and timing
    message("\nStarting simulation of B=", B, " data sets.\n",
            ifelse(!is.null(outFile$sfile),
                   paste("Data will be saved in the file:\n\n", outFile$sfile,
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
      dat <- design$makeData(seeds[b])

      # rename for merging
      names(dat)[2] <- paste('y', b, sep='')

      # return
      return(dat)
    }
    cat("\n\n") # to separate progress bar from warnings and messages
    # stop the cluster
    suppressWarnings( parallel::stopCluster(cl) )

    #
    # data clean up and save
    #
    # merge data in Data
    mergeby <- names(Data[[1]])
    mergeby <- mergeby[mergeby != 'y1']
    Data    <- Reduce(function(df1, df2) merge(df1, df2, by=mergeby, all = TRUE),
                      Data)

    # if requested, save the data
    if(!is.null(outFile$sfile))
    {
      if(outFile$isRData)
      {
        save(Data, file=outFile$sfile)
      }
      if(outFile$iscsv)
      {
        write.csv(Data, file=outFile$sfile, row.names = FALSE)
      }
    }

    # print timing message
    message("Data simulation took: ",
            capture.output(Sys.time() - start), ".\n\n")

    #
    # process inputs in `...` that may be passed to `PersonAlytic`
    #
    # get the ARMA order from `correlation` if it is passed by ...
    if(!exists('correlation'))
    {
      ar          <- length(design$error$parms$ar[design$error$parms$ar!=0])
      ma          <- length(design$error$parms$ma[design$error$parms$ma!=0])
      correlation <- paste('corARMA(p=', ar, ',q=', ma, ')', sep='')
      detectAR    <- FALSE
    }
    if(!exists('detectTO'))   detectTO   <- FALSE
    if(!exists('time_power')) time_power <- design$randFxOrder

    #
    # analyze the data using PersonAlytics, treating y1,...,yB
    # as separate outcomes
    #
    phase <- NULL
    ivs   <- NULL
    int   <- NULL
    if( length(design$phases)>1) phase <- 'phase'
    if( length(design$groups)>1) ivs   <- 'group'
    if( !is.null(ivs) )          int   <- list(c(ivs, phase), c(ivs, 'Time'))


  }

  # non-parametric bootstrap
  if(!is.null(dataFile))
  {
    ext <- tools::file_ext(dataFile)
    if(ext %in% c('csv', 'CSV', 'Csv'))
    {
      Data <- read.csv(dataFile)
    }
    if(ext %in% c('RData', 'rdata', 'Rdata', 'RDATA'))
    {
      load(dataFile)
      if( !exists("Data") )
      {
        stop('The file you provided for `dataFile`:\n',
          dataFile, '\nDoes not contain a data.frame named `Data`. Please ',
          'rename you data.frame to `Data`.')
      }
    }
    else
    {
      stop('`dataFile` has the extension ', ext, ' which is not supported.\n',
           'See the documentation for `dataFile` in ?ICTpower.')
    }

    # data checks
    nameCheck <- function(Data, name)
    {
      if(! name %in% names(Data) )
      {
        stop('The variable`', name, '`` is not in the data set.')
      }
    }
    nameCheck(Data, 'id')
    nameCheck(Data, 'y')
    nameCheck(Data, 'Time')
    # phase and group are optional

    # group check
    if( "group" %in% names(Data) )
    {
      if( length(sampleSizes) != length(unique(Data$group)) )
      {
        stop('The variable `group` in the data has ', length(unique(Data$group)),
             ' groups, but `sampleSizes` has length ', length(sampleSizes), '.')
      }
    }
    if( ! "group" %in% names(Data) )
    {
      Data$group <- 'group1'
    }

    uid    <- table(Data$id, Data$group)
    nchars <- max(nchar(as.character(Data$id)))
    datL   <- list()
    for(b in 1:B)
    {
      datG <- list()
      for(g in seq_along(unique(Data$group)) )
      {
        if(b==1)
        {
          uidg <- length(uid[,g] != 0)
          if( uidg <= sampleSizes[g] )
          {
            stop('There are ', uidg, ' unique ids in ', colnames(uid)[g],
                 'which is less than `sampleSizes` which is ', sampleSizes[g],
                 '.\nSelect `sampleSizes` that is smaller than the number of',
                 'unique ids.')
          }
        }

        set.seed(seeds[b] + g)
        wid <- sample(uid[uid[,g]!=0,g], size = sampleSizes, replace = FALSE)
        dat <- Data[Data$id %in% as.numeric(names(wid)),]
        names(dat)[which(names(dat)=='y')] <- paste('y', b, sep='')

        # create new ids that will be the same across samples for merging, noting
        # that the original ids will be last
        dat$id    <- as.numeric(factor(dat$id, labels=1:sampleSizes[g])) +
                       (10*g)^nchars
        dat       <- dat[order(dat$id, dat$Time),]
        datG[[g]] <- dat
      }
      datL[[b]] <- do.call(rbind, datG)
    }
    mergeby <- names(datL[[1]])
    mergeby <- mergeby[mergeby != 'y1']
    Data    <- Reduce(function(df1, df2) merge(df1, df2, by=mergeby, all = TRUE),
                      datL)

    #
    # process inputs in `...` that may be passed to `PersonAlytic`
    #
    # get the ARMA order from `correlation` if it is passed by ...
    if(!exists('correlation'))
    {
      correlation <- 'corARMA(p=1, q=1)'
    }
    if(!exists('detectAR'))   detectAR   <- FALSE
    if(!exists('detectTO'))   detectTO   <- FALSE
    if(!exists('time_power')) time_power <- 1

    #
    # analyze the data using PersonAlytics, treating y1,...,yB
    # as separate outcomes
    #
    phase <- NULL
    ivs   <- NULL
    int   <- NULL
    if( length(unique(Data$phases))>1 ) phase <- 'phase'
    if( length(unique(Data$groups))>1 ) ivs   <- 'group'
    if( !is.null(ivs) ) int   <- list(c(ivs, phase), c(ivs, 'Time'))
  }


  paout <- PersonAlytic(output       = outFile$file                     ,
                        data         = Data                             ,
                        ids          = 'id'                             ,
                        dvs          = as.list(paste('y', 1:B, sep='')) ,
                        time         = 'Time'                           ,
                        phase        = phase                            ,
                        ivs          = ivs                              ,
                        interactions = int                              ,
                        time_power   = time_power                       ,
                        correlation  = correlation                      ,
                        detectAR     = detectAR                         ,
                        detectTO     = detectTO                         ,
                        cores        = cores                            ,
                        standardize  = standardize                      ,
                        ...
  )

  #
  # power analysis specific summary of results
  #
  powerL <- powerReport(paout, alpha, file=outFile$outFile,
                        saveReport=savePowerReport)


  #
  # distributions of the estimates
  #
  plotDists <- FALSE
  if(plotDists)
  {
    samplingDist(paout)
  }

  # return
  invisible(powerL)


}


#' powerReport - print power results to screen and to a file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
powerReport <- function(paout, alpha, file, saveReport=TRUE)
{
  # power estimates are the proportion of p < alpha
  whichP <- names(paout)[ grepl('p.value', names(paout)) ]
  powerL <- list()
  for(i in whichP)
  {
    powerL[[i]] <- mean(paout[[i]] <= alpha, na.rm = TRUE)
  }

  # effect size estimates (value)
  whichV <- names(paout)[ grepl('Value', names(paout)) ]
  whichV <- whichV[! grepl('statValue', whichV)]
  valueLm <- valueLsd <- list()
  for(i in whichV)
  {
    valueLm[[i]]  <- mean(paout[[i]], na.rm = TRUE)
    valueLsd[[i]] <- sd(paout[[i]], na.rm = TRUE)
  }

  # print the report to the screen
  names(powerL) <- gsub('.p.value', '', names(powerL))
  names(powerL) <- gsub("\\s", " ", format(names(powerL),
                                           width=max(nchar(names(powerL)))) )
  powerOutput <- paste(names(powerL)                                , '\t',
                       sprintf("% 2.2f", round(unlist(valueLm),2))  , '\t\t',
                       sprintf("% 2.2f", round(unlist(valueLsd),2)) , '\t\t',
                       sprintf("% 2.2f", round(unlist(powerL),2))   , '\t\t',
                       '\n' )
  powerOutput <- c(
    paste(gsub("\\s", " ", format("Predictor",
                                  width=max(nchar(names(powerL))))),
      '\tMean Estimates\tSD Estimates\tPower\n'), .hl(),
    powerOutput
  )

  message(.hl(), powerOutput, .hl() )

  # save the report
  if(saveReport)
  {
    powerfile <- paste(file, 'PowerReport.txt', sep='_')
    cat( powerOutput, file = powerfile)
  }

  # return results
  return( data.frame(meanEst = unlist(valueLm),
                     sdEst   = unlist(valueLsd),
                     power   = unlist(powerL ) ) )
}





