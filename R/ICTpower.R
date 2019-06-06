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
#' @param prompt Logical. Default is \code{TRUE}. If TRUE, a check for
#' whether \code{outFile} already exists will be done and if so, the user
#' will be prompted whether they want to overwirte \code{outFile}. The
#' \code{\link{ICTpowerSim}} function turns this off so the user does not
#' have to monitor a series of simulations.
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
#' #  B=3, seed = 24)
#'
#' # safe cloning
#' myPolyICT2 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT2$update(groups=c(group1=20, group2=20))
#' testICTpower20 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT2, B=3, seed = 25)
#'
#' myPolyICT3 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT3$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT3, B=3, seed = 26)
#'
#' myPolyICT4 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT4$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT4, B=3, seed = 27)
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
#'          B          = 3                   ,
#'          dataFile   = "Data.RData"        ,
#'          sampleSizes = c(25,25)           )
#'
#' # with a finite power correction passing `fpc` by ...
#' ICTpower(outFile    = c("npbsFPCtest", "csv") ,
#'          B          = 3                       ,
#'          dataFile   = "Data.RData"            ,
#'          sampleSizes = c(25,25)               ,
#'          fpc        = length(table(Data$id))  )
#'
#' # piecewise growth model example
#' ICTpower(outFile     = c("piecewise", "csv"),
#'          B           = 3                    ,
#'          dataFile    = "Data.RData"         ,
#'          sampleSizes = c(25,25)             ,
#'          alignPhase  = 'piecewise'          )
#'
#' # clean up
#' txts <- dir(getwd(), glob2rx("*.txt"))
#' csvs <- dir(getwd(), glob2rx("*.csv"))
#' file.remove("Data.RData", txts, csvs)
#'  }
#'

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
                                            byids = TRUE )       ,
                     prompt          = TRUE                      ,
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

  # check file name
  if(!is.null(outFile)) outFile <- .checkFile(outFile, prompt)

  # generate seeds
  seeds <- .makeSeeds(seed, B)

  # process dots
  dots <-  list(...)
  dotsNames <- names(dots)
  if(any(dotsNames %in% 'fpc')) fpc <- dots$fpc

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
    if(!file.exists(dataFile))
    {
      stop('The file `', dataFile, '` does not exist. The non-parametric\n',
           'bootsrap estimates of power cannot be estitmated.')
    }
    ext <- tools::file_ext(dataFile)
    if(tolower(ext) == tolower('CSV'))
    {
      Data <- read.csv(dataFile)
    }
    if(tolower(ext) == tolower('RDATA'))
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
      stop('`dataFile` has the extension `', ext, '` which is not supported.\n',
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
    groups <- unique(Data$group)

    # start message
    message("\nStarting bootstrap resampling of n=", sum(sampleSizes),
            " participants repeated B=", B, " times.\n",
            ifelse(!is.null(outFile$sfile),
                    paste("Data will be saved in the file:\n\n", outFile$sfile,
                          "\n\nwhere the outcome for each bootstrap sample will",
                          " be labeled y1, ..., y", B,
                          ".\n\n NOTE: Individual ids are not saved and are",
                          "\narbitrarily relabeled to create a data structure",
                          "\nthat conforms to PersonAlytic() requirements.\n\n",
                          sep=""),
                    '\n\n')
    )

    nchars <- max(nchar(as.character(Data$id)))
    datL   <- list()
    for(b in 1:B)
    {
      datG <- list()
      for(g in seq_along(groups) )
      {
        if(b==1)
        {
          uidg <- length(unique(Data$id[Data$group==groups[g]]))
          if( uidg <= sampleSizes[g] )
          {
            stop('There are ', uidg, ' unique ids in ', colnames(uid)[g],
                 'which is less than `sampleSizes` which is ', sampleSizes[g],
                 '.\nSelect `sampleSizes` that is smaller than the number of',
                 'unique ids.')
          }
        }

        # When sampling from a finite population, we'd never sample the same
        # person twice. To mimic this, we sample without replacement.

        set.seed(seeds[b] + g)
        sid <- unique(Data$id[Data$group==groups[g]])
        wid <- sample(sid, size = sampleSizes[g], replace = FALSE)
        dat <- Data[Data$id %in% wid,]
        names(dat)[which(names(dat)=='y')] <- paste('y', b, sep='')

        # before sorting, nullify row names as they can affect sorting
        row.names(dat) <- NULL
        # create new ids that will be the same across samples for merging, noting
        # that the original ids will be last
        dat$id    <- as.numeric(factor(dat$id, labels=1:sampleSizes[g])) +
                       (10*g)^nchars
        # now sort on the new names
        dat       <- dat[order(dat$id, dat$Time),]
        datG[[g]] <- dat
      }
      datL[[b]] <- do.call(rbind, datG)
    }
    mergeby <- names(datL[[1]])
    mergeby <- mergeby[mergeby != 'y1']
    Data    <- Reduce(function(df1, df2) merge(df1, df2, by=mergeby, all = TRUE),
                      datL)

    # qc
    #any(duplicated(Data$id[Data$group=='group1'], Data$id[Data$group=='group1']))


    #
    # process inputs in `...` that may be passed to `PersonAlytic`
    #
    # get the ARMA order from `correlation` if it is passed by ...
    if(!exists('correlation'))
    {
      correlation <- NULL
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

    wp <- which(tolower(names(Data)) %in% c('phase', 'phases'))
    wg <- which(tolower(names(Data)) %in% c('group', 'groups'))

    if( length(unique((Data[[wp]]))) > 1 ) phase <- names(Data)[wp]
    if( length(unique((Data[[wg]]))) > 1 ) ivs   <- names(Data)[wg]

    if( !is.null(ivs) ) int   <- list(c(ivs, phase), c(ivs, 'Time'))
  }

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

  if(exists("debugforeach"))
  {
    for(i in seq_along(dotsNames))
    {
      print( paste(dotsNames[i], dots[i], sep=': '))
    }
    print(outFile$file)
    print(head(Data))
    print(B)
    print(phase)
    print(ivs)
    print(int)
    print(time_power)
    print(correlation)
    print(detectAR)
    print(detectTO)
    print(cores)
    print(standardize)
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

  # check for failed PersonAlytic calls
  success <- TRUE
  if(nrow(paout)==0)
  {
    message('\nThe call to `PersonAlytic` returned output with 0 rows.',
         '\nPower analysis cannot be completed.')
    success <- FALSE
  }
  if(all(paout$converge=="Model did not converge"))
  {
    message('\nNone of the models converged. Power analysis cannot be completed.')
    success <- FALSE
  }
  if(!success)
  {
    message("\nPower analysis failed. See the `PersonAlytic` output in\n`",
            paste(outFile$file, "csv", sep="."), "` for error messages.\n\n")
  }

  if(success)
  {
    #
    # power analysis specific summary of results
    #
    powerL <- powerReport(paout, alpha, file=outFile$file,
                          saveReport=savePowerReport)

    if(exists("fpc"))
    {
      powerLFPC <- powerReport(paout, alpha, file=outFile$file,
                               saveReport=savePowerReport, fpc=TRUE)
    }

    #
    # distributions of the estimates
    #
    plotDists <- FALSE
    if(plotDists)
    {
      samplingDist(paout)
    }

    # return
    if(!exists("fpc")) invisible( powerL )
    if( exists("fpc")) invisible( list(powerL=powerL, powerLFPC=powerLFPC) )
  }


}


#' powerReport - print power results to screen and to a file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
powerReport <- function(paout, alpha, file, saveReport=TRUE, fpc=FALSE)
{
  # power estimates are the proportion of p < alpha
  whichP <- names(paout)[ grepl('p.value', names(paout)) ]

  # select whichP based on FPC
  whichFPC <- grepl('FPC', whichP)
  if(any(whichFPC))
  {
    if( fpc) whichP <- whichP[ whichFPC]
    if(!fpc) whichP <- whichP[!whichFPC]
  }

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
  powerOutput <- paste(sprintf("%s"    , names(powerL))            , '\t',
                       sprintf("% 2.3f", round(unlist(valueLm),3))  , '\t',
                       sprintf("% 2.3f", round(unlist(valueLsd),3)) , '\t',
                       sprintf("% 2.3f", round(unlist(powerL),3))   , '\t',
                       '\n' )
  powerOutput <- c(
    paste(gsub("\\s", " ", format("Predictor",
                                  width=max(nchar(names(powerL))))),
      '\tMean Estimates\tSD Estimates\tPower\n'), .hl(),
    powerOutput
  )

  if(!fpc) message("Power Report")
  if( fpc) message("Power Report with Finite Population Correction")
  message(.hl(), powerOutput, .hl() )

  # save the report
  if(saveReport)
  {
    fileExt <- 'PowerReport.txt'
    if(fpc) fileExt <- 'PowerReportFPC.txt'
    powerfile <- paste(file, fileExt, sep='_')

    if( file.exists(powerfile) )
    {
      message("\nThe file `", powerfile, "` already exists and will be overwritten.")
    }

    dt  <- format(Sys.time(), format='%Y%m%d_%H.%M%p')
    pap <- paste("PersonAlyticsPower Version", packageVersion("PersonAlyticsPower"))
    pa  <- paste("PersonAlytics Version", packageVersion('PersonAlytics'))

    cat( dt, "\n", pap, "\n", pa, "\n\n", powerOutput, file = powerfile)
  }

  # return results
  return( data.frame(meanEst = unlist(valueLm),
                     sdEst   = unlist(valueLsd),
                     power   = unlist(powerL ) ) )
}





