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
#' \code{outFile=c('n10_medSlopeEffect', 'csv')}
#' or \code{outFile=c('condition2', 'Rdata')}. Only `csv` and `RData` are supported.
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
#'   myPolyICT2, B=3, seed = 25, prompt=FALSE)
#'
#' myPolyICT3 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT3$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT3, B=3, seed = 26, prompt = FALSE)
#'
#' myPolyICT4 <- myPolyICT$clone(deep=TRUE)
#' myPolyICT4$update(groups=c(group1=20, group2=20),
#'   phases=makePhase(c(20,60,20)))
#' testICTpower20t100 <- ICTpower(c('testICTpower20', 'csv'),
#'   myPolyICT4, B=3, seed = 27, prompt = FALSE)
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
#' ICTpower(outFile     = c("npbsTest", "csv") ,
#'          B           = 3                    ,
#'          dataFile    = "Data.RData"         ,
#'          sampleSizes = c(25,25)             ,
#'          prompt      = FALSE                )
#'
#' # with a finite power correction passing `fpc` by ...
#' ICTpower(outFile     = c("npbsFPCtest", "csv") ,
#'          B           = 3                       ,
#'          dataFile    = "Data.RData"            ,
#'          sampleSizes = c(25,25)                ,
#'          fpc         = length(table(Data$id))  ,
#'          prompt      = FALSE                   )
#'
#' # piecewise growth model example
#' ICTpower(outFile     = c("piecewise", "csv"),
#'          B           = 3                    ,
#'          dataFile    = "Data.RData"         ,
#'          sampleSizes = c(25,25)             ,
#'          alignPhase  = 'piecewise'          ,
#'          prompt      = FALSE                )
#'
#' # clean up
#' file.remove( 'piecewise_PowerReport.txt' )
#' file.remove( 'piecewise_PersonAlytic.csv' )
#' file.remove( 'REMLlme.txt' )
#' file.remove( 'piecewise.Data.csv' )
#' file.remove( 'npbsFPCtest_PowerReportFPC.txt' )
#' file.remove( 'npbsFPCtest_PowerReport.txt' )
#' file.remove( 'npbsFPCtest_PersonAlytic.csv' )
#' file.remove( 'npbsFPCtest.Data.csv' )
#' file.remove( 'npbsTest_PowerReport.txt' )
#' file.remove( 'npbsTest_PersonAlytic.csv' )
#' file.remove( 'npbsTest.Data.csv' )
#' file.remove( 'Data.RData' )
#' file.remove( 'testICTpower20_PowerReport.txt' )
#' file.remove( 'testICTpower20_PersonAlytic.csv' )
#' file.remove( 'testICTpower20.Data.csv' )
#' file.remove( 'piecewise_mse.csv' )
#' file.remove( 'npbsFPCtest_mse.csv' )
#' file.remove( 'npbsTest_mse.csv' )
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
                     standardize     = list(dv    = FALSE ,
                                            ivs   = FALSE ,
                                            byids = FALSE )      ,
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
  if(!is.null(outFile)) outFile <- PersonAlyticsPower:::.checkFile(outFile, prompt)

  # generate seeds
  seeds <- PersonAlyticsPower:::.makeSeeds(seed, B)

  # process dots
  dots <-  list(...)
  dotsNames <- names(dots)
  if(any(dotsNames %in% 'fpc')) fpc <- dots$fpc

  # parametric bootstrap ####
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
    }
    if(!exists('time_power')) time_power <- design$randFxOrder
    if(!exists('autoDetect')) autoDetect <- list()

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

  # non-parametric bootstrap ####
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
                          " be labeled y1, ..., y", B, ".\n\n",
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

        # instead of sampling by extraction, we'll sample by deletion, i.e.,
        # if a case is not in the sample, their data are deleted. When we
        # later merge the sampled data with the raw data, the first data set
        # y0 will be the population
        dat <- Data[Data$group==paste("group", g, sep=''),]
        dat$y[!dat$id %in% wid] <- NA
        names(dat)[which(names(dat)=='y')] <- paste('y', b, sep='')

        datG[[g]] <- dat
      }
      datL[[b]] <- do.call(rbind, datG)
    }
    mergeby <- names(datL[[1]])
    mergeby <- mergeby[mergeby != 'y1']
    DataB   <- Reduce(function(df1, df2) merge(df1, df2, by=mergeby, all = TRUE),
                      datL)
    names(Data)[which(names(Data)=='y')] <- 'y0'
    Data    <- merge(Data, DataB); rm(DataB)

    # merging messes up data order, resort
    Data <- Data[order(Data$id, Data$group, Data$phase, Data$Time),]


    # qc
    #any(duplicated(DataB$id[DataB$group=='group1'], DataB$id[DataB$group=='group1']))


    #
    # process inputs in `...` that may be passed to `PersonAlytic`
    #
    # get the ARMA order from `correlation` if it is passed by ...
    if(!exists('correlation'))
    {
      correlation <- NULL
    }
    if(!exists('time_power')) time_power <- 1
    if(!exists('autoDetect')) autoDetect <- list()

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

  # distribution check
  family = gamlss.dist::NO()
  if(!is.null(design$yCut))
  {
    .l <- length(design$yCut)
    if( .l==2 ) family = gamlss.dist::BI()
    if( .l>=3 ) family = eval(parse(text=paste("MULTIN(type = '", .l, "')", sep='')))
    rm(.l)
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
    print(autoDetect)
    print(cores)
    print(standardize)
  }

  if( is.null(dataFile)) dvs   <- as.list(paste('y', 1:B, sep=''))
  if(!is.null(dataFile)) dvs   <- as.list(paste('y', 0:B, sep=''))
  paout <- PersonAlytic(output       = outFile$file ,
                        data         = Data         ,
                        ids          = 'id'         ,
                        dvs          = dvs          ,
                        time         = 'Time'       ,
                        phase        = phase        ,
                        ivs          = ivs          ,
                        interactions = int          ,
                        time_power   = time_power   ,
                        correlation  = correlation  ,
                        autoDetect   = autoDetect   ,
                        cores        = cores        ,
                        standardize  = standardize  ,
                        family       = family       ,
                        ...
  )

  # check for failed PersonAlytic calls
  success <- TRUE
  if(nrow(paout)==0)
  {
    message('\nThe call to `PersonAlytic` returned output with 0 rows.',
         '\nTry a power analysis with different settings.')
    success <- FALSE
  }
  if(all(paout$converge=="Model did not converge"))
  {
    message('\nNone of the models converged. Try a power analysis with different settings.')
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

    # population - sample comparisons for non-parametric bootstrap
    if( !is.null(dataFile) )
    {
      popSampComp(paout, file=outFile$file)
    }

    # return
    if(!exists("fpc")) invisible( powerL )
    if( exists("fpc")) invisible( list(powerL=powerL, powerLFPC=powerLFPC) )
  }


}


#' popSampComp - compare non-parametric bootstrap estimates to population
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @param paout data.frame produced by \code{\link{PersonAlytic}}.
#' @param file character. A file name with no extentsion.
popSampComp <- function(paout, file)
{
  # columns for comparison
  nms    <- names(paout)
  wstats <- nms[which(grepl("statValue", nms))]
  wests  <- nms[which(grepl("..Value", nms))[1]:length(nms)]
  wests  <- wests[!grepl("statName", wests)]
  wclmns <- c(wstats, wests)
  # exclude DF
  wclmns <- wclmns[!grepl(".DF", wclmns)]
  wclmns <- wclmns[!duplicated(wclmns)]

  # which row has y0 and which rows have ! y0
  yr <- which(paout$dv=="y0")
  yhatr <- which(paout$dv!="y0")

  # simple stats (mirroring power report)
  stats <- data.frame(
    pop_value = unlist(paout[yr, wclmns])                        ,
    samp_mean = apply(paout[yhatr, wclmns], 2, mean, na.rm=TRUE) ,
    samp_sd   = apply(paout[yhatr, wclmns], 2, sd  , na.rm=TRUE)
  )

  # column names
  mnms <- c("mse" , "mse.norm.lo" , "mse.norm.hi" , "mse.bca.lo" , "mse.bca.hi" ,
           "rmse", "rmse.norm.lo", "rmse.norm.hi", "rmse.bca.lo", "rmse.bca.hi",
           "mae" , "mae.norm.lo" , "mae.norm.hi" , "mae.bca.lo" , "mae.bca.hi"
           )

  # loop
  mses <- list()
  for(i in seq_along(wclmns))
  {
    w <- wclmns[i]
    y <- paout[[w]][yr]
    yhat <- paout[[w]][yhatr]
    x  <- data.frame(yhat=yhat, y=y)
    b  <- list()
    ci.mse <- ci.rmse <- ci.mae <- list()
    try(b  <- boot::boot(x, mseb, 2000), silent = TRUE)
    try({
    ci.mse  <- boot::boot.ci(b, type=c("norm", "bca"), index=1)
    ci.rmse <- boot::boot.ci(b, type=c("norm", "bca"), index=2)
    ci.mae  <- boot::boot.ci(b, type=c("norm", "bca"), index=3)}, silent = TRUE)
    all     <- c(b$t0[1], ci.mse$normal[2:3] , ci.mse$bca[4:5]  ,
                 b$t0[2], ci.rmse$normal[2:3], ci.rmse$bca[4:5] ,
                 b$t0[3], ci.mae$normal[2:3] , ci.mae$bca[4:5]  )
    if(length(mnms) == length(all)) names(all) <- mnms
    mses[[w]] <- data.frame(t(all))
  }
  mses <- data.frame(stat=wclmns, stats, plyr::rbind.fill(mses))
  row.names(mses) <- NULL

  # statNames
  mses$stat <- as.character(mses$stat)
  # this line may not generalize to multiple columns, but shouldn't have to (20190812)
  statNames <- as.character( unlist( paout[yr, nms[grepl("statName", nms)]] ) )
  mses$stat[1:length(statNames)] <- statNames
  mses$stat <- unlist(mses$stat)

  # save
  file <- paste(file, "_mse.csv", sep='')
  write.csv(mses, file, row.names=FALSE)
}

#' mse - mean squared error & mean absolute error
#' @export
#' @param x data.frame with columns \code{yhat} (the estimates) and \code{y}
#' (the population value which should generally have zero variance).
mse <- function(x)
{
  n <- nrow(x)
  mse <- (1/n) * sum((x$y-x$yhat)^2)
  mae <- (1/n) * sum(abs(x$y-x$yhat))
  c(mse=mse, rmse=sqrt(mse), mae=mae)
}

#' mseb - version for bootstraping
#' @keywords internal
#' @import boot
#' @param x See \code{\link{mse}}
#' @param i The index for resampling
mseb <- function(x,i)
{
  mse(x[i,])
}


#' powerReport - print power results to screen and to a file
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @param paout The *.cvs output from PersonAlytics read in to the object paout.
#'
#' @param alpha The type I error rate.
#'
#' @param file A file name for for saving the results if \code{saveReport} is
#' \code{TRUE}.
#'
#' @param saveReport Logical. Should the report be saved (in addition to
#' being printed to the console).
#'
#' @param fpc Logical. Should finite population corrected p-values be looked for
#' in \code{paout}.
#'
#' @return A data frame with the mean of the parameter estimates, the standard
#' deviation of the parameter estimates, and the power.
#'
powerReport <- function(paout, alpha, file, saveReport=TRUE, fpc=FALSE,
                        printToScreen=TRUE)
{
  # exclude y0
  paout <- paout[paout$dv!='y0',]

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

  # select whichV based on FPC
  whichFPC <- grepl('FPC', whichV)
  if(any(whichFPC))
  {
    if( fpc) whichV <- whichV[ whichFPC]
    if(!fpc) whichV <- whichV[!whichFPC]
  }

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

  if(printToScreen)
  {
    if(!fpc) message("Power Report")
    if( fpc) message("Power Report with Finite Population Correction")
    message(.hl(), powerOutput, .hl() )
  }

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





