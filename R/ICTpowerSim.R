#' ICTpowerSim
#'
#' @export
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param designs A named list of designs.
#' in \code{\link{ICTpower}}. The names in the named list will be passed to
#' the parameter \code{file} in \code{\link{ICTpower}}.
#'
#' @param pReportName Character. A name used to create output files. The default
#' is \code{"ICTpowerSimResults"}
#'
#' @param B Numeric (integer). The number of replicated data sets per \code{design}
#' in \code{designs}.
#'
#' @param Seed Numeric (integer). The random seed for ensuring the simulation
#' can be replicated.
#'
#' @param save \code{NULL} or character. Should the data be saved? If NULL
#' (the default) raw data will not be saved. Other
#' options include \code{save='csv'} and \code{save='RData'}. This will be
#' combined with the names of \code{designs} to pass to the \code{file} option
#' in \code{\link{ICTpower}}.
#'
#' @param alpha See \code{\link{ICTpower}}.
#'
#' @param cores See \code{\link{ICTpower}}.
#'
#' @param standardize See \code{\link{ICTpower}}.
#'
#' @param dotar Logical. The default is \code{FALSE}. Should the data files
#' (if requested, see \code{save}) and analysis results be saved in a *.tar
#' file with the same name as \code{pReportName} and then delete the unzipped
#' copies of these files? This is suggested wheth the number of designs is
#' large. The summary of power across all conditions will no be included in the
#' *.tar archive. It is strongly encourage that a test run with only three
#' designs and B=3 is attempted at first to make sure archiving is done
#' correctly. If the path names resulting from your working directory are
#' too long, archiving will fail with the error "storing paths of more than 100
#' bytes is not portable". See \code{\link{tar}}.
#'
#' @examples
#'
#' \dontrun{
#'
#' # set up two designs
#' designs <- list()
#' designs[["defaultsN20"]] <- polyICT$new()
#' designs[["defaultsN40"]] <- designs[["defaultsN20"]]$clone(deep=TRUE)
#' designs[["defaultsN40"]]$update(groups = c(group1=40, group2=40))
#'
#' # run the simulation (without saving output, console print only)
#' ICTpowerSim(designs, B=3)
#'
#' }

ICTpowerSim <- function(designs                                  ,
                        pReportName  = "ICTpowerSimResults"      ,
                        B            = 1000                      ,
                        seed         = 1                         ,
                        save         = NULL                      ,
                        alpha        = .05                       ,
                        cores        = parallel::detectCores()-1 ,
                        standardize  = list(dv    = TRUE ,
                                            ivs   = FALSE,
                                            byids = TRUE )        ,
                        dotar        = FALSE                     ,
                        ...
                        )
{
  # check whether `pReportName` is already used
  if( file.exists(pReportName) )
  {
    stop("The `pReportName` ", pReportName, "is already used in this directory.",
         "\nUse a different `pReportName`.")
  }

  # create a analysis directory
  upDir   <- getwd()
  dir.create(pReportName)
  setwd(pReportName)

  start  <- Sys.time()

  fnames <- names(designs)

  set.seed(seed)
  nDesigns <- length(designs)
  seeds    <- matrix(ceiling(runif(2*nDesigns, 0, 9e6)), nDesigns, 2)

  msgDelim <- paste('\n\n', paste(rep('*', 80), collapse = ''), '\n\n')

  powerL <- list()
  for(i in seq_along(designs))
  {
    message(msgDelim, ' Starting design: ', names(designs)[i],
            " (", i, " of ",
            length(designs), " designs)",
            msgDelim)
    powerL[[i]] <-
    ICTpower(outFile         = c(fnames[i]  , save) ,
             design          = designs[[i]]         ,
             B               = B                    ,
             alpha           = alpha                ,
             seed            = seeds[i]             ,
             cores           = cores                ,
             savePowerReport = TRUE                 ,
             ...
             )
  }

  # delete text files
  txts <- dir(getwd(), glob2rx("*.txt"))
  file.remove(txts)

  if(dotar)
  {
    # find csvs created since the start of the run
    csvs       <- dir(getwd(), glob2rx("*.csv"), full.names = TRUE)
    csvsCtime  <- unlist( lapply(as.list(csvs), function(x) file.info(x)$ctime) )
    sinceStart <- csvsCtime > start
    csvs       <- csvs[sinceStart]

    # find Rdatas created since the start of the run
    Rdatas      <- dir(getwd(), glob2rx("*.Rdata"), full.names = TRUE)
    RdatasCtime <- unlist( lapply(as.list(Rdatas), function(x) file.info(x)$ctime) )
    sinceStart  <- RdatasCtime > start
    Rdatas      <- Rdatas[sinceStart]

    # package all the results into a tar file
    tarName    <- paste(pReportName, '.tar.gz', sep='')
    tar(tarName, c(csvs, Rdatas), 'g')

    # remove the other files if requested
    file.remove(txts, csvs, Rdatas)
  }

  # save the results
  power <- data.frame(do.call(rbind, lapply(powerL, function(x) x$power))   )
  mEst  <- data.frame(do.call(rbind, lapply(powerL, function(x) x$meanEst)) )
  sdEst <- data.frame(do.call(rbind, lapply(powerL, function(x) x$sdEst))   )
  row.names(power) <- fnames
  row.names(mEst)  <- fnames
  row.names(sdEst) <- fnames
  power$type <- "power"
  mEst $type <- "meanEst"
  sdEst$type <- "sdEst"
  reportName  <- paste(pReportName, 'PAP', packageVersion('PersonAlyticsPower'),
                      'PA', packageVersion('PersonAlytics'), '.csv', sep='_')
  powerOut <- rbind(power, mEst, sdEst)
  names(powerOut) <- c(row.names(powerL[[1]]), 'type')
  write.csv(powerOut, file=reportName)

  # move back to the parent directory
  setwd(upDir)
}


#' simStatus
#'
#' @export
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param studyName
#'
#' @param studyDirectory Character. The full path to the directory the
#' simulation study is being run in
#'
simStatus <- function(studyName, studyDirectory=getwd())
{
  # save and set wd
  curDir <- getwd()
  setwd(studyDirectory)

  # get the designI file
  if(! file.exists('ICT_Sim_Conditions.Rdata') )
  {
    stop("The file `ICT_Sim_Conditions.Rdata` does not exist in the directory\n",
         studyDirectory)
  }
  designI <- read.csv('ICT_Sim_Conditions.Rdata')
  nConditions <- nrow(designI)

  # get the directories
  simDirs <- dir(studyDirectory, glob2rx(paste("*", studyName, "*", sep='')))

  studyChecker <- function(x)
  {
    # csvs that are output
    csvs <- dir(x, glob2rx("*_PersonAlytics_*.csv"), recursive = TRUE)

    # clvs that are data
    csvd <- dir(x, glob2rx("*.csv"), recursive = TRUE)
    csvd <- csvd[! csvd %in% csvs]

    # results
    c(outputCount=length(csvs), dataCount=length(csvd))
  }

  csvsCounts <-
    data.frame(simDirs    = simDirs,
               do.call(rbind, lapply(as.list(simDirs), studyChecker)) )

  print(csvsCounts)
  message("\n\nThe design file `ICT_Sim_Conditions.Rdata` shows there are ", nConditions,
          " conditions in the study.\n",
          sum(csvsCounts$outputCount), "/",
          nConditions, " output *.csv files are completed, (",
          round(100*(sum(csvsCounts$outputCount)/nConditions),2), "%)",
          " and\n", sum(csvsCounts$dataCount), "/", nConditions,
          " data files have been saved (",
          round(100*(sum(csvsCounts$dataCount)/nConditions),2), "%).")

  # reset the wd
  setwd(curDir)
}



# this evolved from
# R:\PaCCT\02 Clients\NIH\0216822.000.001 R21 ICT\03 Simulation Studies\Study1\simulationStudy1.R
#' ICTsimSetup - function to automate setting up \code{designs} to be passed to
#' \code{\link{ICTpowerSim}}, including splitting up the designs over sevar folders
#' so the simulation can be run in chunks or on multiple computers.
#'
#' @export
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @param seed Numeric. A random seed for this simulation.
#'
#' @param nL Numeric list. A list of sample sizes. See \code{n} in
#' \code{\link{polyICT}}.
#'
#' @param phasesL A list of phases. See \code{n} in
#' \code{\link{polyICT}}.
#'
#' @param propErrVarL Numeric List. See \code{n} in
#' \code{\link{polyICT}}.
#'
#' @param randFxMeanL A list of effect sizes. See \code{n} in
#' \code{\link{polyICT}}.
#'
#'  @param randFxCorMatL A list of numeric correlation matrices. See \code{n} in
#' \code{\link{polyICT}}.
#'
#' @param randFxVarL A list of numeric vectors. See \code{n} in
#' \code{\link{polyICT}}.
#'
#' @param nFolders Numeric. The number of folders to divide the simulation into.
#' The default is \code{NULL}, in which case the user is prompted for input at
#' the console after the total number of conditions is computed. To prevent
#' folders that are too large, it is reccomended that \code{nFolders} not
#' exceed 100, especially on Windows systems.
#'
#' @param ICTpowerSimOptions Character. Other options to be passed to the
#' function \code{\link{ICTpowerSim}}. For example,
#'
#' @param muFUN
#'
#' @examples
#' \dontrun{
#' seed <- 1234
#' nL <- list(10, 20)
#' phasesL <- list( makePhase(c(10,10), c(0,1)), makePhase(c(5,15), c(0,1)) )
#' propErrVarL <- list(.5, .75)
#' randFxMeanL <- list( list(randFx=list(intercept=0, slope=.25),
#'                       fixdFx=list(phase=.5, phaseTime=.25)),
#'                       list(randFx=list(intercept=0, slope=.5),
#'                       fixdFx=list(phase=.5, phaseTime=.25)))
#' randFxCorMatL <- list(matrix(c(1,.2,.2,1), 2, 2), matrix(c(1,.6,.6,1), 2, 2))
#' randFxVarL <- list(c(1, .1), c(1, .2))
#' B <- 3
#' nFolders <- 4
#' ICTpowerSimOptions <-"cores = 4"
#' ICTsimSetup(seed    ,
#'   nL                ,
#'   phasesL           ,
#'   propErrVarL       ,
#'   randFxMeanL      ,
#'   randFxCorMatL           ,
#'   randFxVarL        ,
#'   B                 ,
#'   nFolders          ,
#'   ICTpowerSimOptions)
#' }

ICTsimSetup <- function(seed                      ,
                        nL                        ,
                        phasesL                   ,
                        propErrVarL               ,
                        randFxMeanL              ,
                        randFxCorMatL                   ,
                        randFxVarL                ,
                        B                  = 1000 ,
                        nFolders           = NULL ,
                        ICTpowerSimOptions = ""   )
{
  #TODO: add some validation for ICTpowerSimOptions

  # if the Rdata file already exists, error
  if( file.exists('ICT_Sim_Conditions.Rdata') )
  {
    stop("The file `ICT_Sim_Conditions.Rdata` already exists in the current",
         " working directory,\n\n",
         getwd(), "\n\nPlease move or delete `ICT_Sim_Conditions.Rdata`.",
         "\nThis file creates all possible combinations\nof the conditions input",
         " to the function `ICTsimSetup`.")
  }

  # check inputs
  if( any(lapply(phasesL, length) > 2) )
  {
    #TODO: only works for a 2-phase design
    stop("`ICTsimSetup` currently only works for 2-phase designs.")
  }

  # restructure the lists for expand.grid
  nL            <- unlist( nL )
  nObservations <- unlist( lapply(phasesL, function(x) length(unlist(x))) )
  nBaseline     <- unlist( lapply(phasesL, function(x) length(x[[1]])) )
  propErrVarL   <- unlist( propErrVarL )
  randFxMean   <- seq_along(randFxMeanL)
  randFxCorMats       <- seq_along(randFxCorMatL)
  randFxVars    <- seq_along(randFxVarL)

  # set up the designs
  designI <- expand.grid(nL           ,
                         nObservations,
                         nBaseline    ,
                         propErrVarL  ,
                         randFxMean  ,
                         randFxCorMats      ,
                         randFxVars   )

  names(designI) <- c('n'             ,
                      'nObservations' ,
                      'nBaseline'     ,
                      'propErrVar'    ,
                      'effectSize'    ,
                      'randFxCorMat'        ,
                      'randFxVar'     )


  designI$ConditionName <- paste(
    'n',      designI$n             ,
    'o',      designI$nObservations ,
    'nbl',    designI$nBaseline     ,
    'pev',    designI$propsErrVar   ,
    'fx',     designI$effectSize    ,
    'r',      designI$randFxCorMat        ,
    'V',      designI$randFxVar     ,
    sep='')

  # send a message about the simulation design size
  #TODO: estimate run time
  message("\n\nThe current simulation study will have ", nrow(designI),
          " conditions.\n\n")

  # get nFolders
  if(is.null(nFolders))
  {
    prompt   <- "Enter the number of folders to break the simulation into: "
    nFolders <- as.integer( readline(prompt = prompt) )

    if(nrow(designI)/nFolders < 10)
    {
      stop('The choice of `nFolders` will result in fewer than 10 conditions',
           '\nper folder. Use a smaller value for `nFolders`.')
    }
  }

  # break up designI into nFolders subsets and set up directories
  designI$dir <- sort( rep(1:nFolders, length.out = nrow(designI)) )

  # set seeds
  designI$seeds <- .makeSeeds(seed, nrow(designI))

  # save objects for calling by the scripts produced in the loop
  save(designI                           ,
       nL                                ,
       phasesL                           ,
       propErrVarL                       ,
       randFxMeanL                      ,
       randFxCorMatL                           ,
       randFxVarL                        ,
       nFolders                          ,
       ICTpowerSimOptions                ,
       file = 'ICT_Sim_Conditions.Rdata' )

  for(i in unique(designI$dir))
  {
    d <- paste('simStudy1dir', i, sep='')
    dir.create( d )

    cat(file = paste(d, '/simStudy1dir', i, '.r', sep=''),
        paste("setwd('", paste(getwd(), d, sep="/"), "')", "\n", sep=''),
        "load('../ICT_Sim_Conditions.Rdata')\n",
        "designs <- list()\n",
        "library(PersonAlyticsPower)\n",
        paste("for(i in which(designI$dir==",i,"))\n"),
        "{\n",
        "  designs[[designI$ConditionName[i]]] <- polyICT$new(\n",
        "    n = designI$n[i],\n",
        "    phases = makePhase(c(designI$nBaseline[i], \n",
        "                         designI$nObservations[i]-designI$nBaseline[i]),\n",
        "                       c(0,1)),\n",
        "    propErrVar = designI$propErrVar[i],\n",
        "    randFxMean = randFxMeanL[[designI$effectSize[i]]],\n",
        "    randFxCorMat = randFxCorMatL[[designI$randFxCorMat[i]]],\n",
        "    randFxVar = randFxVarL[[designI$randFxVar[i]]]\n",
        "  )\n",
        "}\n",
        paste("ICTpowerSim(designs,\n",
              "pReportName = paste('dir', ",i,", sep=''),\n",
              "B=", B, ",\n",
              "seed=", designI$seeds[i], ",\n",
              ICTpowerSimOptions,
              ")\n")
    )
  }
}




