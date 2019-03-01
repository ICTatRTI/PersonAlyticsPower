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
#' The remaining parameters are as explained in \code{\link{ICTpower}}.
#'
#' @examples
#'
#' designs <- list()
#' designs[["defaults"]] <- polyICT$new()
#' designs[["defaultsN20"]] <- designs[["defaults"]]$clone(deep=TRUE)
#' designs[["defaultsN20"]]$n <- 20
#' ICTpowerSim(designs, B=3)
#'

#TODO: move items explicitly passed by ICTsimSetup() to top
ICTpowerSim <- function(designs                                  ,
                        pReportName  = "ICTpowerSimResults"      ,
                        B            = 1000                      ,
                        seed         = 1                         ,
                        save         = NULL                      ,
                        alpha        = .05                       ,
                        randFxFamily = "qNO"                     ,
                        randFxParms  = list(mu=.5 , sigma=1)     ,
                        errorParms   = list(ar=c(.5), ma=c(0))   ,
                        errorFUN     = arima.sim                 ,
                        errorFamily  = NULL                      ,
                        cores        = parallel::detectCores()-1 ,
                        ...
                        )
{
  # check whether `pReportName` is already used
  pReports <- dir(getwd(), glob2rx(paste("*", pReportName, "*")))
  if( length(pReports) > 0 )
  {
    stop("The `pReportName` ", pReportName, "is already used in this directory.",
         "\nUse a different `pReportName`.")
  }

  # create a temporory directory
  upDir   <- getwd()
  tempDir <- tempfile(pattern='', getwd())
  dir.create(tempDir)
  setwd(tempDir)

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
            " (", round(100*(i/length(designs)),2), "% of ",
            length(designs), " designs)",
            msgDelim)
    powerL[[i]] <-
    ICTpower(file            = c(fnames[i]  , save) ,
             design          = designs[[i]]         ,
             B               = B                    ,
             checkDesign     = 'no'                 ,
             alpha           = alpha                ,
             randFxFamily    = randFxFamily         ,
             randFxParms     = randFxParms          ,
             randFxSeed      = seeds[i, 1]          ,
             errorParms      = errorParms           ,
             errorFUN        = errorFUN             ,
             errorFamily     = errorFamily          ,
             errorSeed       = seeds[i, 2]          ,
             cores           = cores                ,
             savePowerReport = FALSE                , # if there is a crash, this is bad
             ...
             )
  }

  #TODO: if savePowerReport = TRUE, the glob2rx must include txt files

  # package all the results into a tar file
  csvs       <- dir(getwd(), glob2rx("*.csv"), full.names = TRUE)
  csvsCtime  <- unlist( lapply(as.list(csvs), function(x) file.info(x)$ctime) )
  sinceStart <- csvsCtime > start
  csvs       <- csvs[sinceStart]
  tarName    <- paste(upDir, '/', pReportName, '.tar.gz', sep='')
  setwd(upDir)
  tempDir <- unlist( strsplit(tempDir, '\\\\') )[2]
  tar(tarName, tempDir, 'g')
  #TODO: this only empties the directory
  unlink(dir(getwd(), glob2rx(paste("*",tempDir,"*",sep="")), full.names = TRUE),
         TRUE, TRUE)

  # save the results
  powerL <- do.call(rbind, powerL)
  row.names(powerL) <- fnames
  reportName <- paste(pReportName, 'PAP', packageVersion('PersonAlyticsPower'),
                      'PA', packageVersion('PersonAlytics'), '.csv', sep='_')
  write.csv(powerL, file=reportName)

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
    #TODO: make this a unique extension that the user can name, eg
    #mySimStudy.dsgn
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
#' \code{\link{ployICT}}.
#'
#' @param phasesL A list of phases. See \code{n} in
#' \code{\link{ployICT}}.
#'
#' @param propErrVarL Numeric List. See \code{n} in
#' \code{\link{ployICT}}.
#'
#' @param effectSizesL A list of effect sizes. See \code{n} in
#' \code{\link{ployICT}}.
#'
#'  @param corMatL A list of numeric correlation matrices. See \code{n} in
#' \code{\link{ployICT}}.
#'
#' @param randFxVarL A list of numeric vectors. See \code{n} in
#' \code{\link{ployICT}}.
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
#' effectSizesL <- list( list(randFx=list(intercept=0, slope=.25),
#'                       fixdFx=list(phase=.5, phaseTime=.25)),
#'                       list(randFx=list(intercept=0, slope=.5),
#'                       fixdFx=list(phase=.5, phaseTime=.25)))
#' corMatL <- list(matrix(c(1,.2,.2,1), 2, 2), matrix(c(1,.6,.6,1), 2, 2))
#' randFxVarL <- list(c(1, .1), c(1, .2))
#' B <- 3
#' nFolders <- 4
#' ICTpowerSimOptions <-"cores = 4"
#' ICTsimSetup(seed    ,
#'   nL                ,
#'   phasesL           ,
#'   propErrVarL       ,
#'   effectSizesL      ,
#'   corMatL           ,
#'   randFxVarL        ,
#'   B                 ,
#'   nFolders          ,
#'   ICTpowerSimOptions)
#' }
ICTsimSetup <- function(seed                      ,
                        nL                        ,
                        phasesL                   ,
                        propErrVarL               ,
                        effectSizesL              ,
                        corMatL                   ,
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
  effectSizes   <- seq_along(effectSizesL)
  corMats       <- seq_along(corMatL)
  randFxVars    <- seq_along(randFxVarL)

  # set up the designs
  designI <- expand.grid(nL           ,
                         nObservations,
                         nBaseline    ,
                         propErrVarL  ,
                         effectSizes  ,
                         corMats      ,
                         randFxVars   )

  names(designI) <- c('n'             ,
                      'nObservations' ,
                      'nBaseline'     ,
                      'propErrVar'    ,
                      'effectSize'    ,
                      'corMat'        ,
                      'randFxVar'     )


  designI$ConditionName <- paste(
    'n',      designI$n             ,
    'o',      designI$nObservations ,
    'nbl',    designI$nBaseline     ,
    'pev',    designI$propsErrVar   ,
    'fx',     designI$effectSize    ,
    'r',      designI$corMat        ,
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
  set.seed(seed)
  designI$seeds <- ceiling(runif(nrow(designI), 0, 9e6) )

  # save objects for calling by the scripts produced in the loop
  save(designI                           ,
       nL                                ,
       phasesL                           ,
       propErrVarL                       ,
       effectSizesL                      ,
       corMatL                           ,
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
        "    effectSizes = effectSizesL[[designI$effectSize[i]]],\n",
        "    corMat = corMatL[[designI$corMat[i]]],\n",
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




