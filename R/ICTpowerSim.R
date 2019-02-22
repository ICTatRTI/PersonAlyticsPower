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
#' @param save Should the data be saved? If NULL (the default) raw data will not be
#' saved. Other options include \code{save='csv'} and \code{save='RData'}. This will be
#' combined with the names of \code{designs} to pass to the \code{file} option in
#' \code{\link{ICTpower}}.
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

ICTpowerSim <- function(designs                                  ,
                        pReportName  = 'ICTpowerSimResults'      ,
                        save         = NULL                      ,
                        B            = 1000                      ,
                        seed         = 1                         ,
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
  seeds <- ceiling(runif(2, 0, 9e6) )

  msgDelim <- paste('\n\n', paste(rep('*', 80), collapse = ''), '\n\n')

  powerL <- list()
  for(i in seq_along(designs))
  {
    message(msgDelim, ' Starting design: ', names(designs)[i], msgDelim)
    powerL[[i]] <-
    ICTpower(file            = c(fnames[i]  , save) ,
             design          = designs[[i]]         ,
             B               = B                    ,
             checkDesign     = 'no'                 ,
             alpha           = alpha                ,
             randFxFamily    = randFxFamily         ,
             randFxParms     = randFxParms          ,
             randFxSeed      = seeds[1]             ,
             errorParms      = errorParms           ,
             errorFUN        = errorFUN             ,
             errorFamily     = errorFamily          ,
             errorSeed       = seeds[2]             ,
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
  unlink(tempDir, TRUE, TRUE)

  # save the results
  powerL <- do.call(rbind, powerL)
  row.names(powerL) <- fnames
  reportName <- paste(pReportName, 'PAP', packageVersion('PersonAlyticsPower'),
                      'PA', packageVersion('PersonAlytics'), '.csv', sep='_')
  write.csv(powerL, file=reportName)

}
