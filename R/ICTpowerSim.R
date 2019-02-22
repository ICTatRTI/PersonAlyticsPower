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
  fnames <- names(designs)

  set.seed(seed)
  seeds <- ceiling(runif(2, 0, 9e6) )

  powerL <- list()
  for(i in seq_along(designs))
  {
    powerL[[i]] <-
    ICTpower(file         = c(fnames[i]  , save) ,
             design       = designs[[i]]         ,
             B            = B                    ,
             checkDesign  = 'no'                 ,
             alpha        = alpha                ,
             randFxFamily = randFxFamily         ,
             randFxParms  = randFxParms          ,
             randFxSeed   = seeds[1]             ,
             errorParms   = errorParms           ,
             errorFUN     = errorFUN             ,
             errorFamily  = errorFamily          ,
             errorSeed    = seeds[2]             ,
             cores        = cores                ,
             ...
             )
  }
  powerL <- do.call(rbind, powerL)
  row.names(powerL) <- fnames
  reportName <- paste(pReportName, 'PAP', packageVersion('PersonAlyticsPower'),
                      'PA', packageVersion('PersonAlytics'), '.csv', sep='_')
  write.csv(powerL, file=reportName)
}
