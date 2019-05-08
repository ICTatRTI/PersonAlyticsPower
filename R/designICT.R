#' active - active bindings for designICT R6 class definitions and its children
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
.active <- function()
{
  list(

    # fields that can be updated by the user
    inputMat = function(value)
    {
      if( missing(value) ){ private$.inputMat }
      else
      {
        private$.inputMat <- value
      }
    },

    randFxVar = function(value)
    {
      if( missing(value) ){ private$.randFxVar }
      else
      {
        private$.randFxVar <- value
      }
    },

    randFxCorMat = function(value)
    {
      if( missing(value) ){ private$.randFxCorMat }
      else
      {
        private$.randFxCorMat <- value
      }
    },

    randFxCovMat = function(value)
    {
      if( missing(value) ){ private$.randFxCovMat }
      else
      {
        private$.randFxCovMat <- value
      }
    },

    propErrVar = function(value)
    {
      if( missing(value) ){ private$.propErrVar }
      else
      {
        private$.propErrVar <- value
      }
    },

    error = function(value)
    {
      if( missing(value) ){ private$.error }
      else
      {
        private$.error <- value
      }
    },

    merror = function(value)
    {
      if( missing(value) ){ private$.merror }
      else
      {
        private$.merror <- value
      }
    },

    yMean = function(value)
    {
      if( missing(value) ){ private$.yMean }
      else
      {
        private$.yMean <- value
      }
    },

    ySD  = function(value)
    {
      if( missing(value) ){ private$.ySD }
      else
      {
        private$.ySD <- value
      }
    },


    # fields than cannot be updated by the user
    n = function(value)
    {
      if( missing(value) ){ private$.n }
      else
      {
        private$.n <- value
      }
    },

    nObs = function(value)
    {
      if( missing(value) ){ private$.nObs }
      else
      {
        private$.nObs <- value
      }
    },

    groups = function(value)
    {
      if( missing(value) ){ private$.groups }
      else
      {
        private$.groups <- value
      }
    },

    phases= function(value)
    {
      if( missing(value) ){ private$.phases }
      else
      {
        private$.phases <- value
      }
    },

    designMat = function(value)
    {
      if( missing(value) ){ private$.designMat }
      else
      {
        private$.designMat <- value
      }
    },

    unStdRandFxMean = function(value)
    {
      if( missing(value) ){ private$.unStdRandFxMean }
      else
      {
        private$.unStdRandFxMean <- value
      }
    },

    phaseNames = function(value)
    {
      if( missing(value) ){ private$.phaseNames }
      else
      {
        private$.phaseNames <- value
      }
    },

    groupNames = function(value)
    {
      if( missing(value) ){ private$.groupNames }
      else
      {
        private$.groupNames <- value
      }
    },

    randFxOrder = function(value)
    {
      if( missing(value) ){ private$.randFxOrder }
      else
      {
        private$.randFxOrder <- value
      }
    },

    meanNames = function(value)
    {
      if( missing(value) ){ private$.meanNames }
      else
      {
        private$.meanNames <- value
      }
    },

    varNames = function(value)
    {
      if( missing(value) ){ private$.varNames }
      else
      {
        private$.varNames <- value
      }
    },


    # fields not currently implementd
    randFxFam = function(value)
    {
      if( missing(value) ){ private$.randFxFam }
      else
      {
        private$.randFxFam <- value
      }
    },

    randFxFamParms = function(value)
    {
      if( missing(value) ){ private$.randFxFamParms }
      else
      {
        private$.randFxFamParms <- value
      }
    }
  )
}


#' \code{designICT} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @examples
#' # create a basic ABA design
#' basicICT <- designICT$new(n=10, phases=makePhase(), propErrVar=.5)
#' basicICT
#' basicICT$n <- 20
#' basicICT$n <- 'A'
#' basicICT$propErrVar <- 2

designICT <- R6::R6Class("designICT",

                         # this list must contain all child slots, children cannot update
                         # private (to my knowledge)
                         private = list(

                           # editable without a new call to $new
                           .inputMat          = NULL,
                           .randFxVar         = NULL,
                           .randFxCorMat      = NULL,
                           .randFxCovMat      = NULL,
                           .propErrVar        = NULL,
                           .error             = NULL,
                           .merror            = NULL,
                           .yMean             = NULL,
                           .ySD               = NULL,

                           # not editable
                           .n                 = NULL,
                           .nObs              = NULL,
                           .groups            = NULL,
                           .phases            = NULL,
                           .designMat         = NULL,
                           .unStdRandFxMean   = NULL,
                           .phaseNames        = NULL,
                           .groupNames        = NULL,
                           .randFxOrder       = NULL,
                           .meanNames         = NULL,
                           .varNames          = NULL,

                           # not implemented
                           .randFxFam         = NULL,
                           .randFxFamParms    = NULL
                         ),

                         # Child classes of designICT will inherit active bindings and
                         # all active bindings for children will stay in active()
                         active = .active(),

                         public = list(

                           # initialize is child-specific, so designICT does not have an
                           # initialize function

                           print = function(...)
                           {
                             phases <- paste(lapply(self$phases, function(x) shQuote(x[1])),
                                             'for',
                                             lapply(self$phases, length), 'time points\n')
                             phases  <- paste(c('', rep(' ', length(phases)-1)), phases)
                             nGroups <- length(self$groupNames)
                             if(nGroups==0) nGroups <- 1
                             cat(.hl(),
                                 "An ICT with ", sum(self$n), " participants, ",
                                 length(c(unlist(self$phases))), " time points, and ",
                                 ifelse(!is.null(self$groupNames), nGroups, 1),
                                 ifelse(nGroups==1, " group.", " groups.\n"), .hl(),
                                 "\nPhases:\n",
                                 paste(phases, collapse = ''),
                                 "\nGroups:\n", paste(self$groupNames, 'n=', self$n, '\n'),
                                 .hl(),
                                 sep=""
                             )
                             cat("Input Matrix:\n")
                             cat(.hl())
                             print( self$inputMat )
                             cat(.hl())
                             cat("The input matrix can be edited by typing `edit(x$inputmat)`\n",
                                 "where `x` is the name of your ICT design object.\n")
                             invisible(self)
                           },

                           designCheck = function(file=NULL, ylim=NULL, fitMod=FALSE,
                                                  seed=123, npg=5000)
                           {

                             # save and reset n, due to inheritance design will get overwritten, fix
                             # n below
                             originaln <- self$groups
                             tempn <- rep(npg, length(originaln))
                             names(tempn) <- names(self$groups)
                             self$groups <- tempn; rm(tempn)

                             # simulate data
                             dat <- self$makeData(seed=seed)

                             # needs to be reimplemented
                             if(1==2)
                             {
                               # get expected means and variances
                               expectedVar <- self$expectedVariances

                               # compare expected to observed variances
                               expObsVar   <- cbind( expectedVar$Variances,
                                                     aggregate(dat$y, by=list(dat$Time), var)$x)
                               # get the correlation of expected and observed variances
                               expObsVcor <- round(cor(expObsVar)[1,2],4)
                               cat("The correlation between the expected variances and the observed\n",
                                   "variances are ", expObsVcor, "\n\n")

                               # compare expected to observed means
                               expObsMean <- cbind( expectedVar$Means,
                                                    aggregate(dat$y, by=list(dat$Time), mean)$x)
                               expObsMcor <- round(cor(expObsMean)[1,2],4)
                               cat("The correlation between the expected means and the observed\n",
                                   "means are ", expObsMcor, "\n\n")

                               # compare expected to observed total
                               TotalMean <- mean(dat$y)
                               TotalVar  <- var(dat$y)
                               cat("The observed total mean is ", TotalMean,
                                   "\nThe expected total mean is ", expectedVar$TotalMean,
                                   "\nThe observed total variance is ", TotalVar,
                                   "\nThe expected total variance is ", expectedVar$TotalVar, "\n\n")

                               #TODO this only works for linear models
                               #TODO needs better matching, force name matching in polyICT
                               Estimates <- round(rbind(summary(mod0)$tTable[,1]),3 )
                               Inputs    <- round(c(unlist(self$unStdRandFxMean))[c(1,3,2,4)],3)
                               cat('\n\nCheck the effect size estimates against inputs:\n')
                               print( data.frame(Inputs=Inputs, Estimates=t(Estimates)) )

                               cat("\n\n\n")
                             }

                             # get the data and, if requested, fit the model
                             correlation <- paste("corARMA(p=", length(self$error$parms$ar), ", ",
                                                  "q=", length(self$error$parms$ma), ")", sep="")
                             pa   <- Palytic$new(data=dat, ids='id', dv='y', time='Time',
                                                 phase=unlist(ifelse(length(self$phaseNames)>1,
                                                                     'phase',list(NULL))),
                                                 ivs=unlist(ifelse(length(self$groupNames)>1,
                                                                   'group',list(NULL))),
                                                 time_power = self$maxRandFx,
                                                 correlation = correlation)
                             if(fitMod) # runs slow with some examples, qc why
                             {
                               mod0 <- pa$lme()
                               print( summary( mod0 ) )
                             }

                             # save data if requested
                             if(!is.null(file)) save(mod0, file=paste(file[1],
                                                                      'designCheck.RData',
                                                                      sep='_'))


                             # plot
                             if( length( self$groupNames ) == 1 ) print( pa$plot(ylim=ylim) )
                             if( length( self$groupNames ) >= 2 ) print( pa$plot(groupvar = 'group',
                                                                                 ylim=ylim) )

                             # restore the original sample sizes
                             self$groups <- originaln

                             # if the model was fit, return it
                             # TODO consider the conseuences of returning pa intsead of self
                             if(fitMod)  invisible(pa)
                           }



                         )
) # end of designICT class definition




# use inheritance to create specific instances of designICT
#
# current:
# - polyICT: polynomial ICT
# - polyICT2: polynomial ICT via distribution subsetting
#
# Planned (see Singer & Willett page 234, Table 6.7):
# - polyICT
# - hyperICT
# - invPolyICT
# - expICT
# - negExpICT
# - logisticICT





