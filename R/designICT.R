#' active - active bindings for designICT R6 class definitions and it heirs
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
.active <- function()
{
  list(

    n = function(value)
    {
      if( missing(value) ){ private$.n }
      else
      {
        isNum <- is.numeric(value)
        isInt <- TRUE
        if(isNum) isInt <- round(value)==value
        if( !all(isNum) | !all(isInt) )
        {
          stop("`n`, the number of participants, ",
               "must be a positive integer\ninstead of the ",
               "provided value `", value, "`")
        }
        if(length(value) != length(private$.n))
        {
          value <- rep(value, length(private$.n))
          names(value) <- names(private$.n)
        }
        private$.n <- value
        self
      }
    },

    nObservations = function(value)
    {
      if( missing(value) ){ private$.nObservations }
      else
      {
        stop('`nObservations` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

    phases = function(value)
    {
      if( missing(value) ){ private$.phases }
      else
      {
        if( private$.nObservations < 10 )
        {
          warning("`desgin` has ", nObservations, " total time points.\n",
               "At least 10 time points are reccomended.")
        }
        private$.phases <- value
        self
      }
    },

    phaseNames = function(value)
    {
      if( missing(value) ){ private$.phaseNames }
      else
      {
        if( length(private$.phases) != length(value) )
        {
          stop("The `desgin` has ", length(private$.phases), " phases.\n",
               "Only ", length(value), " `phaseNames` were provided:\n",
               paste(value, sep='\n'))
        }
        private$.phaseNames <- value
        self
      }
    },

    groupNames = function(value)
    {
      if( missing(value) ){ private$.groupNames }
      else
      {
        if( length(private$.groupNames) != length(value) )
        {
          stop("The `desgin` has ", length(private$.groupNames), " groups.\n",
               "Only ", length(value), " `groupNames` were provided:\n",
               paste(value, sep='\n'))
        }
        private$.groupNames <- value
        self
      }
    },

    designMat = function(value)
    {
      if( missing(value) ){ private$.designMat }
      else
      {
        stop('`designMat` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

    propErrVar = function(value)
    {
      if( missing(value) ){ private$.propErrVar }
      else
      {
        if( any(value < .01) | any(value > .99) )
        {
          stop("The three values in `propErrVar` must be > .01 and < .99.\n",
               "You provided (", paste(value, collapse=', '), ").")
        }
        if( sum(value) != 1 )
        {
          stop("The three values in `propErrVar` should sum to one.\n",
               "The sum of (", paste(value, collapse=', '), ") is ", sum(value), ".")
        }
        private$.propErrVar <- value
        self
      }
    },

    randFxMean = function(value)
    {
      if( missing(value) ){ private$.randFxMean }
      else
      {
        if( !is.list(value) )
        {
          stop('`randFxMean` should be a named list, not the provided value\n',
               value)
        }
        if("polyICT" %in% class(self))
        {
          checkRandFxMeanPoly(value)
        }
        private$.randFxMean <- value
        self
      }
    },

    unStdRandFxMean = function(value)
    {
      if( missing(value) ){ private$.unStdRandFxMean }
      else
      {
        stop('`unStdRandFxMean` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

    randFxVar = function(value)
    {
      if( missing(value) ){ private$.randFxVar }
      else
      {
        if(( !is.list(value) & (length(value) != (private$.maxRandFx + 1)) ))
        {
          stop('This design has ', private$.maxRandFx + 1, ' random effects,\n',
               'while you provided ', length(value))
        }
        private$.randFxVar <- randFxVarPop(private$.phaseNames,
                                           private$.groupNames,
                                           value, 'n')
        self
      }
    },

    randFxCorMat = function(value)
    {
      if( missing(value) ){ private$.randFxCorMat }
      else
      {
        if("polyICT" %in% class(self))
        {
          temp <- checkPolyICT(private$.n, private$.randFxMean,
                               private$.phases, private$.randFxCorMat,
                               private$.randFVar)

        }
        private$.randFxCorMat <- value
        self
      }
    },

    randFxCovMat = function(value)
    {
      if( missing(value) ){ private$.randFxCovMat }
      else
      {
        stop('`randFxCovMat` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

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
    },

    maxRandFx = function(value)
    {
      if( missing(value) ){ private$.maxRandFx }
      else
      {
        stop('`maxRandFx` is constructed from `randFxMean` during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

    error = function(value)
    {
      if( missing(value) ){ private$.error }
      else
      {
        if( !class(value) %in% "err" ) stop('`error` must be an `err` class object.')
        private$.error <- value
      }

    },

    merror = function(value)
    {
      if( missing(value) ){ private$.merror }
      else
      {
        if( !class(value) %in% "err" ) stop('`merror` must be an `err` class object.')
        private$.merror <- value
      }

    },

    variances = function(value)
    {
      if( missing(value) ){ private$.variances }
      else
      {
        stop('`variances` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
      }
    },

    expectedVariances = function(value)
    {
      if( missing(value) ){ private$.expectedVariances }
      else
      {
        stop('`expectedVariances` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
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

              # design
              .n                 = NULL,
              .nObservations     = NULL,
              .phases            = NULL,
              .phaseNames        = NULL,
              .groupNames        = NULL,
              .designMat         = NULL,
              .propErrVar        = NULL,
              # random effects
              .randFxMean        = NULL,
              .unStdRandFxMean   = NULL,
              .randFxVar         = NULL,
              .randFxCorMat      = NULL,
              .randFxCovMat      = NULL,
              .randFxFam         = NULL,
              .randFxFamParms    = NULL,
              .maxRandFx         = NULL,
              # errors
              .error             = NULL,
              .merror            = NULL,
              # variances
              .variances         = NULL,
              .expectedVariances = NULL
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
                 message(hl(),
                   "An ICT with ", sum(self$n), " participants, ",
                   length(c(unlist(self$phases))), " time points, and ",
                   ifelse(!is.null(self$groupNames), nGroups, 1),
                   ifelse(nGroups==1, " group.", " groups.\n"), hl(),
                   "\nPhases:\n",
                   paste(phases, collapse = ''),
                   "\nGroups:\n", paste(self$groupNames, 'n=', self$n, '\n'),
                   "\nThe variance at the first time point is partitioned as\n",
                   100*self$propErrVar[1], "% random effects variance,\n",
                   100*self$propErrVar[2], "% residual autocorrelation variance,\n",
                   100*self$propErrVar[3], "% measurement error variance.\n",
                   hl(),
                   sep=""
                 )
                 invisible(self)
               },

               designCheck = function(file=NULL, ylim=NULL, fitMod=FALSE)
               {

                 # save and reset n, due to inheritance design will get overwritten, fix
                 # n below
                 originaln <- self$n
                 tempn <- rep(5000, length(originaln))
                 names(tempn) <- names(self$n)
                 self$n <- tempn; rm(tempn)

                 # simulate data
                 dat <- self$makeData()

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
                 self$n <- originaln

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





