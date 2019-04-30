#' active - active bindings for ICTdesign R6 class definitions and ALL children
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @keywords internal
#'
active <- function()
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
        if( !isNum | !isInt )
        {
          stop("`n`, the number of participants, ",
               "must be a positive integer\ninstead of the ",
               "provided value `", value, "`")
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
        if( value < .01 | value > .99 )
        {
          stop("`propErrVar`, the proportion of the total variance ",
               "that is error variance,\nmust be > .01 and < .99")
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
        private$.randFxVar <- value
      }

    },

    randFxCorMat = function(value)
    {
      if( missing(value) ){ private$.randFxCorMat}
      else
      {
        if("polyICT" %in% class(self))
        {
          checkCorMat(value)
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


#' \code{ICTdesign} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @examples
#' # create a basic ABA design
#' basicICT <- ICTdesign$new(n=10, phases=makePhase(), propErrVar=.5)
#' basicICT
#' basicICT$n <- 20
#' basicICT$n <- 'A'
#' basicICT$propErrVar <- 2

ICTdesign <- R6::R6Class("ICTdesign",

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

             # Child classes of ICTdesign will inherit active bindings and
             # all active bindings for children will stay in active()
             active = active(),

             public = list(

               # initialize is child-specific, so ICTdesign does not have an
               # initialize function

               print = function(...)
               {
                 phases <- paste(lapply(self$phases, function(x) shQuote(x[1])),
                                 'for',
                                 lapply(self$phases, length), 'time points\n')
                 phases  <- paste(c('', rep(' ', length(phases)-1)), phases)
                 nGroups <- length(self$groupNames)
                 if(nGroups==0) nGroups <- 1
                 cat(
                   "An ICT with ", sum(self$n), " participants, ",
                   length(c(unlist(self$phases))), " time points, and ",
                   ifelse(!is.null(self$groupNames), nGroups, 1),
                   ifelse(nGroups==1, " group.", " groups.\n"),
                   "\nPhases:\n",
                   paste(phases, collapse = ''),
                   "\nGroups:\n", paste(self$groupNames, 'n=', self$n, '\n'),
                   "\nThe variance at the first time point is partitioned as\n",
                   100*self$propErrVar[1], "% random effects variance,\n",
                   100*self$propErrVar[2], "% residual autocorrelation variance,\n",
                   100*self$propErrVar[3], "% measurement error variance.",
                   sep=""
                 )
                 invisible(self)
               }

             )
) # end of ICTdesign class definition




# use inheritance to create specific instances of ICTdesign
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





