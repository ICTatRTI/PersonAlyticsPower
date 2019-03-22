#' active - active bindings for ICTdesign R6 class definitions and its children
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

    phases = function(value)
    {
      if( missing(value) ){ private$.phases }
      else
      {
        nObservations <- length(c(unlist((value))))
        if( nObservations < 10 )
        {
          stop("`desgin` has ", nObservations, " total timepoints.\n",
               "At least 10 are required")
        }
        private$.phases <- value
        self
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

    effectSizes = function(value)
    {
      if( missing(value) ){ private$.effectSizes }
      else
      {
        if( !is.list(value) )
        {
          stop('`effectSizes` should be a named list, not the provided value\n',
               value)
        }
        if("polyICT" %in% class(self))
        {
          checkEffectSizesPoly(value)
        }
        private$.effectSizes <- value
        self
      }
    },

    effectSizes2 = function(value)
    {
      if( missing(value) ){ private$.effectSizes }
      else
      {
        if( !is.list(value) )
        {
          stop('`effectSizes2` should be a named list, not the provided value\n',
               value)
        }
        if("polyICT" %in% class(self))
        {
          checkEffectSizesPoly(value)
        }
        private$.effectSizes <- value
        self
      }
    },

    corMat = function(value)
    {
      if( missing(value) ){ private$.corMat}
      else
      {
        if("polyICT" %in% class(self))
        {
          checkCorMat(value)
        }
        private$.corMat <- value
        self
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

    muFUN = function(value)
    {
      if( missing(value) ){ private$.muFUN }
      else
      {
        if( !is.function(value) ) stop('`muFUN` must be a function.')
        private$.muFUN <- value
      }

    },

    SigmaFun = function(value)
    {
      if( missing(value) ){ private$.SigmaFun }
      else
      {
        if( !is.function(value) ) stop('`SigmaFun` must be a function.')
        private$.SigmaFun <- value
      }

    },

    polyOrder = function(value)
    {
      if( missing(value) ){ private$.polyOrder }
      else
      {
        stop('`polyOrder` is constructed from `effectSizes` during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
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

    covMat = function(value)
    {
      if( missing(value) ){ private$.covMat }
      else
      {
        stop('`covMat` is constructed from other inputs during the',
             ' creation of a\n`', class(self)[1], '` object and cannot',
             ' be updated.')
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
    },

    unStdEffects = function(value)
    {
      if( missing(value) ){ private$.unStdEffects }
      else
      {
        stop('`unStdEffects` is constructed from other inputs during the',
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
             private = list(
              .n             = NULL,
              .phases        = NULL,
              .propErrVar    = NULL
             ),

             # Note: put all current and future active bindings in 'active()'
             # the child classes of ICTdesign will inherit them
             active = active(),

             public = list(
                initialize = function
                (
                  n          = NULL,
                  phases     = NULL,
                  propErrVar = NULL
                )
             {
               # input checks here

               # populate private
               private$.n             <- n
               private$.phases        <- phases
               private$.propErrVar    <- propErrVar
             },

             print = function(...)
             {
               phases <- paste(lapply(self$phases, function(x) shQuote(x[1])),
                               'for',
                               lapply(self$phases, length), 'time points\n')
               phases <- paste(c('', rep(' ', length(phases)-1)), phases)
               cat(
                 paste("An ICT with", self$n, "participants and",
                       length(c(unlist(self$phases))), "time points.\nPhases:\n",
                       paste(phases, collapse = ''),
                       "Error Variance is ", 100*self$propErrVar,
                       "% of the total variance.")
               )
               invisible(self)
             }


             )
) # end of ICTdesign class definition




# use inheritance to create specific instances of ICTdesign




