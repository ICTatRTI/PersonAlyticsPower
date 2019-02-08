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

    design = function(value)
    {
      if( missing(value) ){ private$.design }
      else
      {
        nObservations <- length(c(unlist((value))))
        if( nObservations < 10 )
        {
          stop("`desgin` has ", nObservations, " total timepoints.\n",
               "At least 10 are required")
        }
        private$.design <- value
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
          checkEffectSizesPoly(effectSizes)
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

    }


  )
}
