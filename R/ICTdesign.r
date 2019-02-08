
#' ICTdesign
#'
#' @examples
#' # create a basic ABA design
#' basicICT <- ICTdesign$new(n=10, design=makePhase(), propErrVar=.5)
#' basicICT
#' basicICT$n <- 20
#' basicICT$n <- 'A'
#' basicICT$propErrVar <- 2

ICTdesign <- R6::R6Class("ICTdesign",
             private = list(
              .n             = NULL,
              .design        = NULL,
              .propErrVar    = NULL
             ),

             # Note: put all current and future active bindings in 'active()'
             # the child classes of ICTdesign will inherit them
             active = active(),

             public = list(
                initialize = function
                (
                  n          = NULL,
                  design     = NULL,
                  propErrVar = NULL
                )
             {
               # input checks here

               # populate private
               private$.n             <- n
               private$.design        <- design
               private$.propErrVar    <- propErrVar
             },

             print = function(...)
             {
               phases <- paste(lapply(self$design, function(x) shQuote(x[1])),
                               'for',
                               lapply(self$design, length), 'time points\n')
               phases <- paste(c('', rep(' ', length(phases)-1)), phases)
               cat(
                 paste("An ICT with", self$n, "participants and",
                       length(c(unlist(self$design))), "time points.\nPhases:\n",
                       paste(phases, collapse = ''),
                       "Error Variance is ", 100*self$propErrVar,
                       "% of the total variance.")
               )
               invisible(self)
             }


             )
) # end of ICTdesign class definition




# use inheritance to create specific instances of ICTdesign

#' @field effectSizes A length 2 list with the names `randFx` and `fixdFx`. Each
#' are themselves a list with the coefficients for a polynomial growth model.
#' The length of `randFx` and `fixdFx` corresponds to the order of the model,
#' where an order of 2 would be a linear model and an order of 3 would be a
#' quadratic model. An example of a quadratic model would be
#' list(randFx = list(intercept = 0, slope = .5, quad = .2),
#'      fixdFx = list(phase = 1, phaseSlope = .2, phaseQuad = .1))
#' @field corMat a symmetric correlation matrix with a dimension equal to the
#' order of the model. For example, a quadratic model with correpspond to a 3x3
#' matrix.
#' @field randFxVar a vector of the same length as the order of the polynomial
#' model containing the variances of the random effects. For example, in a
#' quadratic model, `randFxVar` would be length three, the first element would
#' be the intercept variance, the second element would be the slope variance,
#' and the third element would be the variance of the quadratic term.

polyICT <- R6::R6Class("polyICT",

             inherit = ICTdesign,

             private = list(
               .n           = NULL,
               .design      = NULL,
               .propErrVar  = NULL,
               .effectSizes = NULL,
               .corMat      = NULL,
               .randFxVar   = NULL,
               .muFUN       = NULL,
               .SigmaFun    = NULL
             ),

             public  = list(
               initialize = function
               (
                 n           = NULL                       ,
                 design      = NULL                       ,
                 propErrVar  = NULL                       ,
                 effectSizes = NULL                       ,
                 corMat      = matrix(c(1,.2,.2,1), 2, 2) ,
                 randFxVar   = c(1, .1)                   ,
                 muFUN       = function(x) x              ,
                 SigmaFun    = cor2cov
               )
               {
                 # general input validation
                 checkPolyICT(effecSizes, corMat, randFxVar)

                 # populate private
                 private$.n           <- n
                 private$.design      <- design
                 private$.propErrVar  <- propErrVar
                 private$.effectSizes <- effectSizes
               },

               print = function(...)
               {
                 super$print()
                 cat("\n\neffectSizes:\n",
                     paste(names(self$effectSizes), '=',
                           unlist(self$effectSizes), '\n') )
               }

             )
)

test <- linearICT$new(n=10, design=makePhase(), propErrVar=.5,
                      effectSizes = list(intercept =  0,
                                         phase     = .5,
                                         time      = .5,
                                         phaseTime = .2))


