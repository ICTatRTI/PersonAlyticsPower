#' \code{polyICT} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
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
       .n                 = NULL,
       .phases            = NULL,
       .propErrVar        = NULL,
       .effectSizes       = NULL,
       .corMat            = NULL,
       .randFxVar         = NULL,
       .muFUN             = NULL,
       .SigmaFun          = NULL,
       .polyOrder         = NULL,
       .designMat         = NULL,
       .covMat            = NULL,
       .nObservations     = NULL,
       .variances         = NULL,
       .expectedVariances = NULL,
       .unStdEffects      = NULL
     ),

     public  = list(
       initialize = function
       (
         n                 = 10                                        ,
         phases            = makePhase()                               ,
         propErrVar        = .75                                       ,
         effectSizes       = list(randFx=list(intercept=0, slope=.5),
                                  fixdFx=list(phase=.5, phaseTime=.25)),
         corMat            = matrix(c(1,.2,.2,1), 2, 2)                ,
         randFxVar         = c(1, .1)                                  ,
         muFUN             = function(x) x                             ,
         SigmaFun          = cor2cov                                   ,
         polyOrder         = NULL                                      ,
         designMat         = NULL                                      ,
         covMat            = NULL                                      ,
         nObservations     = NULL                                      ,
         variances         = NULL                                      ,
         expectedVariances = NULL                                      ,
         unStdEffects      = NULL
       )
       {
         # general input validation
         polyOrder <- checkPolyICT(effectSizes, corMat, randFxVar) - 1

         # construct the fixed effects design matrix
         designMat <- getICTdesign(phases, polyOrder, 'polyICT')

         # get the covariance matrix
         covMat <- cor2cov(corMat, randFxVar)

         # rescale the randFx effect sizes (on Cohen's d scale) to raw scale units
         randFxSD <- sqrt(diag(covMat))
         unStdEffects <- effectSizes
         for(i in seq_along(effectSizes[[1]]))
         {
           unStdEffects[[1]][[i]] <- effectSizes[[1]][[i]]*randFxSD[i]
         }

         # get the number of observations
         nObservations <- length(c(unlist(phases)))

         # get the total variance and the error variance @ time = 1
         totalVar  <- sum(covMat)/(1-propErrVar)
         errorVar  <- propErrVar * totalVar
         variances <- list(totalVar = totalVar,
                           errorVar = errorVar,
                           randFxVar = totalVar - errorVar)

         # get the expected varariances
         # initial values
         expectedVariances <- expectedVar(covMat, designMat, variances,
                                          effectSizes, nObservations, n)
         unStdEffects      <- expectedVariances$effectSizes

         # populate private
         private$.n                 <- n
         private$.phases            <- phases
         private$.propErrVar        <- propErrVar
         private$.effectSizes       <- effectSizes
         private$.corMat            <- corMat
         private$.randFxVar         <- randFxVar
         private$.muFUN             <- muFUN
         private$.SigmaFun          <- SigmaFun
         private$.polyOrder         <- polyOrder
         private$.designMat         <- designMat
         private$.covMat            <- covMat
         private$.nObservations     <- nObservations
         private$.variances         <- variances
         private$.expectedVariances <- expectedVariances
         private$.unStdEffects      <- unStdEffects

       },

       print = function(...)
       {
         # use super to call the print function of the parent class
         super$print()


         # add print options for the child class
         cat("\n\nEffect sizes:\n",
             paste(unlist(lapply(self$effectSizes, names)), '=',
                   unlist(self$effectSizes), '\n'),
             "\nRandom Effects\n",
             "\nCorrelation matrix:\n", catMat(self$corMat),
             "\nCovariance matrix:\n",  catMat(round(self$covMat,3))
         )

       },

       makeData = function(randFx, errors, ymean=NULL, yvar=NULL)
       {
         # make each component into a matrix of dimension
         # n by nObservations

         # fixed and random effects
         fe <- re <- list()
         fev <- rev <- list()
         for(i in seq_along(self$unStdEffects[[1]]))
         {
           # fixed effects
           feName  <- names(self$unStdEffects$fixdFx)[i]
           fe[[i]] <- matrix((self$unStdEffects$fixdFx[[feName]] *
                              self$designMat[[feName]]),
                             self$n, self$nObservations, byrow=TRUE)

           # intercepts
           if(i==1)
           {
             re[[i]] <- matrix((self$unStdEffects$randFx[[i]] + randFx[,i]),
                               self$n, self$nObservations)
           }
           # slopes and higher
           if(i>1)
           {
             re[[i]] <- matrix(self$unStdEffects$randFx[[i]] + randFx[,i]) %*%
                         t(self$designMat$Time^(i-1))
             # QC - should be strictly linear for i==2
             #longCatEDA::longContPlot(re[[i]])
           }

           # variance QC
           fev[[i]] <- apply(fe[[i]], 2, var)
           rev[[i]] <- apply(re[[i]], 2, var)
         }

         # sum non-error components
         Ystar <- Reduce(`+`, fe) + Reduce(`+`, re)

         # sum together the effects
         Y <- Ystar + errors

         # QC
         #apply(Ystar, 2, var)
         #apply(Y, 2, var)

         # construct a long dataset
         data <- data.frame(
           id    = sort(rep(1:self$n, self$nObservations)) ,
           y     = c(t(Y))                                 ,
           do.call(rbind, replicate(self$n,
                   self$designMat, simplify = FALSE))
         )

         # QC
         all.equal(aggregate(data$y, list(data$Time), mean)$x,
                   unname(apply(Y, 2, mean)))
         all.equal(aggregate(data$y, list(data$Time), var)$x,
                   unname(apply(Y, 2, var)))

         # make phase a factor (for later analysis and plotting)
         data$phase <- factor(data$phase)

         # rescale the y variance if !is.null
         if(!is.null(yvar) | !is.null(ymean))
         {
           if(is.null(ymean)) ymean <- FALSE
           if(is.null(yvar))  yvar  <- 1
           data$y <- scale(yvar, ymean, TRUE) * sqrt(yvar)
         }

         # TODO: consider having filename option and saving here

         # return the data
         invisible(data)

       }

     )
)


# TODO:need to generalize beyond slopes
#' expectedVar
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
expectedVar <- function(covMat, designMat, variances,
                        effectSizes, nObservations, n)
{
  # 'total' N
  N <- nObservations * n

  # variance at each time point
  vyt <-  covMat[1,1]            +
    designMat$Time^2*covMat[2,2] +
    2*covMat[1,2]                +
    variances$errorVar

  # total variance due to random effects and error variance
  sigmaT <- (N-1)^-1 * ( (n - 1) * sum(vyt) )

  # rescale the fixed effect sizes in terms of the variance when time = 1
  whereIsTeq1 <- which(designMat$Time==1)
  sdYt1 <- sqrt( vyt[whereIsTeq1] )
  effectSizes$fixdFx$phase     <- effectSizes$fixdFx$phase     * sdYt1
  effectSizes$fixdFx$phaseTime <- effectSizes$fixdFx$phaseTime * sdYt1

  # mean at each time point
  mt <-  (effectSizes$randFx$intercept                       +
          effectSizes$randFx$slope     * designMat$Time      +
          effectSizes$fixdFx$phase     * designMat$phase     +
          effectSizes$fixdFx$phaseTime * designMat$phaseTime
  )

  # overall mean
  m = mean(mt)

  # total variance due to fixed effects stacked over time
  muT <- (N-1)^-1 * ( n * sum((mt-m)^2) )

  # overall variance
  vy     <- sigmaT + muT
  #  (N-1)^-1 * ( n * sum((mt-m)^2) + (n - 1) * sum(vyt) )


  return( list(effectSizes = effectSizes ,
               Variances   = vyt         ,
               Means       =  mt         ,
               TotalMean   =   m         ,
               TotalVar    =  vy         ) )
}

#' checkPolyICT
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'

checkPolyICT <- function(effectSizes, corMat, randFxVar)
{
  # check conformity of the elements in `effectSizes`
  effectSizes <<- effectSizes
  nFx <- checkEffectSizesPoly(effectSizes)

  # check that `corMat` is a correlation matrix
  checkCorMat(corMat)

  # check conformity of `corMat` with `effectSizes`
  test3 <- nrow(corMat) == nFx[[1]]
  if(!test3) stop('The number of rows and columns in `corMat` is ', nrow(corMat),
                  '\nbut it should equal the ', nFx[[1]], ' elements in\n',
                  ' `effectSizes$randFx` and `effectSizes$fixdFx`.')

  # check conformity of `randFxVar` with `effectSizes`
  test4 <- length(randFxVar) == nFx[[1]]
  if(!test4) stop('The number of elements in `randFxVar` is ', length(randFxVar),
                  '\nbut it should equal the ', nFx[[1]], ' elements in\n',
                  ' `effectSizes$randFx` and `effectSizes$fixdFx`.')

  # check that resulting covariance matrix is legit
  checkCorMat(cor2cov(corMat, randFxVar), FALSE)

  # return the polynomial order
  invisible( unname(unlist(nFx[1])) )
}

#' checkEffectSizes
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
checkEffectSizesPoly <- function(effectSizes)
{
  # check conformity of the elements in `effectSizes`
  FxNames <- c( "randFx", "fixdFx")
  test1   <- all( names(effectSizes) %in% FxNames)
  if(!test1) stop('`effectSizes` must be a length 2 list with the names ',
                  paste('`', FxNames, '`', collapse=', ', sep=''))
  nFx   <- lapply(effectSizes, length)
  test2 <- nFx[[1]]==nFx[[2]]
  if(!test2) stop('`effectSizes$randFx` must have the same number of elements\n',
                  'as `effectSizes$fixdFx`.')
  invisible(nFx)
}
