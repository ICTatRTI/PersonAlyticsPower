#' @title polyICT2 class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
#'
#' @field n Numeric (integer). The number of participants. Default is 10.
#'
#' @field phases List. Each phase in one item in the list with the phase name
#' repeated for the number of time points in the phase. For example, an "ABA"
#' study with 5 time points each would be \code{list(rep("A", 5), rep("B", 5),
#' rep("A", 5))}. See also the function \code{\link{makePhase}}. Default is
#' \code{makePhase()}.
#'
#' @field propErrVar Numeric. The propotion of total variance that is error
#' variance. Default is .75.
#'
#' @field randFxMean List of lists of named numeric vectors. The general form
#' is as follows, where the elipses (...) only illustrate that additional
#' inputs could be given:
#'
#' \code{
#' randFxMean = list(
#'   group1 = list(
#'     phase1 = c(i=0.0, s=0.0, q=0.0, ...),
#'     phase2 = c(i=0.2, s=0.0, q=0.0, ...),
#'     phase3 = c(i=0.0, s=0.0, q=0.0, ...),
#'     ...
#'     ),
#'   group2 = list(
#'     phase1 = c(i=0.0, s=0.0, q=0.0, ...),
#'     phase2 = c(i=0.0, s=0.0, q=0.0, ...),
#'     phase3 = c(i=0.0, s=0.0, q=0.0, ...),
#'     ...
#'     ),
#'   ...
#' )
#' }
#'
#' The length of randFxMean \code{length(randFxMean)} is the number of groups
#' and can take on any arbitrary name without quotes as long as it is a valid
#' variable name, see \code{\link{make.names}}. In the example above, the group
#' names are `group1` and `group2`, given without quotes.
#'
#' Each group is itself a list whose length is the number of phases. The phase
#' names can take on any arbitrary name without quotes as long as they are valid
#' variable names. In the example above, \code{length(randFxMean$group1)} is 3,
#' i.e., there are three phases, and they are named `phase1`, `phase2`,
#' and `phase3`.
#'
#' Each phase is a named numeric vector of effect sizes on the scale of Cohen's
#' d. In the example above, `i` indicates intercepts, `s` are slopes, and `q`
#' are quadratic terms. The elipsis (...) indicates higher order terms could
#' be included. All `q` and `s` terms are zero inticating no change over time.
#' These could be left out and only `i` included. Since \code{i=0.2} in the
#' `phase2` of `group1`, this indicates a small increase of Cohen's d=0.2 during
#' `phase2` relative to `phase1` (and `phase3`) in `group1` and relative to all
#' phases in `group2`. I other words, this is an ABA design with intervention
#' only in `phase2` for `group1`.
#'
#' @field randFxCorMat Numeric matrix. A symmetric correlation matrix with a dimension
#' equal to the order of the model. For example, a quadratic model would
#' correpspond to a 3x3 matrix. The diagonal elements must equal 1, the off
#' diagonal elements must be between -1 and +1, and the matrix must be
#' invertable.
#'
#' @field randFxVar Numeric vector. A vector of the same length as the order of
#' the polynomial
#' model containing the variances of the random effects. For example, in a
#' quadratic model, `randFxVar` would be length three, the first element would
#' be the intercept variance, the second element would be the slope variance,
#' and the third element would be the variance of the quadratic term.
#'
#' @field SigmaFun An R function. A function to convert \code{randFxCorMat} and
#' \code{randFxVar} into the covariance matrix \code{randFxCovMat}. Default is
#' \code{\link{cor2cov}}.
#'
#' @return
#'
#' \item{n               }{The number of participants (see `Fields`).}
#' \item{phases          }{The phases of the study (see `Fields`).}
#' \item{propErrVar      }{The proportion of error variance (see `Fields`).}
#' \item{randFxMean     }{The fixed effects effect sizes (see `Fields`).}
#' \item{randFxCorMat          }{The correlation matrix of the random effects (see `Fields`).}
#' \item{randFxVar       }{The variance of the random effects (see `Fields`).}
#' \item{muFUN           }{A function for transforming random effects means (see `Fields`).}
#' \item{SigmaFun        }{A function for constructing the covariance matrix (see `Fields`).}
#' \item{maxRandFx       }{The order of the study. For example, a linear model
#' would be of order 1, and a quadratic model would be order 2.}
#' \item{designMat       }{The design matrix for one participant showing the structure
#' of study timing and phases.}
#' \item{randFxCovMat          }{The covariance matrix constructed using \code{randFxCorMat},
#' \code{randFxVar}, and \code{SigmaFun}.}
#' \item{nObservations   }{The number of observations per participant.}
#' \item{variances       }{Partition of the total variance into that due to
#' random effects and that due to error variance at time = 1.}
#' \item{expectedVariance}{The expected variances across all time points. This will
#' not match the variance of the simulated data unless n is large. See the \code{checkDesign}
#' parameter in \code{\link{ICTpower}}.}
#' \item{unStdRandFxMean    }{The unstandardized effect sizes constructed from the total
#' \code{expectedVariance} and the \code{randFxMean}.}

#'
#' @examples
#'
#' # produce a simple ICT design
#' defaultPolyICT <- polyICT$new()
#'
#' # print a summary
#' defaultPolyICT
#'
#' # view the fields that are generated by `$new()` but cannot be changed by
#' # the user
#' defaultPolyICT$maxRandFx
#' defaultPolyICT$designMat
#' defaultPolyICT$randFxCovMat
#' defaultPolyICT$nObservations
#' defaultPolyICT$variances
#' defaultPolyICT$expectedVariances
#' defaultPolyICT$unStdRandFxMean
#'
#'
#'

polyICT2 <- R6::R6Class("polyICT",

     inherit = designICT,

     private = list(
       .n                 = NULL,
       .phases            = NULL,
       .propErrVar        = NULL,
       .randFxMean       = NULL,
       .randFxCorMat            = NULL,
       .randFxVar         = NULL,
       .muFUN             = NULL,
       .SigmaFun          = NULL,
       .maxRandFx         = NULL,
       .designMat         = NULL,
       .randFxCovMat            = NULL,
       .nObservations     = NULL,
       .variances         = NULL,
       .expectedVariances = NULL,
       .unStdRandFxMean      = NULL
     ),

     public  = list(
       initialize = function
       (
         n                 = 10                                        ,
         phases            = makePhase()                               ,
         propErrVar        = .75                                       ,
         randFxMean       = NULL                                      ,
         randFxCorMat            = matrix(c(1,.2,.2,1), 2, 2)                ,
         randFxVar         = c(1, .1)                                  ,
         muFUN             = function(x) x                             ,
         SigmaFun          = cor2cov                                   ,
         maxRandFx         = NULL                                      ,
         designMat         = NULL                                      ,
         randFxCovMat            = NULL                                      ,
         nObservations     = NULL                                      ,
         variances         = NULL                                      ,
         expectedVariances = NULL                                      ,
         unStdRandFxMean      = NULL
       )
       {
         if(is.null(randFxMean))
         {
           randFxMean = list(
             group1 = list(
               phase1 = c(i=0.0, s=0.0, q=0.0),
               phase2 = c(i=0.2, s=0.0, q=0.0),
               phase3 = c(i=0.0, s=0.0, q=0.0)
             ),
             group2 = list(
               phase1 = c(i=0.0, s=0.0, q=0.0),
               phase2 = c(i=0.0, s=0.0, q=0.0),
               phase3 = c(i=0.0, s=0.0, q=0.0)
             )
           )
         }

         # general input validation
         # TODO: checkPolyICT2 needs to be updated for the new definition of
         # randFxMean
         #maxRandFx <- checkPolyICT2(randFxMean, randFxCorMat, randFxVar) - 1

         # construct the fixed effects design matrix
         designMat <- getdesignICT(phases, maxRandFx, 'polyICT')

         # get the covariance matrix
         randFxCovMat <- cor2cov(randFxCorMat, randFxVar)

         # get the number of observations
         nObservations <- length(c(unlist(phases)))

         # get the total variance and the error variance @ time = 1
         totalVar  <- sum(randFxCovMat)/(1-propErrVar)
         errorVar  <- propErrVar * totalVar
         variances <- list(totalVar = totalVar,
                           errorVar = errorVar,
                           randFxVar = totalVar - errorVar)

         # get the expected varariances
         # initial values
         expectedVariances <- expectedVar(randFxCovMat, designMat, variances,
                                          randFxMean, nObservations, n)
         unStdRandFxMean      <- expectedVariances$randFxMean

         # populate private
         private$.n                 <- n
         private$.phases            <- phases
         private$.propErrVar        <- propErrVar
         private$.randFxMean       <- randFxMean
         private$.randFxCorMat            <- randFxCorMat
         private$.randFxVar         <- randFxVar
         private$.muFUN             <- muFUN
         private$.SigmaFun          <- SigmaFun
         private$.maxRandFx         <- maxRandFx
         private$.designMat         <- designMat
         private$.randFxCovMat            <- randFxCovMat
         private$.nObservations     <- nObservations
         private$.variances         <- variances
         private$.expectedVariances <- expectedVariances
         private$.unStdRandFxMean      <- unStdRandFxMean

       },

       print = function(...)
       {
         # use super to call the print function of the parent class
         super$print()


         # add print options for the child class
         cat("\n\nEffect sizes:\n",
             paste(unlist(lapply(self$randFxMean, names)), '=',
                   unlist(self$randFxMean), '\n'),
             "\nRandom Effects\n",
             "\nCorrelation matrix:\n", catMat(self$randFxCorMat),
             "\nCovariance matrix:\n",  catMat(round(self$randFxCovMat,3))
         )

       },

       makeData = function(randFx, errors=NULL, y=NULL, ymean=NULL, yvar=NULL)
       {
         # TODO: consider having filename option and saving here

         # determine whether to build or deconstruct

         # we fail
         if( !is.null(y) & !is.null(errors) )
         {
           stop("Only one of `y` or `errors` should be provided, not both.")
         }

         # we build y
         if(!is.null(randFx) & !is.null(errors))
         {

         }

         # we deconstruct y
         if(!is.null(randFx) & !is.null(errors))
         {

         }

         # return the data
         invisible(data)

       }

     )
)

# This approach takes simulated y and subtracts out its components
# but...do we really need to do this? We already have Y. When we build Y,
# we don't save the random effects or errrors, just Y. I guess we could
# deconstruct Y as a check for user inputs. For example, we don't really
# know what the effect sizes are, or the error variances

#' deconstructY - a function to deconstruct Y into
#' random effects and error
deconstructY <- function(randFx, y, ymean, yvar)
{

  # step 1: simulate y

  # step 2: sample from the specified regions of y

  # step 3:
}

# This approach builds y from its components
#'
buildY <- function(randFx, errors, ymean, yvar)
{
  # make each component into a matrix of dimension
  # n by nObservations

  # fixed and random effects
  fe <- re <- list()
  fev <- rev <- list()
  for(i in seq_along(self$unStdRandFxMean[[1]]))
  {
    # fixed effects
    feName  <- names(self$unStdRandFxMean$fixdFx)[i]
    fe[[i]] <- matrix((self$unStdRandFxMean$fixdFx[[feName]] *
                         self$designMat[[feName]]),
                      self$n, self$nObservations, byrow=TRUE)

    # intercepts
    if(i==1)
    {
      re[[i]] <- matrix((self$unStdRandFxMean$randFx[[i]] + randFx[,i]),
                        self$n, self$nObservations)
    }
    # slopes and higher
    if(i>1)
    {
      re[[i]] <- matrix(self$unStdRandFxMean$randFx[[i]] + randFx[,i]) %*%
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

# TODO:need to generalize beyond slopes
#' expectedVar
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
expectedVar <- function(randFxCovMat, designMat, variances,
                        randFxMean, nObservations, n)
{
  # 'total' N
  N <- nObservations * n

  # variance at each time point
  vyt <-  randFxCovMat[1,1]            +
    designMat$Time^2*randFxCovMat[2,2] +
    2*randFxCovMat[1,2]                +
    variances$errorVar

  # total variance due to random effects and error variance
  sigmaT <- (N-1)^-1 * ( (n - 1) * sum(vyt) )

  # rescale the fixed effect sizes in terms of the variance when time = 1
  whereIsTeq1 <- which(designMat$Time==1)
  sdYt1 <- sqrt( vyt[whereIsTeq1] )
  randFxMean$fixdFx$phase     <- randFxMean$fixdFx$phase     * sdYt1
  randFxMean$fixdFx$phaseTime <- randFxMean$fixdFx$phaseTime * sdYt1

  # mean at each time point
  mt <-  (randFxMean$randFx$intercept                       +
          randFxMean$randFx$slope     * designMat$Time      +
          randFxMean$fixdFx$phase     * designMat$phase     +
          randFxMean$fixdFx$phaseTime * designMat$phaseTime
  )

  # overall mean
  m = mean(mt)

  # total variance due to fixed effects stacked over time
  muT <- (N-1)^-1 * ( n * sum((mt-m)^2) )

  # overall variance
  vy     <- sigmaT + muT
  #  (N-1)^-1 * ( n * sum((mt-m)^2) + (n - 1) * sum(vyt) )


  return( list(randFxMean = randFxMean ,
               Variances   = vyt         ,
               Means       =  mt         ,
               TotalMean   =   m         ,
               TotalVar    =  vy         ) )
}

#' checkPolyICT2
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
# TODO: update for the new definition of randFxMean
checkPolyICT2 <- function(randFxMean, randFxCorMat, randFxVar)
{
  # check conformity of the elements in `randFxMean`
  randFxMean <<- randFxMean
  nFx <- checkRandFxMeanPoly(randFxMean)

  # check that `randFxCorMat` is a correlation matrix
  checkCorMat(randFxCorMat)

  # check conformity of `randFxCorMat` with `randFxMean`
  test3 <- nrow(randFxCorMat) == nFx[[1]]
  if(!test3) stop('The number of rows and columns in `randFxCorMat` is ', nrow(randFxCorMat),
                  '\nbut it should equal the ', nFx[[1]], ' elements in\n',
                  ' `randFxMean$randFx` and `randFxMean$fixdFx`.')

  # check conformity of `randFxVar` with `randFxMean`
  test4 <- length(randFxVar) == nFx[[1]]
  if(!test4) stop('The number of elements in `randFxVar` is ', length(randFxVar),
                  '\nbut it should equal the ', nFx[[1]], ' elements in\n',
                  ' `randFxMean$randFx` and `randFxMean$fixdFx`.')

  # check that resulting covariance matrix is legit
  checkCorMat(cor2cov(randFxCorMat, randFxVar), FALSE)

  # return the polynomial order
  invisible( unname(unlist(nFx[1])) )
}

#' checkRandFxMean
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
checkRandFxMeanPoly <- function(randFxMean)
{
  # check conformity of the elements in `randFxMean`
  FxNames <- c( "randFx", "fixdFx")
  test1   <- all( names(randFxMean) %in% FxNames)
  if(!test1) stop('`randFxMean` must be a length 2 list with the names ',
                  paste('`', FxNames, '`', collapse=', ', sep=''))
  nFx   <- lapply(randFxMean, length)
  test2 <- nFx[[1]]==nFx[[2]]
  if(!test2) stop('`randFxMean$randFx` must have the same number of elements\n',
                  'as `randFxMean$fixdFx`.')
  invisible(nFx)
}
