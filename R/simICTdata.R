# TODO: save data for use in SAS

# TODO: to move to the description file imports
library(MASS)
library(tidyr)

#' simICTdata - this is a user interface function for one dataset,
#' switching between different designs so the function can be updated
#' @author Stephen Tueller \email{stueller@@rti.org}

simICTdata <- function(design=polyICT$new(), randFxSeed=1, errorSeed=2,
                       family=NO)
{

  # polynomial ICT
  if('polyICT' %in% class(design))
  {
    data <- simPolyICT(design, randFxSeed, errorSeed)
  }

  # TODO see ALDA pp 234-235 for implementing these, we may need to create a
  # separate function for each ala `simPolyICT`, or we may be able to do
  # some code injection via formula

  # hyperbolic
  if('hyperICT' %in% class(design))
  {

  }

  # inverse polynomial
  if('invPolyICT' %in% class(design))
  {

  }

  # exponential
  if('expICT' %in% class(design))
  {

  }

  # negative exponential
  if('negExpICT' %in% class(design))
  {

  }

  # logistic
  if('logisticICT' %in% class(design))
  {

  }

  return(data)

}


#' simPolyICT
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#'

simPolyICT <- function(design, randFxSeed, errorSeed, family)
{

}



#' mvrFam - a function to simulate multivariate data from a gamlss family
#' distribution
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' follow http://www.econometricsbysimulation.com/2014/02/easily-generate-correlated-variables.html
mvrFam <- function(n, parms, corMat, FUN=qNO)
{
  # simple demo, move to @examples
  n <- 10000
  set.seed(1)
  yMvrnorm <- mvrnorm(n, mu=rep(0, nrow(corMat)), corMat)
  set.seed(2)
  y1rnorm <- rnorm(n, mean=0, sd=1)
  set.seed(3)
  y2rnorm <- corMat[2,1]*y1rnorm + (1-corMat[2,1])*rnorm(n, 0, 1)
  yrnorm  <- cbind(y1rnorm, y2rnorm)

  all.equal(var(yMvrnorm), var(yrnorm))



  FUN <- function(...){ FUN(...) }



  corMat <- matrix(.5, nrow=3, ncol=3) + diag(3)*.5
  corMat[1,3] <- corMat[3,1] <- .25

  n <- 1000

  Y <- mvrnorm(n, rep(0, 3), corMat)
  psych::pairs.panels(Y) # TODO make the plots not run
  all.equal(cor(Y), corMat, tolerance=.03)

  YpNorm <- pnorm(Y)
  psych::pairs.panels(YpNorm)
  all.equal(cor(YpNorm), corMat, tolerance=.08)

  YNO <- do.call(cbind, lapply(data.frame(YpNorm), FUN, mu=3, sigma=1, nu=2))
  psych::pairs.panels(YNO)
  all.equal(cor(YNO), corMat, tolerance=.08)

  YLOGNO <- do.call(cbind, lapply(data.frame(YpNorm), qLOGNO, mu=3, sigma=1))
  psych::pairs.panels(YLOGNO)
  all.equal(cor(YLOGNO), corMat, tolerance=.08)

  YBEINF <- do.call(cbind, lapply(data.frame(YpNorm), qBEINF, mu=.5, sigma=.35))
  psych::pairs.panels(YBEINF)
  all.equal(cor(YBEINF), corMat, tolerance=.08)




}


#TODO: the next step is to get the relative variances set, maybe work this
#out in the main function and the pass, do a QC step
ICTdata <- function(theDesign = getICTdesign() ,
                    errors    = ICTerror()     ,
                    randFx    = ICTrandFx()    )
{
  effectSizes <- theDesign$effectSizes
  n <- nrow(errors)

  # stack the data


  #TODO: somewhere in here, rescale everything to have the right relative
  #proportions for variances

  # construct the non-error part of the model
  yStar <-
    # random intercepts
    matrix((effectSizes$intercept + randFx$Int), n, theDesign$ntimes) +
    # random slopes
    (effectSizes$time + randFx$Slope) %*% t(theDesign$times) +
    # phase
    matrix((effectSizes$phase * theDesign$phases), n, theDesign$ntimes, byrow=TRUE) +
    # phase x time interaction
    matrix(effectSizes$phaseTime * theDesign$phases * theDesign$times,
           n, theDesign$ntimes, byrow=TRUE)

  # add errors
  y <- yStar + errors

  # create the dataset
  data <- data.frame(
    id    = sort(rep(1:n, theDesign$nObservations)) ,
    y     = c(y)                                    ,
    time  = rep(theDesign$times , n)                ,
    phase = rep(theDesign$phases, n)
  )

  # return the data
  return(data)
}

#' ICTrandFx
#' @author Stephen Tueller \email{stueller@@rti.org}
#'

# if n=1, should we be simulating values, or just use mu and rely on thy
# error covariace matrix for person to person variability? No, still simulate
# otherwise it is repeated observations of one person (or exchangeable persons).

# TODO:Document that the mean of the slopes is an effect size and should come
# from the top-level simICTdata() function
# TODO:Document that mu should be left 0 as means are added in ICTdesign()
ICTrandFx <- function(design, randFxSeed)
{
  # create the mu vector - mu currently added in ICTdata() function
  # with values stored in output of getICTdesign() function
  mu <- design$muFUN(design$effectSizes$randFx)

  # convert variance ratio to values summing to 1
  randFxVar <- (randFxVar/sum(randFxVar))

  # create the sigma matrix
  Sigma <- SigmaFUN(corMat, randFxVar)

  # rescale total to unit TOTAL variance
  Sigma <- Sigma/sum(Sigma)

  # and rescale to the varRatio
  Sigma <- theDesign$varRatio[1]*Sigma

  # QC the rescale, does not require asymptotics
  all.equal(sum(Sigma), theDesign$varRatio[1])

  # simulate the random effects
  set.seed(seed)
  randFx        <- mvrnorm(n, mu, Sigma)
  randFx        <- data.frame(randFx)
  names(randFx) <- randFxNms

  # QC - asymptotic, only true for large N
  all.equal(c(Sigma), c(var(randFx)), tolerance = .05)

  # QC
  # definitional - the varance of a sum is sum of covariance matrix
  all.equal(var(apply(randFx, 1, sum)), sum(var(randFx)))
  # should be near 1, only holds asymptotically
  all.equal(sum(var(randFx)), sum(Sigma), tolerance = .05)

  # return
  return(randFx)
}

#' ICTerror
#' @author Stephen Tueller \email{stueller@@rti.org}
#'

ICTerror <- function(n         = 9                       ,
                     parms     = list(ar=c(.5),ma=c(0))  ,
                     theDesign = getICTdesign()          ,
                     errorFUN  = arima.sim               ,
                     seed      = 1                       ,
                     ...                                     )
{
  # get seeds
  set.seed(seed)
  seeds <- as.list( ceiling(runif(n, 0, 9e6) ) )

  # sim errors, use sd = 1 for now, it can be rescaled in other functions
  FUN <- function(x, errorFUN, model, n, sd,
                  ...)
  {
    set.seed(x)
    errorFUN(model = parms, n = n, sd = sd, ...)
  }
  errors <- lapply(seeds, FUN,
                   errorFUN = errorFUN                   ,
                   model    = parms                      ,
                   n        = theDesign$nObservations ,
                   sd       = 1                          ,
                   ...)
  errors <- data.frame(do.call(rbind, errors))

  # rescale to have unit variances within time points then multiply
  # by the proportion of error variance
  errorsr <- scale(errors, FALSE, TRUE) * sqrt(theDesign$varRatio[2])

  # QC, only holds assymptotically
  all.equal(mean(apply(errorsr, 2, var)), theDesign$varRatio[2], tolerance = .05)

  # return
  return(errorsr)
}
















