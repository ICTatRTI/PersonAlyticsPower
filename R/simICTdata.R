# TODO: save data for use in SAS

# TODO: to move to the description file imports
library(MASS)
library(tidyr)

#' simICTdata
#' @author Stephen Tueller \email{stueller@@rti.org}

simICTdata <- function(design, errors)
{

  # get the variances from the design inputs
  theDesign <- getICTdesign()

  # simulate the errors
  errors <- ICTerror()

  # simulate the random effects
  randFx <- ICTrandFx()

  # construct the observed data

}

#' ICTdata
#' @author Stephen Tueller \email{stueller@@rti.org}
#'


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
ICTrandFx <- function(n         = 9                          ,
                      theDesign = getICTdesign()             ,
                      mu        = c(0,0)                     ,
                      corMat    = matrix(c(1,.2,.2,1), 2, 2) ,
                      randFxVar = c(1, .1)                   ,
                      randFxNms = c('Int', 'Slope')          ,
                      muFUN     = function(x) x              ,
                      SigmaFUN  = cor2cov                    ,
                      seed      = 1                          )
{
  # create the mu vector - mu currently added in ICTdata() function
  # with values stored in output of getICTdesign() function
  mu <- muFUN(mu)

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


#' getICTdesign
#' @author Stephen Tueller \email{stueller@@rti.org}

# currently limited to linear model, this may be a beast to generalize
getICTdesign <- function(varRatio      = c(.5,.5)              ,
                         effectSizes   = list(intercept =  0,
                                              phase     = .5,
                                              time      = .5,
                                              phaseTime = .2)  ,
                         propBaseline  = .25                   ,
                         nObservations = 20                    )
{
  # check the variance ration
  if(sum(varRatio)!=1) stop('`varRatio` must sum to 1;',
                            '`varRatio` is the ratio of between variance to error variance.')

  # generate the times
  times <- seq(0, nObservations-1, 1)
  ntimes <- length(times)
  varTime <- var(times)

  # generate phases
  nBaseline <- ceiling(propBaseline*nObservations)
  phases <- list(phase1=rep(0, nBaseline),
                 phase2=rep(1, nObservations - nBaseline))
  phases <- unlist(phases)
  varPhase <- var(phases)

  # generate interaction
  phaseTime <- phases * times
  varPhaseTime <- var(phaseTime)

  # design linear predictors
  Design <- data.frame(times=effectSizes$time * times,
                       phases=effectSizes$phase * phases,
                       phaseTime=effectSizes$phaseTime * phaseTime)
  DesignTotal  <- apply(Design, 1, sum)
  varDesignTotal <- var(DesignTotal)

  # QC - does the variance of the total var(Total) equal the sum of the variances
  # plus 2 X the covariances?
  varDesign <- var(Design)
  all.equal(varDesignTotal, sum(diag(varDesign)) + sum(2 * varDesign[lower.tri(varDesign)]))

  # QC against inputs
  all.equal(varDesign[1,1], effectSizes$time^2  * varTime)
  all.equal(varDesign[2,2], effectSizes$phase^2 * varPhase)
  all.equal(varDesign[3,3], effectSizes$phaseTime^2 * varPhaseTime)
  all.equal(varDesign[2,1], effectSizes$phase * effectSizes$time * cov(phases, times))
  all.equal(varDesign[3,1], effectSizes$phaseTime * effectSizes$time * cov(phaseTime, times))
  all.equal(varDesign[3,2], effectSizes$phaseTime * effectSizes$phase * cov(phaseTime, phases))

  # return
  return(list(varRatio      = varRatio       ,
              effectSizes   = effectSizes    ,
              propBaseline  = propBaseline   ,
              nObservations = nObservations  ,
              nBaseline     = nBaseline      ,
              ntimes        = ntimes         ,
              times         = times          ,
              phases        = phases         ,
              phaseTime     = phaseTime      ,
              varTime       = varTime        ,
              varPhase      = varPhase       ,
              varPhaseTime  = varPhaseTime   ,
              varDesign     = varDesignTotal ,
              DesignTotal   = DesignTotal    ))
}













