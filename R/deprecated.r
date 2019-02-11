


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