# this file will contain the the *Data functions [consider moving to repsective *ICT files]
#
#  polyData
#  hyperData
#  invPolyData
#  expData
#  negExpData
#  logisticData
#
#  that will be used by
#
#  polyICT$makeData()
#  hyperICT$makeData()
#  invPolyICT$makeData()
#  expICT$makeData()
#  negExpICT$makeData()
#  logisticICT$makeData()
#
# These functions do not simulate data, rather they construct
# latent growth curve model data for one group within one phase
# from inputs of random effects and errrors.
#
# Mutligroup and/or multiphase data are achieved by concatenating data
# from multiple calls to *Data functions within *makeData methods.

# level 2 covariates -
# -- categorical: acheived via multiple groups
# -- continuous: implemented directly in lgmSim

# time varying covariates -
# -- continuous: implemented directly in lgmSim
# -- categorical: achieve by categarizing continuous covariates

polyData <- function(seed=123, n, nObs, mu, Sigma, self, dM,
                     rFxVr, propErrVar, yMean, ySD, group=NULL)
{
  # get seeds
  # 1. random effects
  # 2. autocorrelated errors
  # 3. white noise errors
  seeds <- makeSeeds(seed, 3)

  # random effects
  reVars <- rFxVr/sum(rFxVr)
  eta <- mvrFam(n, mu, Sigma, seeds[1])
  .L  <- list()
  for(i in 1:ncol(eta))
  {
    # rescale data by propErrVar and reVars and expand into matrix
    etaM <- matrix(propErrVar[i] * reVars[i] * eta[,i], n, nObs)
    if(i>1)
    {
      etaM <- etaM * c(dM[,i+1])
    }
    .L[[i]] <- etaM
  }

  # TODO implement transforming error distn via makeFam()
  # errors w/ scaling
  .L[[length(.L) + 1]] <- propErrVar[2] * ICTerror(self$error, n, nObs, seeds[2])
  .L[[length(.L) + 1]] <- propErrVar[3] * ICTerror(self$merror, n, nObs, seeds[3])

  # TODO consider rescaling to unit variance first, otherwise propErrVar will
  # be off

  # now put it all together
  Y <- Reduce('+', .L)

  # construct a long dataset
  data <- data.frame(
    id    = sort(rep(1:n, nObs)) ,
    y     = c(t(Y))                                 ,
    do.call(rbind, replicate(n, dM, simplify = FALSE))
  )

  # QC
  all.equal(aggregate(data$y, list(data$Time), mean)$x,
            unname(apply(Y, 2, mean)))
  all.equal(aggregate(data$y, list(data$Time), var)$x,
            unname(apply(Y, 2, var)))

  # make phase a factor (for later analysis and plotting)
  data$phase <- factor(data$phase)

  # rescale the y variance if !is.null
  if(!is.null(ySD) | !is.null(yMean))
  {
    if(is.null(yMean)) yMean <- mean(data$y)
    if(is.null(ySD))   ySD   <- sd(data$y)
    data$y <- scale(data$y, yMean, TRUE) * ySD
  }

  # add group identifiers, phase is already coming in via dM
  if(!is.null(thisg)) data$group <- group

  # return the data
  return(data)
}



