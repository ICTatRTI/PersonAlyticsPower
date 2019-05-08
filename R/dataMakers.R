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


# TODO document @param
#' polyData - function to simulate polynomial growth data, used by the polyICT
#' class
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
polyData <- function(seed=123, n, nObs, mu, Sigma, self, dM,
                     rFxVr, propErrVar, group=NULL)
{
  # get seeds
  # 1. random effects
  # 2. autocorrelated errors
  # 3. white noise errors
  seeds <- .makeSeeds(seed, 3)

  # random effects
  reVars <- rFxVr/sum(rFxVr)
  eta <- makeEta(n, mu, Sigma, seeds[1])

  # rescale random effects by  by propErrVar and reVars, retaining means
  eta <- scale(eta) %*% diag( sqrt(propErrVar[1] * reVars) ) +
    matrix(apply(eta, 2, mean), nrow(eta), ncol(eta), byrow=TRUE)

  # qc - should hold even in small samples
  #all(propErrVar[1] * reVars, apply(eta, 2, var))


  # expand random effects into matrices
  .L  <- list()
  for(i in 1:ncol(eta))
  {

    if(i==1) etaM <- matrix(eta[,i]) %*% rep(1, nObs)
    if(i>1)
    {
      # multiply by time rescaled to 0,1
      dMi  <- dM[,i]-min(dM[,i])
      dMi  <- c(dMi/max(dMi))
      etaM <- matrix(eta[,i]) %*% dMi
    }
    .L[[i]] <- etaM
  }

  # TODO implement transforming error distn via makeFam()
  # errors w/ scaling
  .L[[length(.L) + 1]] <- sqrt(propErrVar[2]) * scale(self$error$makeErrors(n, nObs, seeds[2]))
  .L[[length(.L) + 1]] <- sqrt(propErrVar[3]) * scale(self$merror$makeErrors(n, nObs, seeds[3]))

  # TODO consider rescaling to unit variance first, otherwise propErrVar will
  # be off

  # now put it all together
  Y <- Reduce('+', .L)

  # create id and prepend group number to ensure unique ids
  id <- sort(rep(1:n, nObs))
  if(!is.null(group))
  {
    groupNumber <- as.numeric( gsub("[^[:digit:]]", "", group) )
    groupNumber <- (10*groupNumber)^max(nchar(id))
    id <- groupNumber + id
  }

  # construct a long dataset
  data <- data.frame(
    id    = id                                       ,
    y     = as.vector(t(Y))                          ,
    do.call(rbind, replicate(n, dM, simplify = FALSE))
  )

  # QC
  all.equal(aggregate(data$y, list(data$Time), mean)$x,
            unname(apply(Y, 2, mean)))
  all.equal(aggregate(data$y, list(data$Time), var)$x,
            unname(apply(Y, 2, var)))

  # make phase a factor (for later analysis and plotting)
  data$phase <- factor(data$phase)

  # add group identifiers, phase is already coming in via dM
  if(!is.null(group)) data$group <- group

  # return the data
  return(data)
}



