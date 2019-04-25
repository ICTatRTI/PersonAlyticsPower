# this file will contain the the *Data functions
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

polyData <- function(randFx, errors, ymean=NULL, yvar=NULL)
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