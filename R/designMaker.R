makeDesign <- function(randFxOrder, phases, groups, propErrVar,
                       randFxVar, randFxCor, design = 'polyICT')
{
  # input check
  if(length(randFxOrder) != length(randFxVar))
  {
    stop('`randFxOrder` and `randFxVar` must be the ',
         'same length.')
  }

  # expand groups and phases into a input matrix
  inputMat <- cbind(
    expand.grid(names(phases), names(groups)),
    expand.grid(lapply(phases, length), groups))
  names(inputMat) <- c('Phase', 'Group', 'nObs', 'n')
  for(i in seq_along(randFxOrder))
  {
    inputMat[[paste('Mean', randFxOrder[i], sep='')]] <- 0
  }
  for(i in seq_along(randFxOrder))
  {
    inputMat[[paste('Var' , randFxOrder[i], sep='')]] <- randFxVar[i]
  }
  for(i in seq_along(propErrVar))
  {
    inputMat[[names(propErrVar)[i]]] <- propErrVar[i]
  }
  rm(i)
  meanNames <- names(inputMat)[grepl('Mean', names(inputMat))]
  varNames  <- names(inputMat)[grepl('Var', names(inputMat))]

  # construct randFx correlation matrices that can be
  # edited later
  randFxCorMat <- randFxCovMat <- list()
  CorMat <- matrix(randFxCor,
                   length(randFxOrder),
                   length(randFxOrder))
  diag(CorMat) <- 1
  for(p in seq_along(phases))
  {
    .np <- names(phases)[p]
    for(g in seq_along(groups))
    {
      .ng <- names(groups)[g]
      randFxCorMat[[.np]][[.ng]] <- CorMat
      .V <- inputMat[inputMat$Phase==.np &
                       inputMat$Group==.ng,
                     varNames]
      randFxCovMat[[.np]][[.ng]] <- cor2cov(CorMat,
                                            unlist(.V))
    }
  }
  rm(p, g, .np, .ng, .V)

  # populate uneditable objects
  unStdRandFxMean <- "Implementation Pending"

  # construct the fixed effects design matrix, only one
  # is needed across all combinations of group and phase
  # by using maxRandFx, any unneeded columns will be
  # ignored when constructing data
  designMat <- makeDesignMat(
    phases     = phases              ,
    phaseNames = names(phases)       ,
    maxRandFx  = length(randFxOrder) ,
    design     = design              )

  # get the number of observations
  nObs <- length(c(unlist(phases)))
  n    <- sum(groups)

  return( list(
    inputMat        = inputMat        ,
    randFxVar       = randFxVar       ,
    randFxCorMat    = randFxCorMat    ,
    randFxCovMat    = randFxCovMat    ,
    propErrVar      = propErrVar      ,
    n               = n               ,
    nObs            = nObs            ,
    groups          = groups          ,
    phases          = phases          ,
    designMat       = designMat       ,
    unStdRandFxMean = unStdRandFxMean
  ))
}