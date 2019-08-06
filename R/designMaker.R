#' makeDesign -
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
makeDesign <- function(randFxOrder, phases, groups, propErrVar,
                       randFxVar, randFxCor, design = 'polyICT',
                       makeInputMat = TRUE, self = NULL, isNew = NULL)
{
  # input check
  if(length(randFxOrder) != length(randFxVar))
  {
    stop('`randFxOrder` and `randFxVar` must be the ',
         'same length.')
  }

  # expand groups and phases into a input matrix
  if( makeInputMat                                                    |
      'randFxOrder' %in% isNew                                        |
      ('phases' %in% isNew & (length(phases) != length(self$phases))) |
      ('groups' %in% isNew & (length(groups) != length(self$groups))) |
      ('randFxCor' %in% isNew)
    )
  {
    inputMat <- cbind(
      expand.grid(names(phases), names(groups)),
      expand.grid(unlist(lapply(phases, length)), groups))
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

  }

  else
  {
    if( 'propErrVar' %in% isNew )
    {
      wc <- which( names(self$inputMat) %in% c('randFx', 'res', 'mserr') )
      self$inputMat[,wc] <- matrix(propErrVar, nrow(self$inputMat), 3, byrow=TRUE)
    }
    if( 'randFxVar' %in% isNew )
    {
      wc <- which( names(self$inputMat) %in% self$varNames )
      self$inputMat[,wc] <- matrix(randFxVar, nrow(self$inputMat), length(wc),
                                   byrow=TRUE)
    }
    if( 'phases' %in% isNew )
    {
      .phases <- unlist( lapply(phases, length) )
      for(p in seq_along(.phases))
      {
        self$inputMat[self$inputMat$Phase==names(.phases)[p],'nObs'] <- .phases[p]
      }
    }
    if( 'groups' %in% isNew )
    {
      .groups <- unlist( lapply(groups, length) )
      for(g in seq_along(.groups))
      {
        self$inputMat[self$inputMat$Group==names(.groups)[g],'n'] <- groups[g]
      }
    }
    inputMat <- self$inputMat
  }
  meanNames <- names(inputMat)[grepl('Mean', names(inputMat))]
  varNames  <- names(inputMat)[grepl('Var', names(inputMat))]

  # construct randFx correlation matrices that can be
  # edited later
  randFxCorMat <- list()
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
    }
  }
  rm(p, g, .np, .ng, .V)

  # construct the fixed effects design matrix, only one
  # is needed across all combinations of group and phase
  # by using randFxOrder, any unneeded columns will be
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
    propErrVar      = propErrVar      ,
    n               = n               ,
    nObs            = nObs            ,
    groups          = groups          ,
    phases          = phases          ,
    designMat       = designMat       ,
    meanNames       = meanNames       ,
    varNames        = varNames
  ))
}


#' makeDesignMat - create the fixed effects design matrix for n=1
#' @author Stephen Tueller \email{stueller@@rti.org}
#' @export
makeDesignMat <- function(phases      = makePhase() ,
                          phaseNames  = NULL        ,
                          maxRandFx   = 2           ,
                          design      = 'polyICT'
)
{
  # make names if null
  if(is.null(phaseNames)) phaseNames <- paste('phase', 1:length(phases))

  # get the number of observations
  nObservations <- length(c(unlist((phases))))

  if(design == 'polyICT')
  {
    # generate the times
    times      <- list()
    times[['Time']] <- seq(0, nObservations-1, 1)
    if(maxRandFx>1)
    {
      for(i in 2:maxRandFx)
      {
        times[[paste('Time', i, sep='')]] <- times[['Time']]^i
      }
    }
    time    <- data.frame( do.call(cbind, times) )
    #varTime <- apply(time, 2, var)

    # clean up phases
    phase <- as.numeric(factor( c(unlist(phases)) ) ) - 1
    phase <- factor(phase, labels=phaseNames)

    designMat <- cbind(phase, time)

    return(designMat)
  }

}