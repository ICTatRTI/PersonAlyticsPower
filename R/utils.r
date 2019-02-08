

#' cor2cov - a version of this function was in the deprecated package `stremo`
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

cor2cov <- function(corMat    = matrix(c(1,.2,.2,1), 2, 2) ,
                    variances = c(1, .1)                   )
{
  sds <- sqrt(variances)
  b <- sds %*% t(sds)
  Sigma <- b * corMat
  # if QC is wanted use
  # all.equal(cov2cor(Sigma), corMat)
  return(Sigma)
}

#' makePhase - a function to create the phase variable
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal

makePhase <- function(nObsPerCondition = c(10,20,10) ,
                       conditions = c("A", "B", "A")  )
{
  phases <- list()
  for(i in seq_along(nObsPerCondition))
  {
    phases[[i]] <- rep(conditions[i], nObsPerCondition[i])
  }
  return( phases )
}

#' checkPolyICT
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'

checkPolyICT <- function(effectSizes, corMat, randFxVar)
{
  # check conformity of the elements in `effectSizes`
  checkEffectSizesPoly(effectSizes)

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
}

#' checkCorMat
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
checkCorMat <- function(corMat, cor=TRUE)
{
  mNm <- 'corMat'
  if(!cor) mNm <- 'varMat'
  if(!is.matrix(corMat)) stop('`corMat` must be a numeric matrix.')
  if( is.matrix(corMat))
  {
    isNum <- is.numeric(value)
    isSym <- isSymmetric(value)
    isCor <- all(diag(value)==1)
    isSlv <- class(try(solve(value), silent = TRUE)) %in% 'matrix'
    if(!isNum) stop('`', mNm, '` is not numeric')
    if(!isSym) stop('`', mNm, '` is not symmetric')
    if(cor & !isCor) stop('`corMat` does not have ones on the diagonal')
    if(!isSlv) stop('`', mNm, '` cannot be inverted')
  }
}
