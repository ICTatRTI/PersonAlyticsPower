#' Function to translate input matrices from the PersonAlytics GUI
#'
#' @export
#'
#' @param inputMat An input matrix created by the PersonAlytics GUI, which
#' may or may not have been edited by a user.
#'
#' @param translator The translator version. Currently only \code{gui1} is
#' available, more will be added if needed.
#'
#' @examples
#'
#' data(GUIinputMatExample)
#' inputMat <- GUIinputMat(GUIinputMatExample)
#'
#' myPolyICT <- polyICT$new(
#'   groups            = c(group1=10, group2=20)          ,
#'   phases            = makePhase()                      ,
#'   propErrVar        = c(randFx=.5, res=.25, mserr=.25) ,
#'   randFxOrder       = 1                                ,
#'   randFxCor         = 0.2                              ,
#'   randFxVar         = c(1, 1)                          ,
#'   error             = armaErr$new()                    ,
#'   merror            = armaErr$new(list())              ,
#'   ySD               = 15                               ,
#'   yMean             = 100                              ,
#'   )
#'
#' myPolyICT$inputMat <- inputMat
#'
#' \dontrun{
#' myPolyICT$designCheck()
#' }
#'
#'
GUIinputMat <- function(inputMat, propErrVar, translator=gui1)
{
  newInputMat <- translator(inputMat, propErrVar)

  translatorMsg <- c("Variable Check:", "Passed")
  if(!is.data.frame(newInputMat))
  {
    translatorMsg <- newInputMat
  }
  .checkInputMat(newInputMat, translatorMsg)

  return(newInputMat)
}

#' Translator for GUI output version April 2020
#'
#' @export
#'
#' @param inputMat An inputmatrix from a \code{\link{polyICT}} object.
gui1 <- function(inputMat, propErrVar)
{
  nms <- c("Phase",
           "Subgroup",
           "Number.of.Observations",
           "Sample.Size",
           "Phase.Intercept",
           "Phase.Slope",
           "Intercept.Standard.Deviation",
           "Slope.Standard.Deviation")

  isGui1 <- all(nms %in% names(inputMat))

  msg <- c("Variable Check:", "Passed")
  if(!isGui1)
  {
    msg[2] <- "Failed: One or more required variables are missing."
    return(msg)
  }
  if(isGui1)
  {
    Phase  <- inputMat$Phase
    Group  <- inputMat$Subgroup
    nObs   <- inputMat$Number.of.Observations
    n      <- inputMat$Sample.Size
    Mean0  <- inputMat$Phase.Intercept
    Mean1  <- inputMat$Phase.Slope
    Var0   <- inputMat$Intercept.Standard.Deviation^2
    Var1   <- inputMat$Slope.Standard.Deviation^2
    randFx <- propErrVar[1]
    res    <- propErrVar[2]
    mserr  <- propErrVar[3]

    newInputMat <- data.frame(
      Phase  = Phase  ,
      Group  = Group  ,
      nObs   = nObs   ,
      n      = n      ,
      Mean0  = Mean0  ,
      Mean1  = Mean1  ,
      Var0   = Var0   ,
      Var1   = Var1   ,
      randFx = randFx ,
      res    = res    ,
      mserr  = mserr  )

    return(newInputMat)
  }
}

#' @keywords internal
.checkInputMat <- function(inputMat, translatorMsg)
{
  dfMsg   <- c("Data Type Check:", ifelse(is.data.frame(inputMat), "Passed",
                                         "Failed: inputMat is not a data.frame"))

  sampMsg <- c("Group Size Check:", "Could not be run")
  if(is.data.frame(inputMat))
  {
    sampMsg <- .checkSampSize(inputMat)
  }

  msgs <- data.frame( rbind(translatorMsg, dfMsg, sampMsg) )
  names(msgs) <- c("Check", "Status")

  write.csv(msgs, "inputMat_Checks.csv")

  if(!all(msgs$Status == "Passed"))
  {
    stop("\n\nOne or more checks on the `inputMat` failed, see the file",
         "\n'inputMat_Checks.csv' in the current working directory",
         "\n(type `getwd()` in the console to see the working directory).")
  }
}

#' @keywords internal
.checkSampSize <- function(inputMat)
{
  groups <- levels(factor(inputMat$Group))
  group_size_variance <- list()
  for(i in groups)
  {
    group_size_variance[[i]] <- var(inputMat$n[inputMat$Group==i])
  }
  group_size_variance <- unlist(group_size_variance)

  msg <- c("Group Size Check:", "Passed")

  if(any(group_size_variance != 0))
  {
    msg[2] <- "Failed: One or more groups have sample sizes that vary over phase."
  }

  msg
}

