# these may be used in the future


#' readinteger - helper function for studySetup
#' @author adapted from http://www.rexamples.com/4/Reading%20user%20input
#' @keywords internal

readinteger <- function(prompt, nth="")
{
  n <- readline(prompt=paste(prompt, " ", nth, " : ", sep=""))
  if(!grepl("^[0-9]+$",n))
  {
    return(readinteger())
  }
  return(as.integer(n))
}

readnumeric <- function(prompt, nth="")
{
  n <- readline(prompt=paste(prompt, " ", nth, " : ", sep=""))
  n <- as.numeric(n)
  if(!is.numeric(n))
  {
    return(readinteger())
  }
  return(n)
}

readcharacter <- function(prompt, nth="")
{
  string <- readline(prompt=paste(prompt, " ", nth, " : ", sep=""))
  if(!is.character(string))
  {
    return(readcharacter())
  }
  return(as.character(string))
}
