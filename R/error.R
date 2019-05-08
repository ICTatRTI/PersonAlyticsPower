
#' errActive Active bindings for Err and its heirs
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
errActive <- function()
{
  list(

    parms = function(value)
    {
      if( missing(value) ){ private$.parms }
      else
      {
        private$.parms <- value
        self
      }
    },

    errVar = function(value)
    {
      if( missing(value) ){ private$.errVar }
      else
      {
        if(!length(errVar)==1 | !is.numeric(errVar))
        {
          stop('`errVar` must be a single number instead of ', value, '.')
        }
        private$.errVar <- value
        self
      }
    },

    # I think there is a way to self-referentially get the class type, which
    # we could pass here
    FUN = function(value)
    {
      if( missing(value) ){ private$.FUN }
      else
      {
        stop('The value for `FUN` is ', private$.FUN, ' and cannot be changed.')
      }
    },

    fam = function(value)
    {
      if( missing(value) ){ private$.fam }
      else
      {
        if( ! is(gamlss.dist::gamlss.family(value)) == "gamlss.family")
        {
          stop("The value of `fam`=", value, " is not a `gamlss.family` distribution.")
        }
        private$.fam <- value
        self
      }
    },

    famParms = function(value)
    {
      if( missing(value) ){ private$.famParms }
      else
      {
        if(is.null(names(value))) names(value) <- c('mu', 'sigma', 'nu', 'tau')
        if(length(value) > 4 |
           any( ! names(value) %in%  c('mu', 'sigma', 'nu', 'tau') ))
        {
          stop("`famParms` should be a named list of length 1 to 4 with possible names\n",
               "'mu', 'sigma', 'nu', or 'tau'.")
        }
        private$.famParms <- value
        self
      }
    }

  )
}


#' \code{Err} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
#'
#' @field parms The parameters for \code{FUN}.
#'
#' @field errVar The error variance.
#'
#' @field FUN The function for simulating error terms.
#'
#' @field fam A random variable generator from the
#' \code{\link{gamlss.family}} for the error term distribution.
#' After data are simulated using \code{FUN} and its \code{parms}, it can be
#' transformed to have a distribution implied by \code{fam}.
#'
#' @field famParms Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.
Err <- R6::R6Class("err",

         private = list(
           .parms    = NULL,
           .errVar   = NULL,
           .FUN      = NULL,
           .fam      = NULL,
           .famParms = NULL
         ),

         active = errActive(),

         public = list(

           initialize = function(
              parms    = NULL,
              errVar   = 1   ,
              FUN      = NULL,
              fam      = NULL,
              famParms = NULL
           )
           {
             # validate inputs

             # populate private
             private$.parms    <- parms
             private$.errVar   <- errVar
             private$.FUN      <- FUN
             private$.fam      <- fam
             private$.famParms <- famParms
           },

           print = function(...)
           {
             print(self$parms)
             print(self$errVar)
             print(self$fam)
             print(self$famParms)
           }
         )



       )

#' \code{armaErr} class generator using \code{\link{arima.sim}}
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#'
#'
#' @field parms The parameters for \code{FUN}.
#'
#' @field errVar The error variance.
#'
#' @field FUN The function for simulating error terms.
#'
#' @field fam A \code{\link{gamlss.family}} family for the error term distribution.
#' After data are simulated using \code{FUN} and its \code{parms}, it can be
#' transformed to have a distribution different from those available to \code{FUN}.
#'
#' @field famParms Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.
#'
#' @examples
#'
#' # show that r* functions from gamlss.family distribution functions can be
#' # passed to arima.sim
#' errBeta <- arima.sim(list(ar=c(.5)), 1000, innov = rBE(1000, mu=.5, sigma=.2))
#'
#' # note that the rand.gen option only allows teh defaults for the function
#' # passed to rand.gen
#' #errBeta <- arima.sim(list(ar=c(.5)), 1000, rand.gen = rBE)
#' # that said, rand.gen's parameters can be passed via ...
#' #errBeta <- arima.sim(list(ar=c(.5)), 1000, rand.gen = rBE, mu=.9, sigma=.2)
#'
#' # note that even though beta innovations are used, the resulting data is not
#' # constrained to be in (0,1)
#' hist(errBeta)
#' plot(errBeta)
#' auto.arima(ts(errBeta))
armaErr <- R6::R6Class("errARMA",

  inherit = Err,

  public = list(
    initialize = function(
      parms    = list(ar=c(.5), ma=c(0)) ,
      fam      = "NO"                    ,
      famParms = list(mu=0, sigma=1)
    )
    {
      # validation
      if( !is.list(parms) )
      {
        if( ! all(names(parms) %in% c('ar', 'ma')) )
        {
          stop('\n`parms` must be a named list of length one or two with the\n',
             'names `ar` and/or `ma`; `parms` can also be an empty list\n',
             'to generate white noise. See ?arima.')
        }
      }
      if(length(parms)>0)
      {
        for(p in 1:length(parms))
        {
          if( !is.numeric(parms[[p]]) | !is.vector(parms[[p]]) )
          {
            stop('`', names(parms)[p], '` must be a numeric vector.')
          }
        }
      }
      # check family and famParms
      .checkFam(fam, famParms)

      # populate private
      private$.parms    <- parms
      private$.errVar   <- 1
      private$.FUN      <- arima.sim
      private$.fam      <- paste('r', fam, sep='')
      private$.famParms <- famParms

    },

    # TODO documentation for this method
    makeErrors = function(n, nObservations, seed=123)
    {
      # get seeds
      seeds <- as.list( .makeSeeds(seed, nObservations) )

      # sim errors
      errors <- .doLapply(seeds, self$fam, n=n,
               mu = self$famParms$mu, sigma = self$famParms$sigma,
               nu = self$famParms$nu, tau = self$famParms$tau)

      # TODO move this to a validation test or examples for this function.
      # rescale to have unit variances within time points then multiply
      # by the proportion of error variance; need to transpose twice so
      # that the rescaling applies within person (which also gets the
      # between person scaling right, but not vice versa in initial tests)
      #errorsr <- t(scale(t(errors), FALSE, TRUE)) * sqrt(design$variances$errorVar)
      #
      #doQC <- FALSE
      #if(doQC)
      #{
      #  # QC, only holds asymptotically
      #  # within persons:
      #  all.equal(mean(apply(errorsr, 1, var)), errorVar, tolerance = .05)
      #  # between persons:
      #  all.equal(mean(apply(errorsr, 2, var)), errorVar, tolerance = .05)
      #  tseries <- lapply(data.frame(t(errorsr)), ts)
      #  aFun <- function(x, order=c(1,0,0)) try(arima(x, order), silent = TRUE)
      #  sigma2s <- lapply(lapply(tseries, aFun), function(x) x$sigma2)
      #  # tends toward underestimation when nObservations are small, but very precise
      #  # when large
      #  hist(unlist(sigma2s))
      #}

      # return
      return(errors)

    }


  )
)


# note that in ?arima.sim, the rand.gen option defualts to rnorm, so if we follow
# Ty's model of E = a + e, we simply have a mixture of normal distributions, one
# that is autocorrelated, one that is not. So I don't really need a second
# child of err called measureErr, but the error structures are extensible
#n <- 1000
#hist(ya <- arima.sim(list(), n)) # white noise
#hist(yn <- rnorm(n))




