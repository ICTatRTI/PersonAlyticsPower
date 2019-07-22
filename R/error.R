
#' .arima.sim - function to run arima.sim with potentially extraneous arguments
#' intended for a gamlss.family distribution passed to rand.gen
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
.arima.sim <- function(seed, model, n, rand.gen = rnorm,
                       mu=NULL, sigma=NULL, nu=NULL, tau=NULL,
                       doPlot=TRUE, doStats=TRUE, debug=FALSE)
{
  set.seed(seed)

  if(debug) cat('mu=',mu,'; sigma=', sigma, '; nu=', nu, '; tau=', tau, '\n\n')

  # R.utils::doCall() isn't working, so we'll be explicit
  if( !is.null(mu) & is.null(sigma) & is.null(nu) & is.null(tau) )
  {
    if(debug) cat('mu\n\n')
    testData <-arima.sim( model=model, n = n,
                          rand.gen = rand.gen,
                          mu = mu )
  }
  if( !is.null(mu) & !is.null(sigma) & is.null(nu) & is.null(tau) )
  {
    if(debug) cat('mu, sigma\n\n')
    testData <-arima.sim( model=model, n = n,
                          rand.gen = rand.gen,
                          mu = mu, sigma = sigma )
  }
  if( !is.null(mu) & !is.null(sigma) & !is.null(nu) & is.null(tau) )
  {
    if(debug) cat('mu, sigma, nu\n\n')
    testData <-arima.sim( model=model, n = n,
                          rand.gen = rand.gen,
                          mu = mu, sigma = sigma,
                          nu = nu )
  }
  if( !is.null(mu) & !is.null(sigma) & is.null(nu) & !is.null(tau) )
  {
    if(debug) cat('mu, sigm, nu, tau\n\n')
    testData <-arima.sim( model=model, n = n,
                          rand.gen = rand.gen,
                          mu = mu, sigma = sigma,
                          nu = nu, tau = tau )
  }
  if( is.null(mu) & is.null(sigma) & is.null(nu) & is.null(tau) )
  {
    if(debug) cat('rnorm')
    testData <-arima.sim( model=model, n = n,
                          rand.gen = rand.gen, ... )
  }

  if(doPlot)
  {
    par(mfrow=c(2,2))
    hist( testData )
    qqnorm( testData ); qqline(testData, col='red')
    plot( testData )
    #acf( testData )
    pacf( testData )
    par(mfrow=c(1,1))
  }

  if(doStats) dstats(testData)

  invisible( testData )
}

#' errActive Active bindings for Err and its heirs
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
errActive <- function()
{
  list(

    model = function(value)
    {
      if( missing(value) ){ private$.model }
      else
      {
        private$.model <- value
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
#' @field model The parameters for \code{FUN}.
#'
#' @field errVar The error variance.
#'
#' @field FUN The function for simulating error terms.
#'
#' @field fam A random variable generator from the
#' \code{\link{gamlss.family}} for the error term distribution.
#' After data are simulated using \code{FUN} and its \code{model}, it can be
#' transformed to have a distribution implied by \code{fam}.
#'
#' @field famParms Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.
Err <- R6::R6Class("err",

         private = list(
           .model    = NULL,
           .errVar   = NULL,
           .FUN      = NULL,
           .fam      = NULL,
           .famParms = NULL
         ),

         active = errActive(),

         public = list(

           initialize = function(
              model    = NULL,
              errVar   = 1   ,
              FUN      = NULL,
              fam      = NULL,
              famParms = NULL
           )
           {
             # validate inputs

             # populate private
             private$.model    <- model
             private$.errVar   <- errVar
             private$.FUN      <- FUN
             private$.fam      <- fam
             private$.famParms <- famParms
           },

           print = function(...)
           {
             print(self$model)
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
#' @field model The parameters for \code{FUN}.
#'
#' @field errVar The error variance.
#'
#' @field FUN The function for simulating error terms.
#'
#' @field fam A \code{\link{gamlss.family}} family for the error term distribution.
#' After data are simulated using \code{FUN} and its \code{model}, it can be
#' transformed to have a distribution different from those available to \code{FUN}.
#'
#' @field famParms Parameters to be passed to \code{fam}, see \code{\link{gamlss.family}}.
#'
#' @examples
#'
#' # set up a stationary arma model
#' testErr <- armaErr$new(model = list(ar=c(.5), ma=c(.2)))
#' testErr$checkModel()
#'
#' # show testing for an unstationary ARMA model
#' testFail <- armaErr$new(model = list(ar=c(-.8, .5), ma=c(.2)))
#' testFail$checkModel()
#'
#' # show that r* functions from gamlss.family distribution functions can be
#' # passed to arima.sim
#' errBeta <- arima.sim(list(ar=c(.5)), 1000, innov = rBE(1000, mu=.5, sigma=.2))
#'
#' # note that the rand.gen option only allows the defaults for the function
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
      model    = list(ar=c(.5), ma=c(0)) ,
      fam      = "NO"                    ,
      famParms = list(mu=0, sigma=1)
    )
    {
      # validation
      if( !is.list(model) )
      {
        if( ! all(names(model) %in% c('ar', 'ma')) )
        {
          stop('\n`model` must be a named list of length one or two with the\n',
             'names `ar` and/or `ma`; `model` can also be an empty list\n',
             'to generate white noise. See ?arima.')
        }
      }
      if(length(model)>0)
      {
        for(p in 1:length(model))
        {
          if( !is.numeric(model[[p]]) | !is.vector(model[[p]]) )
          {
            stop('`', names(model)[p], '` must be a numeric vector.')
          }
        }
      }

      # check family and famParms
      .checkFam(fam, famParms)

      # define random generator for fam
      fam <- paste('r', fam, sep='')

      # populate private
      private$.model    <- model
      private$.errVar   <- 1
      private$.FUN      <- arima.sim
      private$.fam      <- fam
      private$.famParms <- famParms

    },

    checkModel = function(seed=1234, n = 1000, doPlot=TRUE, doStats=TRUE)
    {
      # the following will yield an error if the ar and/or ma parameters yield a
      # model that is not stationary

      # note that we use eval parse without security concerns as .checkFam has
      # already validated the input
      check <- try(
      .arima.sim(seed=seed, doPlot=doPlot, doStats = doStats,
                 model=self$model, n = n,
                 rand.gen = eval(parse(text=self$fam)),
                 mu = self$famParms$mu, sigma = self$famParms$sigma,
                 nu = self$famParms$nu, tau = self$famParms$tau))
      if( is(check)[1] == 'try-error' )
      {
        cat(0, file='stationary.arma') # for chris's gui
      }
      else invisible(check)
    },

    makeErrors = function(n, nObservations, seed=123)
    {
      # get seeds
      seeds <- as.list( .makeSeeds(seed, n) )

      # sim errors
      errors <- lapply(seeds, .arima.sim, doPlot=FALSE, doStats=FALSE,
                       model=self$model, n = nObservations,
                       rand.gen =  eval(parse(text=self$fam)),
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
      return( t( do.call(cbind, errors) ) )

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




