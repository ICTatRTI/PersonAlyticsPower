# function to populate error generators

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

armaErr <- R6::R6Class("errARMA",

  inherit = Err,

  public = list(
    initialize = function(
      parms    = list(ar=c(.5), ma=c(0)) ,
      errVar   = 1                       ,
      fam      = "qNO"                   ,
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
          if( !is.numeric(parms[[i]]) | !is.vector(parms[[i]]) )
          {
            stop('`', names(parms)[i], '` must be a numeric vector.')
          }
        }
      }


      # populate private
      private$.parms    <- parms
      private$.errVar   <- errVar
      private$.FUN      <- arima.sim
      private$.fam      <- fam
      private$.famParms <- famParms

    },

    makeErrors = function()
    {
      # move contents of ICTerror() here, currently in dataMakerUtils.R
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




