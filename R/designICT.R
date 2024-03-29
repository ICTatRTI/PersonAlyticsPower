#' active - active bindings for designICT R6 class definitions and its children
#'
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
.active <- function()
{
  list(

    # fields that can be updated by the user
    inputMat = function(value)
    {
      if( missing(value) ){ private$.inputMat }
      else
      {
        # make factors if not already factors
        value$Phase <- factor(value$Phase)
        value$Group <- factor(value$Group)

        # update phase objects
        phases <- makePhase(value$nObs[value$Group==levels(value$Group)[1]],
                            levels(value$Phase))
        private$.phases <- phases
        private$.phaseNames <- names(phases)
        for(i in seq_along(levels(value$Group)))
        {
          if(i>1)
          {
            value$nObs[value$Group==levels(value$Group)[i]] <-
              value$nObs[value$Group==levels(value$Group)[1]]
          }
        }
        private$.nObs <- sum(value$nObs[value$Group==levels(value$Group)[1]])

        # update group objects
        private$.groups <- value$n[value$Phase==levels(value$Phase)[1]]
        private$.groupNames <- levels(value$Group)
        for(i in seq_along(levels(value$Phase)))
        {
          if(i>1)
          {
            value$n[value$Phase==levels(value$Phase)[i]] <-
              value$n[value$Phase==levels(value$Phase)[1]]
          }
        }
        private$.n <- sum(value$n[value$Phase==levels(value$Phase)[1]])

        # update designMat
        maxRandFx <- length(names(value)[grepl("Mean", names(value))])
        private$.designMat <- makeDesignMat(phases, names(phases), maxRandFx)

        # warning
        warning("\n`inputMat` may have been updated to ensure consistent phase ",
                "\nlengths across groups and consistent group sizes across ",
                "\nphases. Check the input matrix is as intented:\n\n")
        print(value)
        cat("\n\n")

        # update cleaned up inputMat
        private$.inputMat <- value

        self
      }

    },

    randFxVar = function(value)
    {
      if( missing(value) ){ private$.randFxVar }
      else
      {
        if(any(value) <= 0)
        {
          stop("\nRandom effect variances must be greater than 0.")
        }

        if(length(value)==3)
        {
          message("\nAll values in `$inputMat` will be updated.",
                  "\nIf you want group and/or phase specific variances edit",
                  "\nthe `Var` columns in `$inputMat` directly.")
        }

        nms <- names(private$.inputMat)[which(grepl("Var", names(private$.inputMat)))]
        private$.inputMat[,nms] <- matrix(value, nrow(private$.inputMat), length(value),
                                          byrow = TRUE)

        private$.randFxVar <- value
        self
      }
    },

    randFxCor = function(value)
    {
      if( missing(value) ){ private$.randFxCor }
      else
      {
        stop("\n`$randFxCorMat` and `#randFxCor` will not be updated. Update",
             "\n`$randFxCorMat` directly instead.")
        private$.randFxCor <- value
        self
      }
    },

    randFxCorMat = function(value)
    {
      if( missing(value) ){ private$.randFxCorMat }
      else
      {
        private$.randFxCorMat <- value
        self
      }
    },

    propErrVar = function(value)
    {
      if( missing(value) ){ private$.propErrVar }
      else
      {
        if(length(value)!=3 | sum(value) != 1)
        {
          stop("\n`propErrVar` must be a length 3 vector that sums to 1.")
        }

        nms <- c('randFx', 'res', 'mserr')

        if(length(value)==3)
        {
          message("\nAll values in `$inputMat` will be updated with\n",
                  paste(paste(nms, value, sep = "="), collapse=", "),
                  "\nIf you want group and/or phase specific values edit",
                  "\n`$inputMat` directly.")
        }

        if(names(value)!=nms) names(value) <- nms

        private$.inputMat[,nms] <- matrix(value, nrow(private$.inputMat), length(value),
                                          byrow = TRUE)

        names(value) <- nms
        private$.propErrVar <- value
        self
      }
    },

    error = function(value)
    {
      if( missing(value) ){ private$.error }
      else
      {
        if(! "err" %in% class(value))
        {
          stop("\nThe object supplied is not an `err` object. See `?armaErr`.")
        }
        private$.error <- value
        self
      }
    },

    merror = function(value)
    {
      if( missing(value) ){ private$.merror }
      else
      {
        if(! "err" %in% class(value))
        {
          stop("\nThe object supplied is not an `err` object. See `?armaErr`.")
        }
        private$.merror <- value
        self
      }
    },

    yMean = function(value)
    {
      if( missing(value) ){ private$.yMean }
      else
      {
        private$.yMean <- value
        self
      }
    },

    ySD  = function(value)
    {
      if( missing(value) ){ private$.ySD }
      else
      {
        private$.ySD <- value
        self
      }
    },

    yMin  = function(value)
    {
      if( missing(value) ){ private$.yMin }
      else
      {
        private$.yMin <- value
        self
      }
    },

    yMax  = function(value)
    {
      if( missing(value) ){ private$.yMax }
      else
      {
        private$.yMax <- value
        self
      }
    },

    yCut  = function(value)
    {
      if( missing(value) ){ private$.yCut }
      else
      {
        if( sum(yCut) != 1 )
        {
          stop("\n`yCut` mut be a vector of proportions summing to 1.")
        }
        private$.yCut <- value
        self
      }
    },


    # fields than cannot be updated by the user
    randFxOrder = function(value)
    {
      if( missing(value) ){ private$.randFxOrder }
      else
      {
        stop("\n`randFxOrder` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    meanNames = function(value)
    {
      if( missing(value) ){ private$.meanNames }
      else
      {
        stop("\n`meanNames` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    varNames = function(value)
    {
      if( missing(value) ){ private$.varNames }
      else
      {
        stop("\n`varNames` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    groups = function(value)
    {
      if( missing(value) ){ private$.groups }
      else
      {
        stop("\n`groups` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    n = function(value)
    {
      if( missing(value) ){ private$.n }
      else
      {
        stop("\n`n` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    phases= function(value)
    {
      if( missing(value) ){ private$.phases }
      else
      {
        stop("\n`phases` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    nObs = function(value)
    {
      if( missing(value) ){ private$.nObs }
      else
      {
        stop("\n`nObs` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    designMat = function(value)
    {
      if( missing(value) ){ private$.designMat }
      else
      {
        stop("\n`designMat` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    phaseNames = function(value)
    {
      if( missing(value) ){ private$.phaseNames }
      else
      {
        stop("\n`phaseNames` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    groupNames = function(value)
    {
      if( missing(value) ){ private$.groupNames }
      else
      {
        stop("\n`groupNames` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    randFxFam = function(value)
    {
      if( missing(value) ){ private$.randFxFam }
      else
      {
        stop("\n`randFxFam` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    },

    randFxFamParms = function(value)
    {
      if( missing(value) ){ private$.randFxFamParms }
      else
      {
        stop("\n`randFxFamParms` is read only, make changes to `inputMat` instead.\n\n",
             call. = FALSE)
      }
    }
  )
}


#' \code{designICT} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @export
#'
#' @description
#' The \code{designICT} class is a parent class to \code{\link{polyICT}} and
#' other future planned classes such as \code{expICT} for exponential growth
#' models and \code{logisticICT} for logistic ICT growth models. The
#' \code{designICT} class contains methods (functions common to all child
#' classes) and fields (variables common to all child classes) and is not
#' to be used directly.
#'
#' @section Methods:
#' \describe{
#'
#'   \item{\code{print}}{
#'
#'     A method for printing a the primary inputs for an ICT design.
#'
#'   }
#'
#'   \item{\code{designCheck}}{
#'
#'     This method simulates one dataset using the current ICT design and
#'     plots average trajectories and histograms of the outcome by phase and
#'     by group. This allows the user to visually check whether the inputs they
#'     provided match their intendend design. The default is to do this with
#'     large samples (N=5,000 per group) to show the asypmtotic behavior. For
#'     comparison, smaller sample sizes can also be used to see the possible
#'     deviations from the population. If smaller sample sizes are used, this
#'     should be repeated several times to prevent over-interpretation of one
#'     sample. See the examples in \code{\link{polyICT}}.
#'
#'     \code{file} Character. A file name for saving the resulting data. The
#'     default is NULL, in which case the data are not saved.
#'
#'     \code{ylim} Numeric vector of length 2. Limits for the y-axis. This is
#'     useful when variance is high and/or effect sizes are small, making it
#'     difficult to visualize average change over time. The default is NULL, in
#'     which case default limits are used.
#'
#'     \code{fitMod} Logical. Should the model implied by the design be fit to
#'     the simulated data? The default is FALSE.
#'
#'     \code{seed} Numeric. A random seed for enabling replication. The default
#'     is 123.
#'
#'     \code{npg} Numeric (integer). The number per group. The default is
#'     n=5,000 per group. If a small value is used, the user should repeat the
#'     check with several different seeds. See the examples in
#'     \code{\link{polyICT}}.
#'
#'     \code{return} Character. What should  be returned. The default is the
#'     designCheck plot via \code{return='plot'}. Other options are \code{'model'}
#'     which returns model fit and descriptive statistics; \code{'stats'} which
#'     returns descriptive statistics only; \code{'data'} which returns the
#'     simulated data; \code{'datastat'} which returns the data and the
#'     descriptive statistics; and \code{'all'} which returns everything in a list.
#'
#'     \code{justData} Logical. Should only the data be generated (and not
#'     a model or plot). This speeds up other calls to \code{designCheck}, e.g.,
#'     by \code{\link{setupDist}}.
#'   }
#'
#' }
#'
#'


designICT <- R6::R6Class("designICT",

                         # this list must contain all child slots, children cannot update
                         # private (to my knowledge)
     private = list(

       # editable without a new call to $new
       .inputMat          = NULL,
       .randFxVar         = NULL,
       .randFxCor         = NULL,
       .randFxCorMat      = NULL,
       .propErrVar        = NULL,
       .error             = NULL,
       .merror            = NULL,
       .yMean             = NULL,
       .ySD               = NULL,
       .yMin              = NULL,
       .yMax              = NULL,
       .yCut              = NULL,

       # not editable
       .n                 = NULL,
       .nObs              = NULL,
       .groups            = NULL,
       .phases            = NULL,
       .designMat         = NULL,
       .phaseNames        = NULL,
       .groupNames        = NULL,
       .randFxOrder       = NULL,
       .meanNames         = NULL,
       .varNames          = NULL,

       # not implemented
       .randFxFam         = NULL,
       .randFxFamParms    = NULL
     ),

     # Child classes of designICT will inherit active bindings and
     # all active bindings for children will stay in active()
     active = .active(),

     public = list(

       # initialize is child-specific, so designICT does not have an
       # initialize function

       print = function(...)
       {
         phases <- paste(lapply(self$phases, function(x) shQuote(x[1])),
                         'for',
                         lapply(self$phases, length), 'time points\n')
         phases  <- paste(c('', rep(' ', length(phases)-1)), phases)
         nGroups <- length(self$groupNames)
         if(nGroups==0) nGroups <- 1
         cat(.hl(),
             "Design Information: An ICT with ", sum(self$n), " participants, ",
             length(c(unlist(self$phases))), " time points, and ",
             ifelse(!is.null(self$groupNames), nGroups, 1),
             ifelse(nGroups==1, " group.\n", " groups.\n"), .hl(),
             "\nPhases:\n",
             paste(phases, collapse = ''),
             "\nGroups:\n", paste(self$groupNames, 'n=', self$groups, '\n\n'),
             .hl(),
             sep=""
         )
         cat("Input Matrix:\n")
         cat(.hl())
         print( self$inputMat )
         cat(.hl())
         cat("\nThe input matrix can be edited by typing `edit(x$inputmat)`\n",
             "where `x` is the name of your ICT design object.\n",
             "Phase : the labels of the study phases.\n",
             "Group : the names of the groups.\n",
             "nObs  : the number of observations per phase.\n",
             "n     : the number of participants per group.\n",
             "Mean# : the effect size (Cohen's d) per phase & group, where\n",
             "        `#` is the random effect order (0=intercepts, 1=slopes, etc.)\n",
             "Var#  : the variance per phase & group, where",
             "`#` is as explained above.\n",
             "\nThe proportion of the total variance at time 1:\n",
             "randFx: due to random effects.\n",
             "res   : due to residual autocorrelation.\n",
             "mserr : due to measurement error.\n")
         invisible(self)
       },

       designCheck = function(file=NULL, ylim=NULL, fitMod=FALSE,
                              seed=123, npg=2000, return = 'plot',
                              justData=FALSE, type='histogram',
                              title="Design Check")
       {
         # message
         message("\nRunning a large sample `designCheck` with ", npg,
                 " participants per group.\n",
                 "This may take a few minutes to run.\n\n")

         # set mod0 default
         mod0 <- paste("`fitMod` was set to", as.character(fitMod))

         # save and reset n, due to inheritance design will get overwritten, fix
         # n below
         originaln <- aggregate(self$inputMat$n, list(self$inputMat$Group),
                                function(x) x[1])
         tempn <- rep(npg, length(originaln))
         names(tempn) <- names(self$groups)
         for(i in seq_along(tempn))
         {
           self$inputMat$n[self$inputMat$Group==self$groupNames[i]] <- tempn[i]
         }

         # simulate one big data set using `npg` and print descriptive stats
         dat <- self$makeData(seed=seed)
         descriptives <- dstats(dat$y)

         # needs to be reimplemented
         if(1==2)
         {
           # get expected means and variances
           expectedVar <- self$expectedVariances

           # compare expected to observed variances
           expObsVar   <- cbind( expectedVar$Variances,
                                 aggregate(dat$y, by=list(dat$Time), var)$x)
           # get the correlation of expected and observed variances
           expObsVcor <- round(cor(expObsVar)[1,2],4)
           cat("The correlation between the expected variances and the observed\n",
               "variances are ", expObsVcor, "\n\n")

           # compare expected to observed means
           expObsMean <- cbind( expectedVar$Means,
                                aggregate(dat$y, by=list(dat$Time), mean)$x)
           expObsMcor <- round(cor(expObsMean)[1,2],4)
           cat("The correlation between the expected means and the observed\n",
               "means are ", expObsMcor, "\n\n")

           # compare expected to observed total
           TotalMean <- mean(dat$y)
           TotalVar  <- var(dat$y)
           cat("The observed total mean is ", TotalMean,
               "\nThe expected total mean is ", expectedVar$TotalMean,
               "\nThe observed total variance is ", TotalVar,
               "\nThe expected total variance is ", expectedVar$TotalVar, "\n\n")

           #TODO this only works for linear models
           #TODO needs better matching, force name matching in polyICT
           Estimates <- round(rbind(summary(mod0)$tTable[,1]),3 )
           Inputs    <- round(c(unlist(self$InputMat))[c(1,3,2,4)],3)
           cat('\n\nCheck the effect size estimates against inputs:\n')
           print( data.frame(Inputs=Inputs, Estimates=t(Estimates)) )

           cat("\n\n\n")
         }

         if(!justData)
         {
           # set up the Palytic object and, if requested, fit the model
           correlation <- paste("corARMA(p=", length(self$error$parms$ar), ", ",
                                "q=", length(self$error$parms$ma), ")", sep="")
           pa   <- Palytic$new(data=dat, ids='id', dv='y', time='Time',
                               phase=unlist(ifelse(length(self$phaseNames)>1,
                                                   'phase',list(NULL))),
                               ivs=unlist(ifelse(length(self$groupNames)>1,
                                                 'group',list(NULL))),
                               time_power = self$randFxOrder,
                               correlation = correlation)

           if(fitMod) # runs slow with some examples, qc why
           {
             mod0 <<- pa$lme() # need to add option for gamlss
             print( summary( mod0 ) )
             # save data if requested
             if(!is.null(file)) save(mod0, file=paste(file[1],
                                                      'designCheck.RData',
                                                      sep='_'))
           }

           # autoscale ylim, needs assymetric update for skewed data
           ylim <- c(min(dat$y, na.rm=T), max(dat$y, na.rm=T))
           if(is.null(self$ylim) & !is.null(self$ySD) & is.null(self$yCut))
           {
             ylim <- c(-2.75*self$ySD, 2.75*self$ySD) + self$yMean
           }

           # plot
           if( length( self$groupNames ) == 1 ) gg <- pa$plot(ylim=ylim,
                                                              type=type,
                                                              title=title)
           if( length( self$groupNames ) >= 2 ) gg <- pa$plot(groupvar = 'group',
                                                              ylim=ylim,
                                                              type=type,
                                                              title=title)
         }

         # restore the original sample sizes
         for(i in seq_along(tempn))
         {
           self$inputMat$n[self$inputMat$Group==self$groupNames[i]] <- originaln$x[i]
         }

         # return
         if(return=='plot'    ) return( gg )
         if(return=='model'   ) return( list(model = mod0,
                                             descriptives = descriptives) )
         if(return=='stats'   ) return( descriptives )
         if(return=='data'    ) return( dat )
         if(return=='datstat' ) return( list(data=dat, descriptives = descriptives) )
         if(return=='all'     ) return( list(plot=gg, model=mod0, data=dat,
                                    descriptives = descriptives) )
       }

   )
) # end of designICT class definition




# use inheritance to create specific instances of designICT
#
# current:
# - polyICT: polynomial ICT
# - polyICT2: polynomial ICT via distribution subsetting
#
# Planned (see Singer & Willett page 234, Table 6.7):
# - polyICT
# - hyperICT
# - invPolyICT
# - expICT
# - negExpICT
# - logisticICT





