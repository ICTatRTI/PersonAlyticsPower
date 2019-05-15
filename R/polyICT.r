
#' \code{polyICT} class generator
#'
#' @docType class
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @import MASS
#'
#' @export
#'
#' @keywords data
#'
#' @usage polyICT$new()
#'
#' @details
#' The \code{polyICT} class generator specifies the inputs needed to simulate
#' data from a polynomial growth model. Once a \code{polyICT} object is created,
#' you can called its methods and examine or update its fields.
#'
#' Methods are functions that come packaged within your \code{polyICT} object
#' and include \code{$print()} for printing the inputs,
#' \code{$update()} for changing the inputs,
#' \code{$designCheck()} for visualizing the design, and
#' \code{$makeData} for simulating a single data set. See the section
#' titled \strong{Methods}.
#'
#' Fields are all the data stored within a \code{polyICT} object,
#' some of which are provided by the user when initializing a \code{polyICT}
#' object, and others which are derived from these inputs and cannot be changed
#' by the user directly. These are detailed in the section titled \strong{Fields}.
#' Fields can be accessed using the \code{$} operator. For example, if your
#' \code{polyICT} object is called \code{myPolyICT}, use \code{myPolyICT$inputMat}.
#'
#' @field inputMat A \code{\link{data.frame}} containing the inputs needed for
#' data simulation by phase and by group. Columns include \code{Phase},
#' \code{Group}, \code{nObs} (the number of observations in a given phase),
#' \code{n} (the number of participants in a given group), the means and
#' variances for each random effect (see \code{randFxOrder} under \code{new} in
#' the Methods section), and the variance partioning (see
#' \code{propErrVar} under \code{new} in the Methods section). Instructions
#' for editing this field are given in the Examples.
#'
#' The \code{Mean} columns are standardized effect sizes on the scale of Cohen's
#' \emph{d}. For more detailed information and illustrations, see the Example
#' section.
#'
#' @field randFxVar See \code{randFxVar} in the \code{new} method. Phase and/or
#' group specific variances can be specified by editing \code{inputMat} after
#' initializing a \code{polyICT} object.
#'
#' @field randFxCorMat A list of correlation matrices for each phase and group.
#' The default \code{randFxCorMat}s are ceated using \code{randFxCor} for all
#' off-diagonal entries. These can be edited as illustrated in the Examples.
#'
#' @field propErrVar See \code{propErrVar} in the \code{new} Method. See also
#' \code{inputMat} and the Examples for making these inputs phase and/or group
#' specific.
#'
#' @field error See \code{error} in the \code{new} Method. See also
#' \code{\link{errARMA}}.
#'
#' @field merror See \code{merror} in the \code{new} Method.
#'
#' @field yMean See \code{yMean} in the \code{new} Method.
#'
#' @field ySD See \code{ySD} in the \code{new} Method.
#'
#' @field n The total sample size.
#'
#' @field nObs The total number of observations (i.e., time points).
#'
#' @field groups See \code{groups} in the \code{new} Method.
#'
#' @field phases See \code{phases} in the \code{new} Method.
#'
#' @field designMat The design matrix with phases and timepoints.
#'
#' @field meanNames The columns of \code{InputMat}
#' corresponding to the effect sizes/random effect means.
#'
#' @field varNames The columns of \code{InputMat}
#' corresponding to the random effect variances.
#'
#' @field phaseNames The names of the phases, taken from \code{phases}.
#'
#' @field groupNames The names of the groups, taken from \code{groups}.
#'
#' @field randFxFam A \code{\link{gamlss.family}} distribution for non-normal
#' random effects. Not implemented.
#'
#' @field randFxFamParms The parameters for the \code{\link{gamlss.family}}
#' distribution specified in \code{randFxFam}. Not implemented.
#'
#' @section Methods:
#' \describe{
#'
#'   \item{\code{new}}{Used to initialize a \code{polyICT} object as illustrated
#'   in the \strong{Examples} below. The following
#'   parameters can be passed to \code{$new()}:
#'
#'   \code{groups} Named numeric vector. The default is
#'   \code{c(group1=10, group2=10)}. The values are the number of participants
#'   per group and the names are the group names.
#'
#'   \code{phases} Named list. The default is created using the helper function
#'   \code{\link{makePhase}}. Each item in the list replicates the phase name
#'   as many times as there are time points in that phase. Actual time points
#'   are derived during initialization and stored in the field \code{designMat}.
#'
#'   \code{propErrVar} Named numeric vector of length 3. The default is
#'   \code{c(randFx=.5, res=.25, mserr=.25)}. The names must be
#'   \code{randFx}, the proportion of the total variance due to random effects;
#'   \code{res}, the proportion of the total variance due to residual
#'   autocorrelation; and \code{mserr}, the proportion of total variance due to
#'   measurement error. The three values must be proportions and must sum to 1.
#'
#'   \code{randFxOrder} Numeric vector. The default is \code{1}. This is used
#'   to specify the order of the polynomial growth model as follows:
#'   \code{randFxOrder=0} is an intercept only model, \code{randFxOrder=1} adds
#'   random slopes, \code{randFxOrder=2} is a quadratic growth model,
#'   \code{randFxOrder=3} is a cubic growth model, etc.
#'
#'   \code{randFxCor} Numeric. The correlation(s) between all of the random
#'   effects. This can be edited later and made group and/or phase specific. See
#'   the \strong{Examples}.
#'
#'   \code{randFxVar} Numeric vector. The default is \code{c(1, 1)}. This
#'   parameter is used to specify the variances of the random effects. The
#'   number of variances should equal \code{randFxOrder + 1}. The random effects
#'   variances are created as
#'   \code{Var[o] = propErrVar$randFx * (randFxVar[o]/sum(randFxVar))} where
#'   \code{o=0:randFxOredr}. Hence \code{randFxVar} specifies the ratios of
#'   the variance that is partitioned among the random effects and rescaled by
#'   \code{propErrVar$randFx}. See \code{makeData} in the \strong{Methods}
#'   section for more details.
#'
#'   \code{error} An error object for the residual autocorrelation. The default
#'   is \code{armaErr$new(list(ar=c(.5), ma=c(0)))}, a first order AR process
#'   with \eqn{phi_1=.5}. See \code{\link{armaErr}}. See \code{makeData} in the
#'   \strong{Methods} section for more details.
#'
#'   \code{merror} An error object for the measurement error. The default
#'   is \code{armaErr$new(list()}, a white noise process. See
#'   \code{\link{armaErr}}. Also see \code{makeData} in the \strong{Methods}
#'   section for more details.
#'
#'   \code{ySD} Numeric. The default is \code{15}. The standard deviation of the
#'   final data. See \code{makeData} in the \strong{Methods}
#'   section for more details.
#'
#'   \code{ySD} Numeric. The default is \code{100}. The mean of the
#'   final data. See \code{makeData} in the \strong{Methods}
#'   section for more details.
#'
#'   }
#'
#'   \item{\code{print}}{See \code{\link{designICT}}.}
#'
#'   \item{\code{designCheck}}{See \code{\link{designICT}}.}
#'
#'   \item{\code{update}}{A method to update editable field in a \code{polyICT}
#'   object. Fields that can be updated are those listed in \code{new}. New
#'   values are passed by name using, for example,
#'   \code{$update(groups=c(group1=25, group2=25), randFxOrder=2)}. The are no
#'   defaults and any number of fields can be updated at once.
#'   }
#'
#'   \item{\code{makeData}}{A method to simulate one data set from the settings
#'   in a given \code{polyICT} object. This method is not intended to be used
#'   directly by the user, who should instead use \code{\link{ICTpower}} to
#'   automate a power analysis for one ICT design, or \code{\link{ICTpowerSim}} for
#'   conducting a full power analysis simulation study. The parameters are
#'
#'     \code{seed} Numeric. The default is \code{123}. A random seed for
#'     reproducibility. If multiple calls are made to \code{makeData}, the seed
#'     should change for each call as is done automatically by
#'     \code{\link{ICTpower}}.
#'   }
#' }
#'
#'
#' @examples
#' # Set up a new polyICT object using the default parameter settings
#'
#' myPolyICT <- polyICT$new(
#'   groups            = c(group1=10, group2=10)          ,
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
#' # print the object
#'
#' myPolyICT
#' # equivalent to:
#' #print(myPolyICT)
#' #myPolyICT$print()
#'
#' # fields can be accessed directly using $
#' myPolyICT$inputMat
#' myPolyICT$designMat
#'
#' # edit the means in inputMat so that
#' #
#' # 1. group1 is left with no change in all three phases
#' # 2. group2 has no change at phase1 (i.e., baseline), phase2 has a phase
#' #    jump to d=.3 with no within phase change, and phase3 starts where it
#' #    left off at d=.3 and decreases by d=-.6 through the remained of the
#' #    phase.
#' #
#' # Note that editing inputMat may be easier using
#' #edit(myPolyICT)
#' # but this is more diffucult to replicate.
#' myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase2' &
#'   myPolyICT$inputMat$Group=='group2', 'Mean0'] <- .3
#' myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
#'   myPolyICT$inputMat$Group=='group2', 'Mean0'] <- .3
#' myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
#'   myPolyICT$inputMat$Group=='group2', 'Mean1'] <- -.6
#' myPolyICT$inputMat
#'
#' # now do a large sample (n=5,000/group) check of the design using $designCheck.
#' # Notice that in group 2, there is a jump of about .3*15=4.5 between phases
#' # one and two, and that within phase 3 there is -.6*15=-9 point reduction. This
#' # is not run when calling example(polyICT) due to it taking several seconds to
#' # simulate the large sample data.
#'
#'
#' \dontrun{
#'
#' myPolyICT$designCheck(ylim=c(75,125))
#'
#' # for comparison, see what can happen under finite samples (n=10/group).
#' # Notice that there are substantial changes in group1 even though none are
#' # specified, and that the changes in group2 can be different from what was
#' # specified. This illustrates possible. finite sample behaviors.
#'
#' myPolyICT$designCheck(seed=1, npg=10, ylim=c(75,125))
#' myPolyICT$designCheck(seed=2, npg=10, ylim=c(75,125))
#' myPolyICT$designCheck(seed=3, npg=10, ylim=c(75,125))
#'
#' }

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# start of polyICT class ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
polyICT <- R6::R6Class("polyICT",

                       inherit = designICT,

                       public  = list(
                         initialize = function
                         (
                           groups            = c(group1=10, group2=10)                   ,
                           phases            = makePhase()                               ,
                           propErrVar        = c(randFx=.5,res=.25,mserr=.25)            ,
                           randFxOrder       = 1                                         ,
                           randFxCor         = 0.2                                       ,
                           randFxVar         = c(1, 1)                                   ,

                           error             = armaErr$new()                             ,
                           merror            = armaErr$new(list())                       ,
                           ySD               = 15                                        ,
                           yMean             = 100                                       ,

                           # these should be hidden from users until we can model
                           # non-normal random effects
                           randFxFam         = "qNO"                                     ,
                           randFxFamParms    = list(mu=.5 , sigma=1)

                         )
                         {
                           # makeDesign
                           design <- makeDesign(0:randFxOrder, phases, groups,
                                                propErrVar, randFxVar,
                                                randFxCor, design = 'polyICT')


                           #
                           # populate private
                           #
                           #private$.edit              <- TRUE

                           # editable without a new call to $new or via $update
                           private$.inputMat          <- design$inputMat
                           private$.randFxVar         <- design$randFxVar
                           private$.randFxCor         <- randFxCor
                           private$.randFxCorMat      <- design$randFxCorMat
                           private$.propErrVar        <- design$propErrVar
                           private$.error             <- error
                           private$.merror            <- merror
                           private$.yMean             <- yMean
                           private$.ySD               <- ySD

                           # not editable
                           private$.n                 <- design$n
                           private$.nObs              <- design$nObs
                           private$.groups            <- design$groups
                           private$.phases            <- design$phases
                           private$.designMat         <- design$designMat
                           private$.meanNames         <- design$meanNames
                           private$.varNames          <- design$varNames
                           private$.phaseNames        <- names(phases)
                           private$.groupNames        <- names(groups)
                           private$.randFxOrder       <- randFxOrder

                           # not implemented
                           private$.randFxFam         <- randFxFam
                           private$.randFxFamParms    <- randFxFamParms

                           # variances
                           # TODO repopulate these after moving their function to ??
                           #private$.variances         <- variances
                           #private$.expectedVariances <- expectedVariances
                           # internals
                           #private$.edit              <- FALSE


                         },

                         print = function(...)
                         {
                           # use super to call the print function of the parent
                           # class, then you can add anything specific to this
                           # class
                           super$print()
                           invisible(self)
                         },

                         # update function to repeate QC steps in $new (which
                         # calls $initialize) rather than repeating this code
                         # in each active binding
                         update = function(...)
                         {
                           # Turn off warnings so that active binding messages,
                           # which should be implemented as warnings, are
                           # turned off (checkPolyICT will handle via errors)
                           options(warn=-1)
                           #self$edit <- TRUE

                           # This will call each active binding in self for
                           # passed parameters, hence warn=-1
                           dots <-  list(...)#match.call(expand.dots = FALSE)$...
                           #print(dots)
                           dotsNames <- names(dots)
                           for(i in seq_along(dots))
                           {
                             if( dotsNames[i] %in% names(self) )
                             {
                               if( dotsNames[i] == 'randFxOrder' )
                               {
                                 warning('Updating `randFxOrder` will overwrite',
                                         '`inputMat`.')
                                 ans <- readline(prompt =
                                          'Do you want to do this? y/n: ')
                                 if(ans=='n') stop('Call to $update() canceled.')
                               }
                               if( dotsNames[i] == 'phases' )
                               {
                                 pNew <- length(dots[[i]])
                                 if(pNew != length(self$phases))
                                 {
                                   warning('Updating the number of `phases` ',
                                           'will overwrite `inputMat`.')
                                   ans <- readline(prompt =
                                            'Do you want to do this? y/n: ')
                                   if(ans=='n') stop('Call to $update() canceled.')
                                 }
                               }
                               if( dotsNames[i] == 'groups' )
                               {
                                 gNew <- length(dots[[i]])
                                 if(gNew != length(self$groups))
                                 {
                                   warning('Updating the number of `groups` ',
                                           'will overwrite `inputMat`.')
                                   ans <- readline(prompt =
                                                     'Do you want to do this? y/n: ')
                                   if(ans=='n') stop('Call to $update() canceled.')
                                 }
                               }
                               if( dotsNames[i] == 'propErrVar' )
                               {
                                 warning('Updating `propErrVar` will apply the',
                                         ' new values to all groups and phases\n')
                                 ans <- readline(prompt =
                                                   'Do you want to do this? y/n: ')
                                 if(ans=='n') stop('Call to $update() canceled.')
                                 if(ans=='y')
                                 {
                                   message('Use `edit(myPolyICT$inputMat)` for\n',
                                           ' group and/or phase specific values',
                                           ' of `propErrVar`.')
                                 }
                               }
                               if( dotsNames[i] == 'randFxVar' )
                               {
                                 warning('Updating `randFxVar` will apply the',
                                         ' new values to all groups and phases\n')
                                 ans <- readline(prompt =
                                                   'Do you want to do this? y/n: ')
                                 if(ans=='n') stop('Call to $update() canceled.')
                                 if(ans=='y')
                                 {
                                   message('Use `edit(myPolyICT$inputMat)` for\n',
                                           ' group and/or phase specific values',
                                           ' of `randFxVar`.')
                                 }
                               }
                               if( dotsNames[i] == 'randFxCor' )
                               {
                                 warning('Updating `randFxCor` will apply the',
                                         ' new values to all groups and phases\n')
                                 ans <- readline(prompt =
                                                   'Do you want to do this? y/n: ')
                                 if(ans=='n') stop('Call to $update() canceled.')
                                 if(ans=='y')
                                 {
                                   message('Use `edit(myPolyICT$randFxCorMat)`\n',
                                           ' or see ?polyICT for',
                                           ' group and/or\nphase specific values',
                                           ' of `randFxCorMat`.')
                                 }
                               }
                               self[[dotsNames[i]]] <- dots[[i]]
                             }
                             if( ! dotsNames[i] %in% names(self) )
                             {
                               message("The parameter `", dotsNames[i], "` is not ",
                                    "a valid parameter for a `polyICT`\n",
                                    "object and will be ignored.")
                             }
                           }
                           # only run this if
                           if(names(dots) %in% c('randFxOrder', 'phases',
                                                 'groups', 'propErrVar',
                                                 'randFxCor', 'randFxVar'))
                           {
                             design <- makeDesign(
                               randFxOrder  = 0:self$randFxOrder ,
                               phases       = self$phases        ,
                               groups       = self$groups        ,
                               propErrVar   = self$propErrVar    ,
                               randFxVar    = self$randFxVar     ,
                               randFxCor    = self$randFxCor     ,
                               design       = 'polyICT'          ,
                               makeInputMat = FALSE              ,
                               self         = self               ,
                               isNew        = dotsNames          )

                             self$inputMat      <- design$inputMat
                             self$randFxVar     <- design$randFxVar
                             self$randFxCorMat  <- design$randFxCorMat
                             self$propErrVar    <- design$propErrVar
                             self$n             <- design$n
                             self$nObs          <- design$nObs
                             self$groups        <- design$groups
                             self$phases        <- design$phases
                             self$designMat     <- design$designMat
                             self$meanNames     <- design$meanNames
                             self$varNames      <- design$varNames
                             rm(design)

                             if(!all(names(self$phases)==self$phaseNames))
                             {
                               names(self$phases) <- self$phaseNames
                             }

                             self$designMat <- makeDesignMat(self$phases,
                                                             self$phaseNames,
                                                             self$randFxOrder,
                                                             'polyICT')

                             self$nObs <- length(c(unlist(self$phases)))
                           }

                           # reset warnings
                           options(warn=0)
                           #self$edit <- FALSE

                           # return self to allow for chaining of method calls,
                           # as well as lazy cloning (though this must be
                           # tested)
                           invisible(self)
                         },

                         # makeData method
                         makeData = function(seed=123)
                         {
                           seeds <- .makeSeeds(seed, length(self$phaseNames) *
                                                length(self$groupNames))
                           seeds <- matrix(seeds, length(self$phaseNames),
                                           length(self$groupNames))
                           data <- list(); d <- 1
                           for(p in seq_along(self$phaseNames))
                           {
                             thisp <- self$phaseNames[[p]]
                             for(g in seq_along(self$groupNames))
                             {
                               thisg <- self$groupNames[[g]]
                               # Sigma should be CorMat, not CovMat, otherwhise
                               # the slope variance (and higher polynomial terms)
                               # get scaled twice
                               Sigma <- self$randFxCorMat[[thisp]][[thisg]]
                               rFxVr <- self$inputMat[
                                 self$inputMat$Phase==thisp &
                                   self$inputMat$Group==thisg,
                                 self$varNames]
                               # multiply mu by sqrt(rFxVr) to get unstandardized
                               # means
                               mu    <- self$inputMat[
                                 self$inputMat$Phase==thisp &
                                   self$inputMat$Group==thisg,
                                 self$meanNames]*sqrt(rFxVr)
                               n     <- self$groups[[thisg]]
                               nObs  <- length(self$phases[[thisp]])
                               dM    <- self$designMat[self$designMat$phase==thisp,]
                               propErrVar <- self$inputMat[
                                 self$inputMat$Phase==thisp &
                                   self$inputMat$Group==thisg,
                                 c('randFx', 'res', 'mserr')]

                               # calls to .polyData() here
                               data[[d]] <- .polyData(seed       = seeds[p,g]         ,
                                                     n          = n                  ,
                                                     nObs       = nObs               ,
                                                     mu         = unlist(mu)         ,
                                                     Sigma      = Sigma              ,
                                                     self       = self               ,
                                                     dM         = dM                 ,
                                                     rFxVr      = unlist(rFxVr)      ,
                                                     propErrVar = unlist(propErrVar) ,
                                                     group      = thisg              )
                               d <- d + 1
                             }
                           }
                           data <- do.call(rbind, data)

                           # rescale y
                           data$y <- scale(data$y)*self$ySD + self$yMean

                           # the makeData method is terminal, self is not returned
                           # which means chaining is not possible past this point
                           return(data)

                         }

                       )
)


# TODO:need to generalize beyond slopes
#' expectedVar
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
# NOTE: randFxCovMat IS DEPRECATED
expectedVar <- function(randFxCovMat, designMat, variances,
                        randFxMean, nObs, n)
{
  # 'total' N
  N <- nObs * n

  # variance at each time point
  vyt <-  randFxCovMat[1,1]            +
    designMat$Time^2*randFxCovMat[2,2] +
    2*randFxCovMat[1,2]                +
    variances$errorVar

  # total variance due to random effects and error variance
  sigmaT <- (N-1)^-1 * ( (n - 1) * sum(vyt) )

  # rescale the fixed effect sizes in terms of the variance when time = 1
  whereIsTeq1 <- which(designMat$Time==1)
  sdYt1 <- sqrt( vyt[whereIsTeq1] )
  randFxMean$fixdFx$phase     <- randFxMean$fixdFx$phase     * sdYt1
  randFxMean$fixdFx$phaseTime <- randFxMean$fixdFx$phaseTime * sdYt1

  # mean at each time point
  mt <-  (randFxMean$randFx$intercept                         +
            randFxMean$randFx$slope     * designMat$Time      +
            randFxMean$fixdFx$phase     * designMat$phase     +
            randFxMean$fixdFx$phaseTime * designMat$phaseTime
  )

  # overall mean
  m = mean(mt)

  # total variance due to fixed effects stacked over time
  muT <- (N-1)^-1 * ( n * sum((mt-m)^2) )

  # overall variance
  vy     <- sigmaT + muT
  #  (N-1)^-1 * ( n * sum((mt-m)^2) + (n - 1) * sum(vyt) )


  return( list(randFxMean = randFxMean ,
               Variances   = vyt         ,
               Means       =  mt         ,
               TotalMean   =   m         ,
               TotalVar    =  vy         ) )
}

#' checkRandFxMean
#' @author Stephen Tueller \email{stueller@@rti.org}
#'
#' @keywords internal
#'
checkRandFxMeanPoly <- function(randFxMean)
{
  if(! var( unlist( lapply(randFxMean, length) ) ) == 0)
  {
    stop('All phases must have the same number of groups.')
  }

  if(! var(unlist( lapply(randFxMean, function(x) lapply(x, length)) ) ) == 0 )
  {
    stop('All combinations of phase and group must have the same number of\n',
         'random effects.')
  }
}






