if(1==2)
{

  # load the library
  library(PersonAlyticsPower)

  # set up the error first and do a check, if the check fails a text file
  # will be produce called 'stationary.arma' with the value 0 in it, as an
  # example, this code yields an error:

  #error <- armaErr$new(list(ar=c(.5,.5),
  #                          ma=c(.3,.3)))
  #error$checkModel(doPlot = FALSE)

  error <- armaErr$new(list(ar=c(.5,-.2),
                            ma=c(.3,.3)))
  error$checkModel(doPlot = FALSE)

  # create the base design
  design <- polyICT$new(
    groups            = c(group1=10, group2=10)                   ,
    phases            = makePhase(c(5,10,5))                      ,
    propErrVar        = c(randFx=.5,res=.5,mserr=.0)              ,
    randFxOrder       = 1                                         ,
    randFxCor         = 0.2                                       ,
    randFxVar         = c(1, 1)                                   ,

    error             = error                                     ,
    merror            = armaErr$new(list())                       ,
    ySD               = 15                                        ,
    yMean             = 100                                       ,
    yMin              = NULL                                      ,
    yMax              = NULL                                      ,
    yCut              = NULL                                      ,
  )

  # now the effect sizes need to be added
  # first save in a format the gui can read in to display and allow the
  # user to edit values
  write.csv(design$inputMat, 'inputMat.csv')
  # then overwrite with the user's edits
  design$inputMat <- read.csv('inputMat.csv')

  # provide the user with a visual check of their design
  png("designCheck.png")
  design$designCheck(return = "plot")
  dev.off()

  # if they do not like the design, open inputMat again
  # if they do like the design, run ICTpower

  # run the power analysis which returns
  ICTpower <- function(outFile         = "GUItest"                 ,
                       design          = design                    ,
                       B               = 10                        ,
                       alpha           = .05                       ,
                       seed            = 123                       ,
                       cores           = parallel::detectCores()-1 )


}
