if(1==2)
{
  library(PersonAlyticsPower)

  design <- polyICT$new()

  ICTpower <- function(outFile         = "GUItest"                      ,
                       design          = NULL                      ,
                       B               = 100                       ,
                       dataFile        = NULL                      ,
                       sampleSizes     = NULL                      ,
                       alpha           = .05                       ,
                       seed            = 123                       ,
                       cores           = parallel::detectCores()-1 ,
                       savePowerReport = TRUE                      ,
                       standardize     = list(dv    = FALSE ,
                                              ivs   = FALSE ,
                                              byids = FALSE )      ,
                       prompt          = TRUE                      )



  myPolyICT <- polyICT$new(
    groups            = c(group1=10, group2=10)          ,
    phases            = makePhase()                      ,
    propErrVar        = c(randFx=.5, res=.25, mserr=.25) ,
    randFxOrder       = 1                                ,
    randFxCor         = 0.2                              ,
    randFxVar         = c(1, 1)                          ,
    error             = armaErr$new()                    ,
    merror            = armaErr$new(list())              ,
    yCut              = c(.5,.5)                         ,
  )
  myPolyICT$designCheck()
  dat <- myPolyICT$makeData()
  head(dat)
  table(dat$y)

  ggplot(data = dat, aes(x=y)) + geom_histogram()
}
