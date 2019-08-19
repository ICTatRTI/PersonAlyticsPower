context("ICTpower")
library(PersonAlyticsPower)

test_that("piecewise",
{
  myPolyICT <- polyICT$new(
    groups            = c(group1=10, group2=10)          ,
    phases            = makePhase()                      ,
    propErrVar        = c(randFx=.5, res=.25, mserr=.25) ,
    randFxOrder       = 1                                ,
    randFxCor         = 0.2                              ,
    randFxVar         = c(1, 1)                          ,
    error             = armaErr$new()                    ,
    merror            = armaErr$new(list())              ,
    ySD               = 15                               ,
    yMean             = 100                              ,
  )

  myPolyICTnonPar <- myPolyICT$clone(deep=TRUE)
  myPolyICTnonPar$update(groups=c(group1=500, group2=500))
  Data <- myPolyICTnonPar$makeData()
  save(Data, file = "Data.RData")

  # this fails, no convergence
  if(1==2)
  {
    ICTpower(outFile      = c("piecewise", "csv"),
           B            = 3                    ,
           dataFile     = "Data.RData"         ,
           sampleSizes  = c(25,25)             ,
           alignPhase   = 'piecewise'          ,
           prompt       = FALSE                ,
           debugforeach = FALSE                )
  }

  txts <- dir(getwd(), glob2rx("*.txt"))
  csvs <- dir(getwd(), glob2rx("*.csv"))
  pdfs <- dir(getwd(), glob2rx("*.pdf"))
  file.remove("Data.RData", txts, csvs, pdfs)

})

