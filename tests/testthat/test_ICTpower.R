context("ICTpower")
library(PersonAlyticsPower)

test_that("piecewise",
{
  example(polyICT)
  myPolyICTnonPar <- myPolyICT$clone(deep=TRUE)
  myPolyICTnonPar$update(groups=c(group1=500, group2=500))
  Data <- myPolyICTnonPar$makeData()
  save(Data, file = "Data.RData")

  # this fails, no convergence
  ICTpower(outFile     = c("piecewise", "csv"),
           B           = 3                    ,
           dataFile    = "Data.RData"         ,
           sampleSizes = c(25,25)             ,
           alignPhase  = 'piecewise'          )

  # we do get convergence on the full data set
  pweg <- Palytic$new(data = Data, ids = 'id', dv = 'y', time = 'Time',
                      phase = 'phase', alignPhase = 'piecewise', ivs='group')
  pweg$lme()

  # this fails, no convergence
  ICTpower(outFile     = c("piecewise", "csv"),
           B           = 3                    ,
           dataFile    = "Data.RData"         ,
           sampleSizes = c(400,400)           ,
           alignPhase  = 'piecewise'          )



  file.remove("Data.RData")

})

