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
  ICTpower(outFile      = c("piecewise", "csv"),
           B            = 3                    ,
           dataFile     = "Data.RData"         ,
           sampleSizes  = c(25,25)             ,
           alignPhase   = 'piecewise'          ,
           prompt       = FALSE                ,
           debugforeach = TRUE                 )



  file.remove("Data.RData")

})

