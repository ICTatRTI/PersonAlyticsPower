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
           alignPhase  = 'piecewise'          ,
           prompt      = FALSE                ,
           debugforeach = TRUE                )

  # we do get convergence on the full data set
  pweg <- Palytic$new(data = Data, ids = 'id', dv = 'y', time = 'Time',
                      phase = 'phase', alignPhase = 'piecewise', ivs='group')
  pweg$lme()

  # this fails, no convergence
  ICTpower(outFile     = c("piecewise", "csv"),
           B           = 3                    ,
           dataFile    = "Data.RData"         ,
           sampleSizes = c(100,100)           ,
           seed        = 302                  ,
           alignPhase  = 'piecewise'          ,
           prompt      = FALSE                )

  # manual palytic works...
  rm(Data); Data <- read.csv("piecewise.Data.csv")

  y1 <- Palytic$new(data = Data , ids = "id", dv = "y1", time = "Time",
                    phase = "phase", ivs = "group", alignPhase = "piecewise")
  y1$lme()$PalyticSummary$formula
  y1$correlation
  y1$correlation <- "corARMA(p=1,q=1)"
  y1$correlation
  y1$formula
  y1lme <- y1$lme()
  y1lme$PalyticSummary$formula
  y1lme$whichPalyticMod
  y1lme$PalyticSummary$correlation # why is correlation being overwritten??
  #y1$time <- "Time"
  #y1$plot()

  y2 <- Palytic$new(data = Data , ids = "id", dv = "y2", time = "Time",
                    phase = "phase", ivs = "group", alignPhase = "piecewise")
  y2$lme()
  #y2$time <- "Time"
  #y2$plot()

  y3 <- Palytic$new(data = Data , ids = "id", dv = "y3", time = "Time",
                    phase = "phase", ivs = "group", alignPhase = "piecewise")
  y3$lme()
  #y3$time <- "Time"
  #y3$plot()

  file.remove("Data.RData")

})

