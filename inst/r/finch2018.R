finch2018 <- polyICT$new(
  groups            = c(group1=25, group2=25)          ,
  phases            = makePhase(c(5,4))                ,
  propErrVar        = c(randFx=.5, res=.5, mserr=.0)   ,
  randFxOrder       = 1                                ,
  randFxCor         = 0.1                              ,
  randFxVar         = c(1, 1)                          ,
  error             = armaErr$new()                    ,
  merror            = armaErr$new(list())              ,
  ySD               = 1                                ,
  yMean             = 0                                ,
)

finch2018$inputMat[finch2018$inputMat$Phase=='phase1' &
                     finch2018$inputMat$Group=='group1', 'Mean1'] <- .75/5

finch2018$inputMat[finch2018$inputMat$Phase=='phase2' &
                     finch2018$inputMat$Group=='group1', 'Mean0'] <- .75
finch2018$inputMat[finch2018$inputMat$Phase=='phase2' &
                     finch2018$inputMat$Group=='group1', 'Mean1'] <- .375/4

finch2018$inputMat[finch2018$inputMat$Phase=='phase1' &
                     finch2018$inputMat$Group=='group2', 'Mean1'] <- .25/4

finch2018$inputMat[finch2018$inputMat$Phase=='phase2' &
                     finch2018$inputMat$Group=='group2', 'Mean0'] <- .25/5
finch2018$inputMat[finch2018$inputMat$Phase=='phase2' &
                     finch2018$inputMat$Group=='group2', 'Mean1'] <- .125


finch2018$inputMat
debug(PersonAlyticsPower:::.polyData)
#debug(finch2018$designCheck)
debug(finch2018$makeData)
finch2018$designCheck()

interactions = list()
ints <- c("group1_phase1_slope",
          "group1_phase2_int"  , "group1_phase2_slope",
          "group2_phase1_slope",
          "group2_phase2_int"  , "group2_phase2_slope")
finch2018formula <- list()
finch2018formula$fixed <- formula(paste("y1 ~ -1 +", paste(ints, collapse = "+")))
finch2018formula$random <- formula(~Time | id)

alignPhase = "align"

ICTpower(outfile = c("Finch_et_al._2018", "csv"),
         design  = finch2018,
         userFormula = finch2018formula,
         interactions = interactions,
         alignPhase = alignPhase,
         B = 100,
         seed = 90423)


Y <- Reduce('+', .L)



