# extended example illustrating phase and group specific options

# load a base example
library(PersonAlyticsPower)
example(polyICT)
myPolyICT$inputMat
myPolyICT$yMean <- 0
myPolyICT$ySD <- 1
myPolyICT$inputMat[6,6] <- -.3
myPolyICT$inputMat[3,3] <- 10

# check the design
myPolyICT$designCheck()

# set the number of bootstrap replications
B <- 1000

# no userFormula, no interactions, default growth model
nofrm <- ICTpower(c('nofrm', 'csv'),
                  myPolyICT,
                  B=B,
                  seed = 25,
                  prompt=FALSE,
                  interactions = list())

# two methods for adding interaction
# method 1: manually via userFormula
frm <- list()
frm$fixed <- y1 ~ group+phase+Time + group*Time + phase*Time + group*phase + group*phase*Time
frm$random <- ~Time | id
manualfrm <- ICTpower(c('manualfrm', 'csv'),
                     myPolyICT,
                     B=3,
                     seed = 25,
                     prompt=FALSE,
                     userFormula=frm,
                     interactions = list())

# method 2: with `interactions` option
manualfrm <- ICTpower(c('manualfrm', 'csv'),
                      myPolyICT,
                      B=B,
                      seed = 25,
                      prompt=FALSE,
                      interactions = list(c("group", "Time"),
                                          c("phase", "Time"),
                                          c("group", "phase")))

# two illustrations of adding a userformula with phase and group specific
# intercepts and slopes
# Example 1: all phase and group specific intercepts and slopes
userFormula <- byPhasebyGroup(NULL, 1:myPolyICT$nObs, myPolyICT$phaseNames,
                              myPolyICT$groupNames)
userFormula <- list(fixed=formula(paste("y1~", userFormula$fixed)),
                    random=~Time|id)
userfrm <- ICTpower(c('userfrm', 'csv'),
                    myPolyICT,
                    B=B,
                    seed = 25,
                    prompt=FALSE,
                    userFormula=userFormula,
                    interactions = list())

# Example 2: only the non-zero intercepts and slopes as implied in myPolyICT$inputMat
myPolyICT$inputMat
userFormula2 <- userFormula
userFormula2$fixed <- y1 ~ -1 + group2_phase2_int +
  group2_phase3_int + group2_phase3_slope
userfrm2 <- ICTpower(c('userfrm2', 'csv'),
                    myPolyICT,
                    B=B,
                    seed = 25,
                    prompt=FALSE,
                    userFormula=userFormula2,
                    interactions = list())

