library(PersonAlyticsPower)
example(polyICT)
myPolyICT$inputMat
myPolyICT$yMean <- 0
myPolyICT$ySD <- 1
myPolyICT$inputMat[6,6] <- -.3
myPolyICT$inputMat[3,3] <- 10


temp <- myPolyICT$makeData()

myPolyICT$designCheck(npg = 100)

B <- 10#00

nofrm <- ICTpower(c('nofrm', 'csv'),
                  myPolyICT,
                  B=B,
                  seed = 25,
                  prompt=FALSE,
                  interactions = list())

# two methods for adding interaction
# method 1: manually
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

# method 2: with `interactions`
manualfrm <- ICTpower(c('manualfrm', 'csv'),
                      myPolyICT,
                      B=B,
                      seed = 25,
                      prompt=FALSE,
                      #userFormula=frm,
                      interactions = list(c("group", "Time"),
                                          c("phase", "Time"),
                                          c("group", "phase")))

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

out <- read.csv("./userfrm_PersonAlytic.csv")
out2 <- read.csv("./userfrm2_PersonAlytic.csv")
mean(out$AIC)
mean(out2$AIC)
library(coxed)
bca(out$AIC)
bca(out2$AIC)
