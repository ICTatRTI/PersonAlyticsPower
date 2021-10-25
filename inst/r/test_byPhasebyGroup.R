library(PersonAlyticsPower)
example(polyICT)
myPolyICT$inputMat

head(myPolyICT$makeData())

userFormula <- byPhasebyGroup(NULL, 1:myPolyICT$nObs, myPolyICT$phaseNames,
                              myPolyICT$groupNames)
userFormula <- list(fixed=formula(paste("y1~", userFormula$fixed)),
                    random=~Time|id)

debug(ICTpower)
userfrm <- ICTpower(c('userfrm', 'csv'),
                    myPolyICT,
                    B=1000,
                    seed = 25,
                    prompt=FALSE,
                    userFormula=userFormula,
                    interactions = list())

