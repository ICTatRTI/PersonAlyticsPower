library(PersonAlyticsPower)
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
  yCut = c(.5, .5)
)
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase2' &
                     myPolyICT$inputMat$Group=='group2', 'Mean0'] <- .3
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
                     myPolyICT$inputMat$Group=='group2', 'Mean0'] <- .3
myPolyICT$inputMat[myPolyICT$inputMat$Phase=='phase3' &
                     myPolyICT$inputMat$Group=='groups', 'Mean1'] <- -.6
myPolyICT$inputMat
y <- myPolyICT$makeData(y0atyMean=F) # needed to prevent mean centering of y
hist(y$y)
myPolyICT$designCheck() # this fails


m0 <- PersonAlytic(output = 'BinaryTest0',
                   data = y,
                   ids = 'id',
                   dvs = 'y',
                   phase = 'phase',
                   time = 'Time',
                   ivs = c('group'),
                   package = 'gamlss',
                   family = BI(),
                   autoSelect = list(),
                   individual_mods  = FALSE)
summary(m0)

ICTpower(outFile = c('BinaryTest0', 'csv'),
         design = myPolyICT,
         B = 3,
         seed = 3,
         prompt = F,
         y0atyMean = F, # needed to prevent mean centering of y, BI() is automatically selected
         family = BI()
         )

 # try "count" variable
nb <- rnbinom(nrow(y), mu = 1, size = 1)
countICT <- polyICT$new(
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
  yCut = c(table(nb)/length(nb))
)
ynb <- countICT$makeData()
hist(ynb$y)
table(ynb$y)
ynb$y <- as.numeric(factor(ynb$y)) - 1
table(ynb$y)

nb0 <- PersonAlytic(output = 'NegBinTest0',
                   data = y,
                   ids = 'id',
                   dvs = 'y',
                   phase = 'phase',
                   time = 'Time',
                   ivs = c('group'),
                   package = 'gamlss',
                   family = NBI(),
                   autoSelect = list(),
                   individual_mods  = FALSE)
summary(nb0)

