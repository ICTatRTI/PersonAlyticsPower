library(PersonAlyticsPower)
library(ggplot2)
library(longCatEDA)
library(gamlss)
binICT <- polyICT$new(
  groups            = c(group1=100)                    ,
  phases            = makePhase(c(20, 20))             ,
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
binICT$inputMat[binICT$inputMat$Phase=='phase2', 'Mean0'] <- .1
y <- binICT$makeData(y0atyMean=F) # needed to prevent mean centering of y

# viz
hist(y$y)
# line plot - fails
ggplot(y, aes(x=Time, y=y, group=id, col=phase)) + geom_point() + geom_line()
# see ?sorter from longCatEDA
ggplot(y, aes(x=Time, y=id, fill=factor(y))) +
  geom_tile(colour="transparent") +
  scale_fill_manual(values=1:2) +
  facet_grid(phase ~ ., space="free_y", scales="free_y")
#lc <- longCat()


binICT$designCheck() # this fails

# https://www.escal.site/ says d=1 is OR=1.8
y <- byPhasebyGroup(y, time="Time", phase="phase", group="group")
(frm <- paste("y  ~ - 1 +", paste(y$dummyNames, collapse = "+")))
userFormula <- list()
userFormula$fixed <- formula(frm)
userFormula$random <-  formula(~ Time | id)
m0 <- PersonAlytic(output = 'BinaryTest0',
                   data = y$data,
                   ids = 'id',
                   dvs = 'y',
                   phase = 'phase',
                   time = 'Time',
                   ivs = y$dummyNames,
                   package = 'gamlss',
                   family = BI(),
                   autoSelect = list(),
                   userFormula = userFormula,
                   individual_mods  = FALSE)
summary(m0)

# QC model manually
(frm <- paste("y  ~ - 1 +", paste(y$dummyNames, collapse = "+"),
              "+ re(random=~ Time | id)"))
m0m <- gamlss(formula(frm), data=y$data,
              family=BI())
summary(m0m)
# drop the terms that should be 0
(frm <- paste("y  ~ - 1 + group1_phase2_int" ,
              "+ re(random=~ Time | id)"))
m0m <- gamlss(formula(frm), data=y$data,
              family=BI())
summary(m0m)

# update userFormula
userFormula <- termsToFormula("group1_phase2_int")
ICTpower(outFile = c('BinaryTest0', 'csv'),
         design = binICT,
         B = 100,
         seed = 3,
         prompt = F,
         y0atyMean = F, # needed to prevent mean centering of y, BI() is automatically selected
         family = BI(),
         userFormula = userFormula
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
ynb <- countICT$makeData(y0atyMean=F)
hist(ynb$y)
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

