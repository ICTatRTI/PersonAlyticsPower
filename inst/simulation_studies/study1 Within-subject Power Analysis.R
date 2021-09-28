# code for

# Power Analysis for Idiographic (Within-Subject) Clinical Trials:
# Implications for Treatments of Rare Conditions and Precision Medicine

# Stephen Tueller, Derek Ramirez, Jessica D. Cance, Ai Ye, Anne C. Wheeler,
# Zheng Fan, Christoph Hornik & Ty A. Ridenour

# Correspondence concerning this article should be addressed to Ty Ridenour,
# RTI International, 3040 E. Cornwallis Rd., PO Box 12194, 326 Cox Bldg,
# Research Triangle Park, NC 27709-2194, email: tridenour@rti.org


# download and install PersonAlytics v0.2.6.8 from
# https://github.com/ICTatRTI/PersonAlytics/releases/tag/v0.2.6.8

# download and install PersonAlyticsPower v0.1.6.0 from
# https://github.com/ICTatRTI/PersonAlyticsPower/releases/tag/v0.1.6.0

# download and install PersonAlyticsSim v0.1.2.7 from
# https://github.com/ICTatRTI/PersonAlyticsSim/releases/tag/v0.1.2.8

# ensure you are in the desired directory first, setting up the study
# generates ~7000 subfolders

library(PersonAlyticsSim)

  groups <- list(p3500 = c(group1=3500),
                 p1750 = c(group1=1750),
                  p200 = c(group1=200) )

  sampSizes <- c(20,30,50,75,100)
  names(sampSizes) <- paste('n', sampSizes, sep='')

  nObs   <- c(30,50,75,100)
  propBL <- c(5,10,20)/100
  # check, only set up with a minumum of 5 bl observations
  lapply(as.list(nObs), function(x) x*propBL)

  phases <- list(
    # nObs = 20
    p5_15  = c(5, 15),
    # nObs = 40
    p5_35  = c(5, 35),
    p8_32  = c(8, 32),
    # nObs = 60
    p5_55  = c(5, 55),
    p6_54  = c(6, 54),
    p12_48 = c(12, 48),
    # nObs = 80
    p5_75  = c(5, 75),
    p8_72  = c(8, 72),
    p16_64 = c(16, 64),
    # nObs = 100
    p5_95  = c(5, 95),
    p10_90 = c(10, 90),
    p20_80 = c(20, 80)
  )
  unlist(lapply(phases, sum))

  propErrVar <- list(
    pev10 = c(randFx=.90, res=.05 , mserr=.05),
    pev25 = c(randFx=.75, res=.125, mserr=.125),
    pev50 = c(randFx=.50, res=.25 , mserr=.25)
  )
  unlist(lapply(propErrVar, sum))

  d <- list(
    # jump
    smlJump = matrix(c(0,.2,0,0), 2, 2),
    medJump = matrix(c(0,.5,0,0), 2, 2),
    lrgJump = matrix(c(0,.8,0,0), 2, 2),

    # slope - multiply by 2 so that the average difference is the effect size
    smlSlope = matrix(c(0,0,0,.2*2), 2, 2),
    medSlope = matrix(c(0,0,0,.5*2), 2, 2),
    lrgSlope = matrix(c(0,0,0,.8*2), 2, 2)
  )

  alignPhase  <- c(mmta = 'none', piecewise = 'piecewise')

  # if running on multiple computers, list their names here, otherwise use
  # noneNames <- list()
  nodeNames <- list(#sdesk   = c(name="RTI-102807", cores=8),
                    sgreen  = c(name="RTI-104040", cores=8),
                    sorange = c(name="RTI-102921", cores=8)
  )

  # directory set up
  dirNames <- list(c("d", "alignPhase"), c("groups", "sampSizes", "phases", "propErrVar"))

  # aggregate conditions
  conditions <- list(
    groups     = groups     ,
    sampSizes  = sampSizes ,
    phases     = phases    ,
    propErrVar = propErrVar  ,
    d          = d          ,
    alignPhase = alignPhase
  )

  # set up the simulation with 3 replications per condition to prevent
  # overloading cpu if run accidentally prior to a full run
  t1 <- PersonAlyticsSim:::Sim$new(
    simName    = "t1"      ,
    seed       = 1         ,
    B          = 3         ,
    conditions = conditions,
    dirNames   = dirNames  ,
    type       = 'nonpar'  ,
    nodeNames  = nodeNames )
  t1 # study summary
  save(t1, file='study1.RData')

  # full study with 1000 replications per condition
  library(PersonAlyticsSim)
  library(gridExtra)
  load('study1.RData')
  setwd(t1$wd)
  t1$seed # check seed
  t1$progress()
  t1$B <- 1000
  t1$start()
  t1$summarize()

