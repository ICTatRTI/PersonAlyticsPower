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


compare <- function(fam, cuts=c(.5, .5), seed=42, d=1:10)
{
  des <- polyICT$new(
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
    yCut = cuts
  )

  set.seed(seed)
  seeds <- runif(length(d), 0, 7e7)

  out <- list()
  for(i in d)
  {
    des2 <- des$clone()
    des2$inputMat[binICT$inputMat$Phase=='phase2', 'Mean0'] <- i/10
    y <- des2$makeData(y0atyMean=F, seed = seeds[i]) # needed to prevent mean centering of y
    y <- byPhasebyGroup(y, time="Time", phase="phase", group="group")
    dat <<- data.frame(y$data)

    # intercept only model (to check data generation)
    m00 <- gamlss(y ~ 1, data=dat, family=fam)
    tab0 <- summary(m00)

    # PersonAlytics model
    frm <- paste("y  ~ - 1 + group1_phase2_int" ,
                  "+ re(random=~ Time | id)")
    m0m <- gamlss(formula(frm), data=dat,
                  family=fam, control = gamlss.control(n.cyc = 200,
                                                        trace=F,
                                                        c.crit = 0.05))
    tab <- summary(m0m)
    out[[i]] <- c(d=i/10, est=tab[1,1], OR=exp(tab[1,1]), p=tab[1,4],
                  mu=exp(tab0[1,1]), sigma=exp(tab0[2,1]))
    rm(m0m)
  }
  do.call(rbind, out)
}

# median split
median_split <- compare(BI(), cuts=c(.5, .5))

# rare event
round(split_90_10 <- compare(BI(), cuts=c(.90, .10)), 2)

# "count"
# rbinom approach
nb <- rnbinom(500, mu = 1, size = 1)
yCut = c(table(nb)/length(nb))
round(split_count <- compare(NBI(), cuts=yCut), 2)

# rNBI approach
nbi1 <- rNBI(500, mu = 1)
(yCut1 = c(table(nbi1)/length(nbi1)))
round(split_count1 <- compare(NBI(), cuts=yCut1), 2)

nbi2 <- rNBI(500, mu = 2)
(yCut2 = c(table(nbi2)/length(nbi2))); yCut1
round(split_count2 <- compare(NBI(), cuts=yCut2), 2)

nbi12 <- rNBI(500, mu = 12)
(yCut12 = c(table(nbi12)/length(nbi12)))
round(split_count12 <- compare(NBI(), cuts=yCut12), 2)

set.seed(98742)
nbi5.0.5 <- rNBI(1000, mu=5, sigma=0.5); hist(nbi5.0.5)
(yCut5.0.5 = c(table(nbi5.0.5)/length(nbi5.0.5)))
round(split_count5.0.5 <- compare(NBI(), cuts=yCut5.0.5), 2)
round(split_count5.0.5_2345 <- compare(NBI(), cuts=yCut5.0.5, seed = 2345), 2)



# ICTpower directly
des_nb <- polyICT$new(
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
  yCut = yCut5.0.5
)
frm <- termsToFormula("group1_phase2_int")
ICTpower(outFile = "des_nb",
         design = des_nb,
         B = 3,
         userFormula = frm,
         family = NBI())

