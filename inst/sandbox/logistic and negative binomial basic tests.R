library(PersonAlyticsPower)
library(ggplot2)
library(longCatEDA)
library(gamlss)


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
    frm <- paste("y  ~ - 1 + group1_phase2_int" ,
                  "+ re(random=~ Time | id)")
    m0m <- gamlss(formula(frm), data=y$data,
                  family=fam, control = gamlss.control(n.cyc = 200,
                                                        trace=F,
                                                        c.crit = 0.05))
    tab <- summary(m0m)
    out[[i]] <- c(d=i/10, est=tab[1,1], OR=exp(tab[1,1]), p=tab[1,4])
    rm(m0m)
  }
  do.call(rbind, out)
}

# median split
median_split <- compare(BI(), cuts=c(.5, .5))

# rare event
round(split_90_10 <- compare(BI(), cuts=c(.90, .10)), 2)

# "count"
nb <- rnbinom(500, mu = 1, size = 1)
yCut = c(table(nb)/length(nb))
round(split_count <- compare(NBI(), cuts=yCut), 2)


