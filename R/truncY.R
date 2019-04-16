# try last answer from
# https://stats.stackexchange.com/questions/113230/generate-random-numbers-following-a-distribution-within-an-interval
# noting that this is very similar in structure to the approach taken in mvrFam() for
# transforming data

# note, we use doCall from R.utils to allow for ignore arguments (e.g., nu, tau)
truncY0 <- function(N, parms=list(mu=0, sigma=1, nu=2, tau=2),
                   DIST='NO', a = -Inf, b = Inf, seed=123)
{
  if(a > b) stop('`b` must be > `a`')

  # TODO this check is also in ICTviz, consider moving common checks to one function
  if( ! is(gamlss.dist::gamlss.family(DIST)) == "gamlss.family")
  {
    stop("The value of `DIST`=", DIST, " is not a `gamlss.family` distribution.")
  }

  pdist <- paste('p', DIST, sep='')
  qdist <- paste('q', DIST, sep='')

  set.seed(seed)
  p <- runif(N,
             doCall(pdist, q=a,
                    mu=parms$mu,
                    sigma=parms$sigma,
                    nu=parms$nu, tau=parms$tau),
             doCall(pdist, q=b,
                    mu=parms$mu,
                    sigma=parms$sigma,
                    nu=parms$nu, tau=parms$tau))

  doCall(qdist, p=p,
         mu=parms$mu,
         sigma=parms$sigma,
         nu=parms$nu, tau=parms$tau)

}

# TODO this needs to be a method for a class, we're moving the same
# parameters around too much
# UL <- designMatrix[,4:5]
truncY <- function(N, UL, DIST='NO', parms=list(mu=0, sigma=1, nu=2, tau=2),
                   seed = 123)
{
  # this approach can't handle varying...so we'll have to loop or lapply?
  trCall <- paste("gamlss.tr::gen.trun(par = UL, family = ", DIST,
                  ", type = 'both', varying = TRUE)")
  eval(parse(text=trCall))
}


if(1==2)
{
  #dY <- data.frame(y=dY)
  #library(ggplot2)
  #ggplot(dY, aes(x=y)) + geom_density()
  # produces error
  #library(ggmap)
  #gglocator()

  rm3 <- function(x) round(mean(x), 3)
  pdf('truncYexample_Normal.pdf', height = 11, width = 8.5)
  par(mfrow=c(4,2))
  y0 <- truncY(10000, a=-Inf, b=Inf, seed=1);  hist(y0, main=paste('Mean=', rm3(y0))); qqnorm(y0); qqline(y0, col=2)
  y1 <- truncY(10000, a=-4,   b=-2,  seed=2);  hist(y1, main=paste('Mean=', rm3(y1))); qqnorm(y1); qqline(y1, col=2)
  y2 <- truncY(10000, a=-3,   b=-1,  seed=3);  hist(y2, main=paste('Mean=', rm3(y2))); qqnorm(y2); qqline(y2, col=2)
  y3 <- truncY(10000, a=-2,   b= 0,  seed=4);  hist(y3, main=paste('Mean=', rm3(y3))); qqnorm(y3); qqline(y3, col=2)
  y4 <- truncY(10000, a=-1,   b= 1,  seed=5);  hist(y4, main=paste('Mean=', rm3(y4))); qqnorm(y4); qqline(y4, col=2)
  y5 <- truncY(10000, a= 0,   b= 2,  seed=6);  hist(y5, main=paste('Mean=', rm3(y5))); qqnorm(y5); qqline(y5, col=2)
  y6 <- truncY(10000, a= 1,   b= 3,  seed=7);  hist(y6, main=paste('Mean=', rm3(y6))); qqnorm(y6); qqline(y6, col=2)
  y7 <- truncY(10000, a= 2,   b= 4,  seed=8);  hist(y7, main=paste('Mean=', rm3(y7))); qqnorm(y7); qqline(y7, col=2)
  dev.off()

  par(mfrow=c(1,2))

  pdf('truncYexample_Skew.pdf', height = 11, width = 8.5)
  par(mfrow=c(4,2))
  #y0 <- truncY(100, pdist='pSEP1',           , seed=1);  hist(y0); qqnorm(y0); qqline(y0, col=2)
  hist(rSEP1(1000, mu=0, sigma=1, nu=0, tau=.2))
  y1 <- truncY(1000, pdist='pSEP1', a=-4, b=-2, seed=2);  hist(y1); qqnorm(y1); qqline(y1, col=2)
  y2 <- truncY(1000, pdist='pSEP1', a=-3, b=-1, seed=3);  hist(y2); qqnorm(y2); qqline(y2, col=2)
  y3 <- truncY(1000, pdist='pSEP1', a=-2, b= 0, seed=4);  hist(y3); qqnorm(y3); qqline(y3, col=2)
  y4 <- truncY(1000, pdist='pSEP1', a=-1, b= 1, seed=5);  hist(y4); qqnorm(y4); qqline(y4, col=2)
  y5 <- truncY(1000, pdist='pSEP1', a= 0, b= 2, seed=6);  hist(y5); qqnorm(y5); qqline(y5, col=2)
  y6 <- truncY(1000, pdist='pSEP1', a= 1, b= 3, seed=7);  hist(y6); qqnorm(y6); qqline(y6, col=2)
  y7 <- truncY(1000, pdist='pSEP1', a= 2, b= 4, seed=8);  hist(y7); qqnorm(y7); qqline(y7, col=2)
  dev.off()

}

  ## get user inputs
  ## TODO: we need error handling and restarts here!!
  #
  #cat('\n')
  #nGroups <- readinteger("Enter the number of groups")
  #groupSizes <- list()
  #for(g in 1:nGroups)
  #{
  #  groupName <- readcharacter("Enter the name of group", g)
  #  groupSizes[[groupName]] <- readinteger("Enter the # of participants in group",
  #                                     groupName)
  #}
  #
  #nPhases <- readinteger("Enter the number of phases: ")
  #phaseLengths <- list()
  #for(p in 1:nPhases)
  #{
  #  phaseName <- readcharacter("Enter the name of phase", p)
  #  phaseLengths[[ phaseName ]] <- readinteger( "Enter the # of timepoints in phase",
  #                                              phaseName    )
  #}
  #
  #locator(type = 'p', col = 'red')


