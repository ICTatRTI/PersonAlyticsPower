# 20190129 by Stephen Tueller

# check whether arima.sim is adequate for PersonAlyticsPower
# historical note: the prior work in PaCCT.Power
# R:\PaCCT\09 Repository\PaCCT.Power
# is complex and deeply looped, limited to a small number of implied designs,
# and difficult to communicate to the user how it can be used

# Following ALDA Chapter 7 (see esp. 7.3.4 & Table 7.3) can we get assymptotically
# equivalent data from
# 1. the autoregressive error covariance matrix
# 2. arima.sim
#
# why is this important? one or the other may be easier to implement




# currently just AR(1)
sigmaMaker <- function(n=20, phi=.5, sigma2=.75)
{
  Sigma <- matrix(as.numeric(NA), n, n)
  diag(Sigma) <- sigma2
  Sigma <- data.frame(Sigma)

  # define a matrix of diagonal bands
  delta <- row(Sigma) - col(Sigma)
  # make delta symmetric
  delta <- abs(delta)

  for(i in seq_along(data.frame(Sigma)))
  {
    Sigma[delta==i] <- sigma2 * phi^i
  }

  return(Sigma)
}

# use sigmaMaker and simulate some MVN data
library(MASS)


ar1sim <- function(B=20, n=20, phi=.5, sigma2=.75)
{
  Out <- list()

  for(b in 1:B)
  {
    set.seed(b)
    method1 <- arima.sim(model = list(ar=c(phi)), n = n, sd = sigma2)

    Sigma <- sigmaMaker(n = n, phi = phi, sigma2 = sigma2)
    set.seed(b^2)
    method2 <- mvrnorm(1, rep(0,nrow(Sigma)), Sigma)

    try(
      {
        method1 <- arima(method1, order = c(1,0,0))
        method2 <- arima(method2, order = c(1,0,0))

        out <- c(
          s1=method1$sigma2,
          s2=method2$sigma2,
          ar1=method1$coef[1],
          ar2=method2$coef[1],
          arse1=sqrt(diag(method1$var.coef)[1]),
          arse2=sqrt(diag(method2$var.coef)[1]))

        Out[[b]] <- out
      }  , silent = TRUE
    )

  }

  Out <- as.data.frame(do.call(rbind, Out))
  print(lapply(Out, mean))
  invisible(Out)
}

Out <- ar1sim(20000)
