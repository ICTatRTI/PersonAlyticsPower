context("Err")
library(PersonAlyticsPower)

test_that("dists",
{
  wei2_0_0 <- armaErr$new(model = list(ar=c(.5), ma=c(.5)), fam = 'WEI2',
                      famParms = list(mu=1, sigma=0.1)
  )

  wei2_2_.5 <- armaErr$new(model = list(ar=c(.5), ma=c(.5)), fam = 'WEI2',
                        famParms = list(mu=2, sigma=.5)
  )

  d1 <- wei2_0_0$checkModel()
  d2 <- wei2_2_.5$checkModel()

  testthat::expect_false(identical(d1, d2))

  #file.remove('Rplots.pdf')
  file.remove('stationary.arma')
})
