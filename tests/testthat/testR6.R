context("R6")
library(PersonAlyticsPower)

test_that("polyICT updates",
          {
            # starting here we can mix regular code with tests

            # test class
            test1 <- polyICT$new()
            testthat::expect_equal(class(test1), c("polyICT", "designICT", "R6" ))

            # test updateability of test1
            .groups <- test1$inputMat$n[c(1,4)] <- c(group1=201, group2=202)
            names(.groups) <- names(test1$groups)
            testthat::expect_equal(test1$groups, .groups)


            .phases <- makePhase(c(10,30,20))
            test1$inputMat$nObs[1:3] <- c(10,30,20)
            testthat::expect_equal(test1$phases, .phases)

            #QC updating and reopen these tests
            #.propErrVar <- c(randFx=.4,res=.3,mserr=.3)
            #test1$update(propErrVar=.propErrVar)
            #testthat::expect_equal(test1$propErrVar, .propErrVar)

            #.randFxOrder <- 2
            #.randFxVar <- c(1,1,1)
            #test1$update(randFxOrder=.randFxOrder, randFxVar=.randFxVar)
            #testthat::expect_equal(test1$randFxOrder, .randFxOrder)
            #testthat::expect_equal(test1$randFxVar, .randFxVar)

            #.randFxCor <- .5
            #test1$update(randFxCor=.randFxCor)
            #testthat::expect_equal(test1$randFxCor, .randFxCor)

            .error <- armaErr$new(list(ar=c(.5, -.2)))
            test1$update(error=.error)
            testthat::expect_equal(test1$error, .error)

            .merror <- armaErr$new(list(ar=c(0)))
            test1$update(merror=.merror)
            testthat::expect_equal(test1$merror, .merror)

            .ySD <- 1
            test1$update(ySD=.ySD)
            testthat::expect_equal(test1$ySD, .ySD)

            .yMean <- 0
            test1$update(yMean=.yMean)
            testthat::expect_equal(test1$yMean, .yMean)


          })




