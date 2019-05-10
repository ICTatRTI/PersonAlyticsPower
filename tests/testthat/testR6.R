context("R6")
library(PersonAlyticsPower)

test_that("polyICT updates",
          {
            # starting here we can mix regular code with tests

            # test class
            test1 <- polyICT$new()
            expect_equal(class(test1), c("polyICT", "designICT", "R6" ))

            # test updateability of test1
            .groups <- test1$groups <- c(group1=201, group2=202)
            names(.groups) <- names(test1$groups)
            expect_equal(test1$groups, .groups)


            .phases <- makePhase(c(10,30,20))
            test1$update(phases=.phases)
            expect_equal(test1$phases, .phases)

            .propErrVar <- c(randFx=.4,res=.3,mserr=.3)
            test1$update(propErrVar=.propErrVar)
            expect_equal(test1$propErrVar, .propErrVar)

            .randFxOrder <- 2
            .randFxVar <- c(1,1,1)
            test1$update(randFxOrder=.randFxOrder, randFxVar=.randFxVar)
            expect_equal(test1$randFxOrder, .randFxOrder)
            expect_equal(test1$randFxVar, .randFxVar)

            .randFxCor <- .5
            test1$update(randFxCor=.randFxCor)
            expect_equal(test1$randFxCor, .randFxCor)

            .error <- armaErr$new(list(ar=c(.5, -.2)))
            test1$update(error=.error)
            expect_equal(test1$error, .error)

            .merror <- armaErr$new(list(ar=c(0)))
            test1$update(merror=.merror)
            expect_equal(test1$merror, .merror)

            .ySD <- 1
            test1$update(ySD=.ySD)
            expect_equal(test1$ySD, .ySD)

            .yMean <- 0
            test1$update(yMean=.yMean)
            expect_equal(test1$yMean, .yMean)


          })




