context("R6")
library(PersonAlyticsPower)

test_that("polyICT updates",
          {
            # starting here we can mix regular code with tests

            # test class
            test1 <- polyICT$new()
            expect_equal(class(test1), c("polyICT", "designICT", "R6" ))

            # test updateability of test1
            .n <- test1$n <- c(201,202)
            expect_equal(test1$n, .n)
            test1$update(n=.n)
            expect_equal(test1$n, .n)

          })