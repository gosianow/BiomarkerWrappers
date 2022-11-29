context("duplicated2")

x <- c(1, 2, 3, 2)


test_that("duplicated2 finds all duplicated values", {
  expect_true(all(duplicated2(x) == c(FALSE, TRUE, FALSE, TRUE)))
})

