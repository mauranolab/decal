test_that("estimate_theta is consistent with seed", {
  Y <- matrix(rpois(1e5, 100), ncol=1000)
  g <- which(rowMeans(Y) > 1)
  th1 <- withr::with_seed(42, {
    estimate_theta(Y, rowMeans(Y), log(colSums(Y)), genes = g)
  })
  th2 <- withr::with_seed(42, {
    estimate_theta(Y, rowMeans(Y), log(colSums(Y)), genes = g)
  })
  expect_equal(th1[g], th2[g], tolerance = 1e-5)
})
