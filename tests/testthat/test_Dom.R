context("Domenico 1987")
library(AnTrans)

# wrapper with convenient defaults
cD_def <- function(C0 = 1, x = 0, y = 0, z = 0, t = 1, vx = 1,
                   ax = 1, ay = 1, az = 1, lambda = 0, sY = 1, sZ = 1,
                   Rf = 1, decaysorbed = FALSE){
  const_Dom(C0, x, y, z, t, vx, ax, ay, az, lambda, sY, sZ,
            Rf, decaysorbed)
}

test_that("Known values", {
  # at the source, plug flow
  expect_equal(cD_def(ax = 0), 1)

  # at the source, very progressed plume (source at steady state)
  expect_equal(cD_def(t = 1e5), 1)

  # plug flow, before and after vt
  expect_equal(cD_def(t = 1, vx = 1, ax = 0, ay = 0, az = 0,
                      x = c(.9, 1.1)), c(1, 0))

  # longitudinal dispersion zone
  expect_gt(cD_def(t = 1, vx = 1, ax = 0, ay = 0, az = 0, x = .9), .5)
  expect_lt(cD_def(t = 1, vx = 1, ax = 0, ay = 0, az = 0, x = 1.1), .5)

  # no transverse dispersion and not downstream of source
  expect_equal(cD_def(t = 50, ay = 0, az = 0, x = rep(c(1, 10), each = 4L),
                      y = rep(c(.6, 0, -.6, 0), times = 2L),
                      z = rep(c(0, .6, 0, -.6), times = 2L),
                      sY = 1, sZ = 1),
               rep(0, 8L))

  # declining centre line with transverse dispersion
  expect_true(all(diff(cD_def(x = 0:10, t = 20, ax = 0, ay = 1, az = 1)) < 0))

  # degradation reduces concentration
  expect_lt(cD_def(x = 5, y = 3, z = 2, t = 10, lambda = .1),
            cD_def(x = 5, y = 3, z = 2, t = 10, lambda = 0))

  # retardation with plug flow
  expect_equal(cD_def(t = 1, vx = 1, ax = 0, ay = 0, az = 0,
                      x = c(.9, 1.1)/2, Rf = 2), c(1, 0))
  #
  # - decay and decaysorbed
  expect_lt(cD_def(t = 1, vx = 1, ax = .1, ay = 0, az = 0, lambda = .1,
                   x = .9/2, Rf = 2, decaysorbed = FALSE),
            cD_def(t = 1, vx = 1, ax = .1, ay = 0, az = 0, lambda = 0,
                   x = .9/2, Rf = 2, decaysorbed = FALSE))
  expect_gt(cD_def(t = 1, vx = 1, ax = .1, ay = 0, az = 0, lambda = .1,
                   x = .9/2, Rf = 2, decaysorbed = FALSE),
            cD_def(t = 1, vx = 1, ax = .1, ay = 0, az = 0, lambda = .1,
                   x = .9/2, Rf = 2, decaysorbed = TRUE))

  # special case equivalence to O-B
  t <- sample(1:50, 1L)
  ax <- sample(seq(.1, 3, .1), 1L)
  x <- 0:20
  expect_equal(cD_def(t = t, x = x, ax = ax, ay = 0, az = 0, lambda = 0),
               .5*pracma::erfc((x - 1*t)/(2*(ax*1*t)^.5)),
               tolerance = 1e-10)
})
