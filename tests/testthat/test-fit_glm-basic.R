# End-to-end smoke tests of basic GLM functionality for one formula

test_that("linear-basic", {
  expect_error( # expect no error (smoke test)
    lin_fit <- fit_glm(
      formula             = y ~ X1 + X2,
      p_family            = "linear",
      data                = "random",
      iterations          = 10,
      warmup_iterations   = 5,
      chains              = 1, 
      seed                = 123
    ),
    NA
  )

})

test_that("logistic-basic", {
  expect_error( # expect no error (smoke test)
    lin_fit <- fit_glm(
      formula             = y ~ X1 + X2,
      p_family            = "logistic",
      data                = "random",
      iterations          = 10,
      warmup_iterations   = 5,
      chains              = 1, 
      seed                = 123
    ),
    NA
  )
})

test_that("poisson-basic", {
  expect_error( # expect no error (smoke test)
    lin_fit <- fit_glm(
      formula             = y ~ X1 + X2,
      p_family            = "poisson",
      data                = "random",
      iterations          = 10,
      warmup_iterations   = 5,
      chains              = 1, 
      seed                = 123
    ),
    NA
  )
})

test_that("gamma-basic", {
  expect_error( # expect no error (smoke test)
    lin_fit <- fit_glm(
      formula             = y ~ X1 + X2,
      p_family            = "gamma",
      data                = "random",
      iterations          = 10,
      warmup_iterations   = 5,
      chains              = 1, 
      seed                = 123
    ),
    NA
  )
})

# End-to-end smoke test of multiple formulas

test_that("gamma-multi", {
  expect_error( # expect no error (smoke test)
    lin_fit <- fit_glm(
      formula             = list(y1 ~ X1 + X2, y2 ~ X1 + X2 + y1),
      p_family            = "gamma",
      data                = "random",
      iterations          = 10,
      warmup_iterations   = 5,
      chains              = 1, 
      seed                = 123
    ),
    NA
  )
})
