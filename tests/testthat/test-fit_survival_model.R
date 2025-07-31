test_that("Testing works", {
  expect_equal(2 * 2, 4)
})

test_that("Testing works", {
  cov_matrix <- matrix(c(
    2, 1,
    1, 2 
  ), nrow = 2, byrow = TRUE)
  
  n <- 200             # observations
  
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  
  # First response: depends on x1 and x2
  y1 <- 1 + 2 * X1 - 1.5 * X2 + rnorm(n, sd = 1)
  
  # Second response: depends on x1, x2, and y1
  y2 <- -0.5 + 0.8 * X1 + 0.4 * X2 + 1.2 * y1 + rnorm(n, sd = 1)
  
  synthetic_df <- data.frame(X1, X2, y1, y2)
  
  lin_fit <- LinksMixtureModeling::fit_glm(
    formula             = y1 ~ X1 + X2,
    p_family            = "linear",
    data                = "random",
    components          = c("linear", "linear"),
    iterations          = 10,
    warmup_iterations   = 5,
    chains              = 2, 
    seed                = 123
  )
})