# This file will test a likely point of edge case failure, priors string processing
# Most of the other tests will likely be end-to-end given that almost all 
# of the compile/runtime is from the stan model (so its not worth the effort 
# skipping the before/after functions), its the likeliest point of failure, 
# and its most important to functionality

test_that("Testing works", {
  expect_equal(2 * 2, 4)
})
