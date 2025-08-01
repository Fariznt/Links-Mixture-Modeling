# This file will test a likely point of edge case failure, priors string processing
# Most of the other tests will likely be end-to-end given that almost all 
# of the compile/runtime is from the stan model (so its not worth the effort 
# skipping the before/after functions), its the likeliest point of failure, 
# and its most important to functionality

# unit tests make sense for this file

# While this is the set of tests for which unit tests make the most sense, we
# could also do end to end tests with really small iterations and simply check
# that the code runs with no error (smoke tests). 

# delete everything above once read^

# Tests for custom prior functionality

test_that("Testing works", {
  expect_equal(2 * 2, 4)
})

test_that("name here", {
  # no need for expect_equal for very basic tests
  # function call only. if an error happens, the test fails
  # these are appropriate here because if stan compiled and there were values
  # in the posterior, the priors almost certainly worked
})
