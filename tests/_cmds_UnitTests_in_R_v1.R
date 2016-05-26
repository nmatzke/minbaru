
# UNIT TESTS IN R -- THE BARE MINIMUM
# 
# http://www.johnmyleswhite.com/notebook/2010/08/17/unit-testing-in-r-the-bare-minimum/
# 





# COMMANDS:

cd /drives/GDrive/__github/minbaru/
Rscript tests/run_tests.R





library(testthat)

# Passing test:
test_that(desc="trigonometric functions match identities", code={
  expect_equal(sin(pi / 4), 1 / sqrt(2))
  expect_equal(cos(pi / 4), 1 / sqrt(2))
  expect_equal(tan(pi / 4), 1)
})

# Failing test:
test_that(desc="trigonometric functions match identities", code={
  expect_equal(sin(pi / 4), 1)
})


