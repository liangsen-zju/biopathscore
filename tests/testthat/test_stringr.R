# ********************************************
# * Name: sen@uga.edu
# * Date: Thu Jun 15 10:34:03 2017
# * Detail: test
# *
# ********************************************
library(testthat)
library(stringr)
context("String length")

test_that("str_length is number of characters", {
  expect_equal(str_length('a'), 1)
  expect_equal(str_length('ab'), 2)
  expect_equal(str_length('abc'), 3)
})

