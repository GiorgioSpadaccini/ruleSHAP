test_that("Friedman generating function runs smoothly", {
  data=expect_no_error(gendata.friedman1(1e3,p=10,sd=1,binary=F))
  expect_equal(ncol(data), 11)
})

test_that("Friedman generating function can generate binary outcome", {
  data=expect_no_error(gendata.friedman1(1e3,p=10,sd=1,binary=T))
  expect_equal(all(unique(data$y) %in% 0:1),T)
})
