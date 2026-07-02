test_that("optop_make_baseline returns the corpus term distribution", {
  fx <- optop_test_fixture()
  base <- fx$baseline

  expect_named(base, "pi_glob")
  expect_length(base$pi_glob, fx$W)
  expect_identical(names(base$pi_glob), colnames(fx$counts))
  expect_equal(sum(base$pi_glob), 1)
  expect_equal(unname(base$pi_glob),
               unname(colSums(fx$counts) / sum(fx$counts)))
})
