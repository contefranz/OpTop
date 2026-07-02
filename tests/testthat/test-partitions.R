test_that("optop_make_partition returns document lengths and a labelled mask", {
  fx <- optop_test_fixture()
  part <- fx$partition

  expect_named(part, c("rare_mask", "L"))
  expect_identical(dim(part$rare_mask), c(fx$J, fx$W))
  expect_type(part$rare_mask, "logical")
  expect_identical(colnames(part$rare_mask), colnames(fx$counts))
  expect_identical(rownames(part$rare_mask), rownames(fx$counts))
  expect_equal(unname(part$L), unname(Matrix::rowSums(fx$dtm)))

  # the fixture must exercise both branches of the harmonized support
  expect_true(any(part$rare_mask))
  expect_true(any(!part$rare_mask))
})

test_that("rare_mask equals the running min of Theta Phi against tau", {
  fx <- optop_test_fixture()

  L <- Matrix::rowSums(fx$dtm)
  tau <- 5 / pmax(L, 1)
  i_min <- Reduce(pmin, lapply(fx$models, function(m) {
    tp <- ref_theta_phi(m)
    tp$theta %*% tp$phi
  }))
  expected <- i_min < tau  # recycles tau_j down each column (row j)

  expect_equal(unname(fx$partition$rare_mask), unname(expected))
})

test_that("partition is invariant to the block size", {
  fx <- optop_test_fixture()
  part_small <- optop_make_partition(fx$models, fx$dtm, c = 5, block = 7)
  expect_equal(part_small$rare_mask, fx$partition$rare_mask)
  expect_equal(part_small$L, fx$partition$L)
})

test_that("larger c marks more words as rare", {
  fx <- optop_test_fixture()
  part_loose <- optop_make_partition(fx$models, fx$dtm, c = 20)
  expect_true(all(fx$partition$rare_mask <= part_loose$rare_mask))
  expect_gt(sum(part_loose$rare_mask), sum(fx$partition$rare_mask))
})
