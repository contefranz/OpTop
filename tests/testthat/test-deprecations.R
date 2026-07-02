# The legacy stability/matching helpers are deprecated as of 0.9.8 and
# scheduled for removal; each exported one must signal a lifecycle deprecation
# warning on entry (before argument validation stops the call).

test_that("legacy stability helpers warn about their deprecation", {
  expect_warning(try(topic_stability(), silent = TRUE),
                 class = "lifecycle_warning_deprecated")
  expect_warning(try(agg_topic_stability(), silent = TRUE),
                 class = "lifecycle_warning_deprecated")
  expect_warning(try(agg_document_stability(), silent = TRUE),
                 class = "lifecycle_warning_deprecated")
  expect_warning(get_topic_models(envir = new.env()),
                 class = "lifecycle_warning_deprecated")
})
