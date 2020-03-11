test_that("Normalise checks work", {
  m <- matrix(data = sample(0:1, 10 * 20, replace = TRUE, prob = c((1 - 0.3), 0.3)), nrow = 10, ncol = 20)
  mc <- mcount(M = m, six_node = TRUE, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)
  expect_equal(object = sizeclass_check(target_motif_distribution = mc$normalise_sizeclass), expected = TRUE)
  expect_equal(object = sizeclass_check(target_motif_distribution = mc$normalise_sizeclass[1:17]), expected = TRUE)
  expect_equal(object = levelsize_check(target_motif_distribution = mc$normalise_levelsize), expected = TRUE)
  expect_equal(object = levelsize_check(target_motif_distribution = mc$normalise_levelsize[1:17]), expected = TRUE)
})

