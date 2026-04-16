test_that("plotDriversDS works with basic input", {
  # Simulate data
  pcs <- matrix(rnorm(100 * 5), ncol = 5)
  colnames(pcs) <- paste0("PC", 1:5)

  vars <- data.frame(
    age = rnorm(100, 50, 10),
    sex = sample(c("M", "F"), 100, replace = TRUE),
    batch = factor(sample(1:3, 100, replace = TRUE))
  )

  # Run function
  result <- plotDriversDS(
    pcs.name = "pcs",
    vars.name = "vars",
    parametric = TRUE,
    n_pc = 5,
    na_drop_threshold = 4
  )

  # Check output structure
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_true("pc_names" %in% names(result))
  expect_true("var_names" %in% names(result))
  expect_equal(length(result$pc_names), 5)
})
