# Unit tests for internal helper functions of dsPCAdrivers --------------------
#
# These tests exercise pure-R helpers and do NOT require DSLite, DataSHIELD or
# any disclosure-settings infrastructure. They access non-exported functions
# via the triple-colon operator.

# clean_vars_data -------------------------------------------------------------

test_that("clean_vars_data removes Inf, zero-variance and ID-like columns", {

  set.seed(1)
  n <- 20
  vars <- data.frame(
    age          = stats::rnorm(n),
    constant_num = rep(1, n),
    constant_chr = rep("a", n),
    id_chr       = paste0("id_", seq_len(n)),
    has_inf      = c(Inf, -Inf, stats::rnorm(n - 2)),
    sex          = sample(c("M", "F"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  cleaned <- dsPCAdrivers:::clean_vars_data(vars, na_threshold = 4)

  # Constant columns and ID-like character columns should be dropped.
  expect_false("constant_num" %in% colnames(cleaned))
  expect_false("constant_chr" %in% colnames(cleaned))
  expect_false("id_chr"       %in% colnames(cleaned))

  # Informative columns must survive.
  expect_true("age" %in% colnames(cleaned))
  expect_true("sex" %in% colnames(cleaned))

  # Inf values must be converted to NA.
  expect_true("has_inf" %in% colnames(cleaned))
  expect_true(any(is.na(cleaned$has_inf)))
  expect_false(any(is.infinite(cleaned$has_inf)))
})

test_that("clean_vars_data drops columns below the non-NA threshold", {

  n <- 20
  vars <- data.frame(
    full   = stats::rnorm(n),
    sparse = c(stats::rnorm(2), rep(NA, n - 2)),
    stringsAsFactors = FALSE
  )

  cleaned <- dsPCAdrivers:::clean_vars_data(vars, na_threshold = 4)

  expect_true("full" %in% colnames(cleaned))
  expect_false("sparse" %in% colnames(cleaned))
})

test_that("clean_vars_data keeps numeric columns even when all values are unique", {

  # Numeric columns whose values happen to be unique (e.g. continuous age)
  # must NOT be flagged as ID-like.
  n <- 20
  vars <- data.frame(
    unique_num = stats::runif(n),  # all unique numeric values
    sex        = sample(c("M", "F"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  cleaned <- dsPCAdrivers:::clean_vars_data(vars, na_threshold = 4)

  expect_true("unique_num" %in% colnames(cleaned))
})

test_that("clean_vars_data returns a 0-column data.frame when nothing survives", {

  n <- 20
  vars <- data.frame(
    const = rep(1, n),
    id    = paste0("x_", seq_len(n)),
    stringsAsFactors = FALSE
  )

  cleaned <- dsPCAdrivers:::clean_vars_data(vars, na_threshold = 4)

  expect_s3_class(cleaned, "data.frame")
  expect_equal(ncol(cleaned), 0)
})


# compute_single_association --------------------------------------------------

test_that("compute_single_association recovers a known numeric association", {

  set.seed(42)
  n  <- 50
  pc <- stats::rnorm(n)
  # Strong linear association with PC
  feature <- 2 * pc + stats::rnorm(n, sd = 0.2)

  pcs  <- matrix(pc, ncol = 1, dimnames = list(NULL, "PC1"))
  vars <- data.frame(x = feature)

  p_parametric <- dsPCAdrivers:::compute_single_association(
    feature_name = "x", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = TRUE,  threshold = 3
  )
  p_nonparam <- dsPCAdrivers:::compute_single_association(
    feature_name = "x", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = FALSE, threshold = 3
  )

  expect_true(is.finite(p_parametric))
  expect_true(is.finite(p_nonparam))
  expect_lt(p_parametric, 0.01)
  expect_lt(p_nonparam,   0.01)
})

test_that("compute_single_association uses ANOVA / Kruskal for categorical features", {

  set.seed(7)
  n <- 60
  group <- sample(c("a", "b", "c"), n, replace = TRUE)
  pc    <- ifelse(group == "a", 2, ifelse(group == "b", 0, -2)) +
           stats::rnorm(n, sd = 0.3)

  pcs  <- matrix(pc, ncol = 1, dimnames = list(NULL, "PC1"))
  vars <- data.frame(group = group, stringsAsFactors = FALSE)

  p_param <- dsPCAdrivers:::compute_single_association(
    feature_name = "group", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = TRUE,  threshold = 3
  )
  p_nonparam <- dsPCAdrivers:::compute_single_association(
    feature_name = "group", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = FALSE, threshold = 3
  )

  expect_lt(p_param,    1e-10)
  expect_lt(p_nonparam, 1e-10)
})

test_that("compute_single_association returns NA below the disclosure threshold", {

  pcs  <- matrix(stats::rnorm(3), ncol = 1, dimnames = list(NULL, "PC1"))
  vars <- data.frame(x = stats::rnorm(3))

  p <- dsPCAdrivers:::compute_single_association(
    feature_name = "x", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = TRUE, threshold = 5
  )

  expect_true(is.na(p))
})

test_that("compute_single_association returns NA for single-level or under-populated categorical features", {

  # All values identical -> only 1 level -> NA
  pcs  <- matrix(stats::rnorm(20), ncol = 1, dimnames = list(NULL, "PC1"))
  vars <- data.frame(g = rep("only_one", 20), stringsAsFactors = FALSE)

  p1 <- dsPCAdrivers:::compute_single_association(
    feature_name = "g", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = TRUE, threshold = 3
  )
  expect_true(is.na(p1))

  # 2 levels but one of them under threshold -> NA
  vars2 <- data.frame(
    g = c(rep("a", 19), "b"),
    stringsAsFactors = FALSE
  )
  p2 <- dsPCAdrivers:::compute_single_association(
    feature_name = "g", pc_name = "PC1",
    pcs = pcs, vars = vars2, parametric = TRUE, threshold = 3
  )
  expect_true(is.na(p2))
})

test_that("compute_single_association ignores rows with NA in either column", {

  set.seed(3)
  n <- 30
  pc <- stats::rnorm(n)
  x  <- pc + stats::rnorm(n, sd = 0.1)
  # Introduce NAs that must be filtered
  pc[1:3] <- NA
  x[28:30] <- NA

  pcs  <- matrix(pc, ncol = 1, dimnames = list(NULL, "PC1"))
  vars <- data.frame(x = x)

  p <- dsPCAdrivers:::compute_single_association(
    feature_name = "x", pc_name = "PC1",
    pcs = pcs, vars = vars, parametric = TRUE, threshold = 3
  )

  expect_true(is.finite(p))
  expect_lt(p, 0.01)
})
