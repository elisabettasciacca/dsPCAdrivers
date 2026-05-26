# End-to-end tests for plotDriversDS ------------------------------------------
#
# Two layers:
#   1. Input validation (pure R, no DSLite required).
#   2. End-to-end execution via DSLite + dsBase, exercising the full code path
#      including .getMinObsSetting().

# Input validation ------------------------------------------------------------

test_that("plotDriversDS rejects invalid name arguments", {

  expect_error(
    plotDriversDS(pcs.name = NULL, vars.name = "vars"),
    "pcs.name must be a single character string"
  )
  expect_error(
    plotDriversDS(pcs.name = "pcs", vars.name = NULL),
    "vars.name must be a single character string"
  )
  expect_error(
    plotDriversDS(pcs.name = c("a", "b"), vars.name = "vars"),
    "pcs.name must be a single character string"
  )
  expect_error(
    plotDriversDS(pcs.name = 1, vars.name = "vars"),
    "pcs.name must be a single character string"
  )
})

test_that("plotDriversDS rejects unknown p_adj methods", {
  expect_error(
    plotDriversDS(pcs.name = "pcs", vars.name = "vars", p_adj = "not_a_method"),
    "'p_adj' must be one of"
  )
})

test_that("plotDriversDS errors when referenced objects are missing", {

  # Reference a name that doesn't exist in the parent frame
  rm(list = ls())  # ensure no stray "pcs"/"vars" objects
  expect_error(
    plotDriversDS(pcs.name = "nonexistent_pcs", vars.name = "nonexistent_vars"),
    "not found on server"
  )
})

test_that("plotDriversDS errors when objects have incompatible shapes", {

  pcs  <- matrix(stats::rnorm(20), ncol = 2)
  vars <- data.frame(age = stats::rnorm(10))  # different nrow

  expect_error(
    plotDriversDS(pcs.name = "pcs", vars.name = "vars"),
    "must have same number of rows"
  )

  pcs_wrong  <- "not_a_matrix"
  vars_ok    <- data.frame(age = stats::rnorm(10))
  expect_error(
    plotDriversDS(pcs.name = "pcs_wrong", vars.name = "vars_ok"),
    "must be a matrix or data frame"
  )
})


# End-to-end via DSLite -------------------------------------------------------

test_that("plotDriversDS returns correctly structured output via DSLite", {

  setup <- setup_single_site_dslite(n = 50, seed = 42L)
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  cally <- call("plotDriversDS",
                pcs.name          = "pcs",
                vars.name         = "vars",
                parametric        = TRUE,
                n_pc              = 5L,
                na_drop_threshold = 4L,
                p_adj             = NULL,
                verbose           = FALSE)

  out <- DSI::datashield.aggregate(conns, cally)
  expect_named(out, "site1")

  res <- out$site1
  expect_type(res, "list")
  expect_named(res,
               c("results", "pc_names", "var_names",
                 "n_observations", "parameters"),
               ignore.order = TRUE)

  expect_s3_class(res$results, "data.frame")
  expect_named(res$results,
               c("Feature", "PC", "pvalue", "Association"),
               ignore.order = TRUE)

  # 4 vars (age, sex, batch, bmi) x 5 PCs = 20 rows
  expect_equal(nrow(res$results), length(res$var_names) * length(res$pc_names))
  expect_equal(length(res$pc_names), 5L)
  expect_equal(res$n_observations, 50L)

  # All p values must be in [0, 1] or NA
  expect_true(all(is.na(res$results$pvalue) |
                  (res$results$pvalue >= 0 & res$results$pvalue <= 1)))

  # Association = -log10(p) must be non-negative or NA
  expect_true(all(is.na(res$results$Association) |
                  res$results$Association >= 0))
})

test_that("plotDriversDS detects the planted PC1 ~ age signal", {

  setup <- setup_single_site_dslite(n = 80, seed = 99L)
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  cally <- call("plotDriversDS",
                pcs.name = "pcs", vars.name = "vars",
                parametric = TRUE, n_pc = 5L, na_drop_threshold = 4L)

  res <- DSI::datashield.aggregate(conns, cally)$site1$results

  # The helper injects 0.5 * scale(age) into PC1 -> expect a small p value
  age_pc1 <- res$pvalue[res$Feature == "age" & res$PC == "PC1"]
  expect_true(is.finite(age_pc1))
  expect_lt(age_pc1, 0.05)
})

test_that("plotDriversDS truncates n_pc when it exceeds available PCs", {

  setup <- setup_single_site_dslite(n = 40, seed = 1L)
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  cally <- call("plotDriversDS",
                pcs.name = "pcs", vars.name = "vars",
                n_pc = 50L, na_drop_threshold = 4L)

  res <- DSI::datashield.aggregate(conns, cally)$site1

  # PCs in the helper are 5; n_pc was 50 -> must be capped at 5
  expect_equal(length(res$pc_names), 5L)
})

test_that("plotDriversDS applies p_adj when requested", {

  setup <- setup_single_site_dslite(n = 50, seed = 7L)
  conns <- setup$conns
  withr::defer(DSI::datashield.logout(conns))

  cally_raw <- call("plotDriversDS",
                    pcs.name = "pcs", vars.name = "vars",
                    n_pc = 5L, na_drop_threshold = 4L, p_adj = NULL)
  cally_bh  <- call("plotDriversDS",
                    pcs.name = "pcs", vars.name = "vars",
                    n_pc = 5L, na_drop_threshold = 4L, p_adj = "BH")

  res_raw <- DSI::datashield.aggregate(conns, cally_raw)$site1$results
  res_bh  <- DSI::datashield.aggregate(conns, cally_bh)$site1$results

  # BH-adjusted p values must be >= raw p values (no NAs in this scenario)
  ok <- !is.na(res_raw$pvalue) & !is.na(res_bh$pvalue)
  expect_true(all(res_bh$pvalue[ok] >= res_raw$pvalue[ok] - 1e-12))
})

test_that("plotDriversDS gracefully handles NA-heavy variables", {

  skip_if_no_dslite()

  set.seed(13)
  n <- 60
  pcs <- matrix(stats::rnorm(n * 5), ncol = 5,
                dimnames = list(NULL, paste0("comp.", 1:5)))
  vars <- data.frame(
    age          = stats::rnorm(n),
    mostly_na    = c(stats::rnorm(2), rep(NA, n - 2)),  # below threshold
    constant     = rep(1, n),
    sex          = sample(c("M", "F"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )

  dslite_server <- DSLite::newDSLiteServer(
    tables = list(
      pcs_table  = as.data.frame(pcs),
      vars_table = vars
    )
  )
  dslite_server$config(DSLite::defaultDSConfiguration(include = "dsBase"))
  dslite_server$aggregateMethod("plotDriversDS", dsPCAdrivers::plotDriversDS)

  options(datashield.privacyControlLevel = "permissive",
          nfilter.tab = 3, nfilter.subset = 3)

  builder <- DSI::newDSLoginBuilder()
  builder$append(server = "site1", url = "dslite_server",
                 table = "pcs_table", driver = "DSLiteDriver")
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)
  withr::defer(DSI::datashield.logout(conns))

  DSI::datashield.assign.table(conns, "pcs",  c(site1 = "pcs_table"))
  DSI::datashield.assign.table(conns, "vars", c(site1 = "vars_table"))

  cally <- call("plotDriversDS", pcs.name = "pcs", vars.name = "vars",
                n_pc = 5L, na_drop_threshold = 4L)
  res <- DSI::datashield.aggregate(conns, cally)$site1

  # mostly_na and constant must have been filtered out by clean_vars_data
  expect_false("mostly_na" %in% res$var_names)
  expect_false("constant"  %in% res$var_names)
  expect_true("age" %in% res$var_names)
  expect_true("sex" %in% res$var_names)
})
