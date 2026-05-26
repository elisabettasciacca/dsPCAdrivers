# Test helpers for dsPCAdrivers ------------------------------------------------
#
# Shared utilities used across all test files. `testthat` sources this file
# automatically before running any test.
#
# The helpers below cover two needs:
#   1. Deterministic in-memory data generation for pure-R unit tests.
#   2. A minimal DSLite session that registers `plotDriversDS` as an aggregate
#      method, so the end-to-end behaviour can be exercised without a real
#      DataSHIELD server.

# Skip helpers ----------------------------------------------------------------

skip_if_no_dslite <- function() {
  testthat::skip_if_not_installed("DSLite")
  testthat::skip_if_not_installed("DSI")
  testthat::skip_if_not_installed("dsBase")
}

# Synthetic data --------------------------------------------------------------

#' Build a small reproducible PCs/vars dataset for unit tests.
#'
#' @param n Integer. Number of samples.
#' @param n_pc Integer. Number of PCs to simulate.
#' @param seed Integer. Random seed.
#' @return A list with elements `pcs` (matrix) and `vars` (data.frame).
make_test_data <- function(n = 50, n_pc = 5, seed = 1L) {
  set.seed(seed)

  pcs <- matrix(
    stats::rnorm(n * n_pc),
    nrow = n,
    ncol = n_pc,
    dimnames = list(NULL, paste0("comp.", seq_len(n_pc)))
  )

  vars <- data.frame(
    age   = stats::rnorm(n, mean = 50, sd = 10),
    sex   = sample(c("M", "F"), n, replace = TRUE),
    batch = factor(sample(seq_len(3), n, replace = TRUE)),
    bmi   = stats::rnorm(n, mean = 25, sd = 5),
    stringsAsFactors = FALSE
  )

  # Inject a true association between PC1 and `age` so the test can verify
  # that the function detects a known signal (p value below 0.05).
  pcs[, 1] <- pcs[, 1] + 0.5 * scale(vars$age)[, 1]

  list(pcs = pcs, vars = vars)
}

# DSLite session --------------------------------------------------------------

#' Create a single-site DSLite session with `plotDriversDS` registered.
#'
#' The PCs and vars matrices are pre-assigned on the server under the symbols
#' "pcs" and "vars". The caller is responsible for `datashield.logout()` —
#' use `withr::defer(datashield.logout(conns))` inside the test.
#'
#' The DSLite server object is assigned to .GlobalEnv under the name
#' "dslite_server" because the DSLiteDriver looks the server up by name from
#' the calling frames. `withr::defer()` removes it after the test.
#'
#' @return A list with `conns` (DSConnection list) and `data` (raw test data).
setup_single_site_dslite <- function(n = 50, seed = 42L) {
  skip_if_no_dslite()

  data <- make_test_data(n = n, n_pc = 5, seed = seed)

  dslite_server <- DSLite::newDSLiteServer(
    tables = list(
      site1_pcs  = as.data.frame(data$pcs),
      site1_vars = data$vars
    )
  )
  dslite_server$config(DSLite::defaultDSConfiguration(include = "dsBase"))
  dslite_server$aggregateMethod("plotDriversDS", dsPCAdrivers::plotDriversDS)

  assign("dslite_server", dslite_server, envir = globalenv())
  withr::defer(
    suppressWarnings(rm("dslite_server", envir = globalenv())),
    envir = parent.frame()
  )

  options(
    datashield.privacyControlLevel = "permissive",
    nfilter.tab     = 3,
    nfilter.subset  = 3,
    nfilter.glm     = 1.0,
    nfilter.string  = 80,
    nfilter.kNN     = 3
  )

  builder <- DSI::newDSLoginBuilder()
  builder$append(
    server = "site1",
    url    = "dslite_server",
    table  = "site1_pcs",
    driver = "DSLiteDriver"
  )
  conns <- DSI::datashield.login(logins = builder$build(), assign = FALSE)

  DSI::datashield.assign.table(conns, "pcs",  c(site1 = "site1_pcs"))
  DSI::datashield.assign.table(conns, "vars", c(site1 = "site1_vars"))

  list(conns = conns, data = data)
}
