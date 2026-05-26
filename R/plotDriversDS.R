#' @title Server-side function to compute PC-variable associations
#'
#' @description Computes associations between principal components and variables
#' without exposing individual-level data. Enforces DataSHIELD disclosure settings.
#'
#' @param pcs.name Character string specifying the name of the matrix/dataframe containing PC scores
#' @param vars.name Character string specifying the name of the dataframe containing variables
#' @param parametric Logical. Use parametric tests? Default TRUE.
#' @param n_pc Integer. Number of PCs to include. Default 5.
#' @param na_drop_threshold Integer. Minimum non-NA values required. Default 4.
#' @param p_adj Character. P value adjustment method. Default NULL.
#' @param verbose Logical. Print messages? Default FALSE.
#'
#' @return A list containing results, variable names, and parameters.
#' @export
plotDriversDS <- function(pcs.name = NULL,
                          vars.name = NULL,
                          parametric = TRUE,
                          n_pc = 5L,
                          na_drop_threshold = 4,
                          p_adj = NULL,
                          verbose = FALSE) {

  # Validate argument types (function only accepts object names as strings)
  if (!is.character(pcs.name) || length(pcs.name) != 1) {
    stop("pcs.name must be a single character string specifying a server object name", call. = FALSE)
  }

  if (!is.character(vars.name) || length(vars.name) != 1) {
    stop("vars.name must be a single character string specifying a server object name", call. = FALSE)
  }

  valid_p_adj_methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
  if (!is.null(p_adj) && !p_adj %in% valid_p_adj_methods) {
    stop("'p_adj' must be one of: ", paste(valid_p_adj_methods, collapse = ", "), call. = FALSE)
  }

  # Resolve object names in the server environment
  pcs <- tryCatch(
    eval(parse(text = pcs.name), envir = parent.frame()),
    error = function(e) stop("Object '", pcs.name, "' not found on server", call. = FALSE)
  )
  vars <- tryCatch(
    eval(parse(text = vars.name), envir = parent.frame()),
    error = function(e) stop("Object '", vars.name, "' not found on server", call. = FALSE)
  )

  if (is.null(pcs))  stop("Object '", pcs.name,  "' is NULL on server", call. = FALSE)
  if (is.null(vars)) stop("Object '", vars.name, "' is NULL on server", call. = FALSE)

  if (!is.matrix(pcs) && !is.data.frame(pcs)) {
    stop("Object '", pcs.name, "' must be a matrix or data frame", call. = FALSE)
  }

  if (!is.data.frame(vars) && !is.matrix(vars)) {
    stop("Object '", vars.name, "' must be a data frame or matrix", call. = FALSE)
  }

  if (nrow(pcs) != nrow(vars)) {
    stop("Objects '", pcs.name, "' and '", vars.name, "' must have same number of rows", call. = FALSE)
  }

  # truncate to n_pc if required (optional)
  if (n_pc > ncol(pcs)) {
    if (verbose) message("n_pc exceeds available PCs, using all ", ncol(pcs), " PCs")
    n_pc <- ncol(pcs)
  }

  pcs <- pcs[, seq_len(n_pc), drop = FALSE]
  # Rename columns from comp.1, comp.2... to PC1, PC2...
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

  # Clean variables data
  vars_clean <- dsPCAdrivers:::clean_vars_data(
    vars = vars,
    na_threshold = na_drop_threshold,
    verbose = verbose
  )

  if (ncol(vars_clean) == 0) {
    stop("No valid variables remaining after filtering", call. = FALSE)
  }

  # Compute associations
  if (verbose) message("Computing associations for ", ncol(vars_clean),
                       " variables and ", n_pc, " PCs")

  pc_names <- colnames(pcs)
  var_names <- colnames(vars_clean)

  # Retrieve disclosure threshold once, then pass it to each association test
  threshold <- dsPCAdrivers:::.getMinObsSetting()

  # Build a data frame with all variable × PC combinations
  combinations <- expand.grid(
    feature = factor(var_names, levels = var_names),
    pc = factor(pc_names, levels = pc_names),
    stringsAsFactors = FALSE
  )

  combinations$feature <- as.character(combinations$feature)
  combinations$pc <- as.character(combinations$pc)

  # Compute associations for each variable-PC pair via mapply
  pvals <- mapply(
    FUN = dsPCAdrivers:::compute_single_association,
    feature_name = combinations$feature,
    pc_name = combinations$pc,
    MoreArgs = list(
      pcs       = pcs,
      vars      = vars_clean,
      parametric = parametric,
      threshold  = threshold,
      verbose   = verbose
    ),
    SIMPLIFY = TRUE
  )

  # Build results
  results <- data.frame(
    Feature = combinations$feature,
    PC = combinations$pc,
    pvalue = pvals,
    stringsAsFactors = FALSE
  )

  results$Feature <- factor(results$Feature, levels = var_names)
  results$PC <- factor(results$PC, levels = pc_names)

  # multiple testing correction
  if (!is.null(p_adj)) {
    if (verbose) message("Applying ", p_adj, " correction")
    results$pvalue <- p.adjust(results$pvalue, method = p_adj)
  }

  results$Association <- -log10(results$pvalue)
  results$Feature <- as.character(results$Feature)
  results$PC <- as.character(results$PC)

  # Prepare output (only aggregated data - no individual data)
  output <- list(
    results = results,
    pc_names = pc_names,
    var_names = var_names,
    n_observations = nrow(pcs),
    parameters = list(
      parametric = parametric,
      p_adj = p_adj,
      na_drop_threshold = na_drop_threshold
    )
  )

  return(output)
}

# Helper functions -------------------------------------------------------------
#' @title Get minimum observations threshold from DataSHIELD disclosure settings
#' @description Queries the DataSHIELD server disclosure settings and returns
#'   the minimum number of observations allowed, computed as the maximum between
#'   \code{nfilter.tab} and \code{nfilter.subset}. Falls back to 3 if either
#'   setting is missing or \code{NA}.
#' @return A single numeric value representing the minimum observations threshold.
#'
#' @importFrom dsBase listDisclosureSettingsDS
#' @keywords internal
.getMinObsSetting <- function() {
  settings <- dsBase::listDisclosureSettingsDS()

  n_tab <- as.numeric(settings$nfilter.tab)
  n_subset <- as.numeric(settings$nfilter.subset)

  if (is.null(n_tab) || is.na(n_tab)) n_tab <- 3
  if (is.null(n_subset) || is.na(n_subset)) n_subset <- 3

  return(max(n_tab, n_subset))
}

#' Clean and validate a variables data frame
#'
#' Prepares a data frame of clinical or technical variables for association
#' testing by removing columns that are uninformative or unsuitable: infinite
#' values are coerced to \code{NA}, then zero-variance, data-sparse, and
#' ID-like columns are dropped.
#'
#' @param vars A data frame of variables (samples as rows, variables as
#'   columns).
#' @param na_threshold Integer. Minimum number of non-\code{NA} values required
#'   for a column to be retained.
#' @param verbose Logical. If \code{TRUE}, messages are printed when columns
#'   are removed. Default is \code{FALSE}.
#'
#' @return A data frame with unsuitable columns removed.
#'
#' @keywords internal
clean_vars_data <- function(vars, na_threshold, verbose = FALSE) {

  # Convert Inf/-Inf to NA using lapply
  vars[] <- lapply(vars, function(x) {
    if (is.numeric(x)) {
      x[is.infinite(x)] <- NA
    }
    return(x)
  })

  has_variance <- vapply(vars, function(x) {
    length(unique(x[!is.na(x)])) > 1
  }, logical(1))

  if (any(!has_variance) && verbose) {
    message("Removing zero-variance columns")
  }

  vars <- vars[, has_variance, drop = FALSE]

  sufficient_data <- colSums(!is.na(vars)) >= na_threshold

  if (any(!sufficient_data) && verbose) {
    message("Removing columns with insufficient non-NA values")
  }

  vars <- vars[, sufficient_data, drop = FALSE]

  is_numeric <- vapply(vars, is.numeric, logical(1))
  all_unique <- vapply(vars, function(x) {
    n_unique <- length(unique(x[!is.na(x)]))
    n_nonmissing <- sum(!is.na(x))
    # Potential ID columns: character/factor variables where every
    # non-missing value is unique (n_unique == n_nonmissing).
    return((n_unique == n_nonmissing) && (n_nonmissing > 0))
  }, logical(1))

  is_id_column <- all_unique & !is_numeric

  if (any(is_id_column) && verbose) {
    message("Removing ", sum(is_id_column), " ID-like columns containing unique values: ",
            paste(colnames(vars)[is_id_column], collapse = ", "))
  }

  vars <- vars[, !is_id_column, drop = FALSE]

  return(vars)
}



#' @title Compute association between a single variable and a single PC
#'
#' @description Tests the association between one variable and one principal
#'   component score, selecting the appropriate statistical test based on the
#'   variable type. For numeric variables, uses Pearson correlation
#'   (parametric) or Spearman correlation (non-parametric). For categorical
#'   variables, uses one-way ANOVA (parametric) or Kruskal-Wallis
#'   (non-parametric). Returns \code{NA} if the number of complete cases falls
#'   below the DataSHIELD disclosure threshold, or if a categorical variable
#'   has fewer than two levels or any level with too few observations.
#'
#' @param feature_name Character string. Name of the variable column in
#'   \code{vars}.
#' @param pc_name Character string. Name of the PC column in \code{pcs}.
#' @param pcs Matrix or data frame of PC scores (samples × PCs).
#' @param vars Data frame of variables (samples × features).
#' @param parametric Logical. Use parametric tests? Default \code{TRUE}.
#' @param threshold Integer. Minimum number of complete cases required, as
#'   returned by \code{.getMinObsSetting()}. Passed in from the parent call to
#'   avoid querying the disclosure settings once per pair.
#' @param verbose Logical. Print diagnostic messages? Default \code{FALSE}.
#'
#' @return A single numeric p value, or \code{NA_real_} if the association
#'   could not be computed.
#'
#' @keywords internal
compute_single_association <- function(feature_name,
                                       pc_name,
                                       pcs,
                                       vars,
                                       parametric,
                                       threshold,
                                       verbose = FALSE) {

  # Extract data and filter to complete cases
  pc_values <- pcs[, pc_name]
  feature_values <- vars[[feature_name]]

  complete_cases <- !is.na(pc_values) & !is.na(feature_values)

  if (sum(complete_cases) < threshold) {
    return(NA_real_)
  }

  pc_values <- pc_values[complete_cases]
  feature_values <- feature_values[complete_cases]

  # Select and run the appropriate statistical test
  if (is.numeric(feature_values)) {
    method <- if (parametric) "pearson" else "spearman"
    test_result <- cor.test(pc_values, feature_values, method = method)
    return(test_result$p.value)
  } else {
    feature_values <- as.factor(feature_values)

    # Check levels for disclosure risk
    level_counts <- table(feature_values)
    if (length(levels(feature_values)) < 2 || any(level_counts < threshold)) {
      return(NA_real_)
    }

    if (parametric) {
      fit <- lm(pc_values ~ feature_values)
      return(anova(fit)[1, "Pr(>F)"])
    } else {
      test_result <- kruskal.test(pc_values ~ feature_values)
      return(test_result$p.value)
    }
  }
}
