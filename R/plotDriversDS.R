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

  # DataSHIELD security checks ---------------------------------------------------

  if (!is.character(pcs.name) || length(pcs.name) != 1) {
    stop("pcs.name must be a single character string specifying server object name", call. = FALSE)
  }

  if (!is.character(vars.name) || length(vars.name) != 1) {
    stop("vars.name must be a single character string specifying server object name", call. = FALSE)
  }

  pcs <- eval(parse(text = pcs.name), envir = parent.frame())
  vars <- eval(parse(text = vars.name), envir = parent.frame())

  if (is.null(pcs)) stop("Object '", pcs.name, "' not found on server", call. = FALSE)
  if (is.null(vars)) stop("Object '", vars.name, "' not found on server", call. = FALSE)

  if (!is.matrix(pcs) && !is.data.frame(pcs)) {
    stop("Object '", pcs.name, "' must be a matrix or data frame", call. = FALSE)
  }

  if (!is.data.frame(vars) && !is.matrix(vars)) {
    stop("Object '", vars.name, "' must be a data frame or matrix", call. = FALSE)
  }

  if (nrow(pcs) != nrow(vars)) {
    stop("Objects '", pcs.name, "' and '", vars.name, "' must have same number of rows", call. = FALSE)
  }

  if (n_pc > ncol(pcs)) {
    if (verbose) message("n_pc exceeds available PCs, using all ", ncol(pcs), " PCs")
    n_pc <- ncol(pcs)
  }

  pcs <- pcs[, seq_len(n_pc), drop = FALSE]
  # Rename columns from comp.1, comp.2... to PC1, PC2...
  colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))

  # Clean variables data ---------------------------------------------------------

  vars_clean <- clean_vars_data(
    vars = vars,
    na_threshold = na_drop_threshold,
    verbose = verbose
  )

  if (ncol(vars_clean) == 0) {
    stop("No valid variables remaining after filtering", call. = FALSE)
  }

  # Compute associations ---------------------------------------------------------

  if (verbose) message("Computing associations for ", ncol(vars_clean),
                       " variables and ", n_pc, " PCs")

  pc_names <- colnames(pcs)
  var_names <- colnames(vars_clean)

  combinations <- expand.grid(
    feature = factor(var_names, levels = var_names),
    pc = factor(pc_names, levels = pc_names),
    stringsAsFactors = FALSE
  )

  combinations$feature <- as.character(combinations$feature)
  combinations$pc <- as.character(combinations$pc)

  pvals <- mapply(
    FUN = compute_single_association,
    feature_name = combinations$feature,
    pc_name = combinations$pc,
    MoreArgs = list(
      pcs = pcs,
      vars = vars_clean,
      parametric = parametric,
      verbose = verbose
    ),
    SIMPLIFY = TRUE
  )

  # Build results ----------------------------------------------------------------

  results <- data.frame(
    Feature = combinations$feature,
    PC = combinations$pc,
    pvalue = pvals,
    stringsAsFactors = FALSE
  )

  results$Feature <- factor(results$Feature, levels = var_names)
  results$PC <- factor(results$PC, levels = pc_names)

  if (!is.null(p_adj)) {
    if (verbose) message("Applying ", p_adj, " correction")
    results$pvalue <- p.adjust(results$pvalue, method = p_adj)
  }

  results$Association <- -log10(results$pvalue)
  results$Feature <- as.character(results$Feature)
  results$PC <- as.character(results$PC)

  # Prepare output ---------------------------------------------------------------

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

#' @importFrom dsBase listDisclosureSettingsDS
.getMinObsSetting <- function() {
  settings <- dsBase::listDisclosureSettingsDS()

  n_tab <- as.numeric(settings$nfilter.tab)
  n_subset <- as.numeric(settings$nfilter.subset)

  if (is.null(n_tab) || is.na(n_tab)) n_tab <- 3
  if (is.null(n_subset) || is.na(n_subset)) n_subset <- 3

  return(max(n_tab, n_subset))
}

clean_vars_data <- function(vars, na_threshold, verbose = FALSE) {

  # Replace infinite values with NA using lapply
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

compute_single_association <- function(feature_name,
                                       pc_name,
                                       pcs,
                                       vars,
                                       parametric,
                                       verbose = FALSE) {

  pc_values <- pcs[, pc_name]
  feature_values <- vars[[feature_name]]

  complete_cases <- !is.na(pc_values) & !is.na(feature_values)

  threshold <- .getMinObsSetting()

  if (sum(complete_cases) < threshold) {
    return(NA_real_)
  }

  pc_values <- pc_values[complete_cases]
  feature_values <- feature_values[complete_cases]

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
