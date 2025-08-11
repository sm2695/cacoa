#' @keywords internal
estimateGraphVarianceSignificance <- function(adj.mat, signal, n.permutations=5000) {
  if (!is.numeric(signal)) {
    signal <- as.factor(signal)
    signal[is.na(signal)] <- table(signal) %>% which.max() %>% names()
    comp.op <- "!="
  } else {
    signal[is.na(signal)] <- median(signal, na.rm=TRUE)
    comp.op <- "-"
  }

  obs.var <- signal %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  perm.vars <- sapply(1:n.permutations, function(i) {
    sample(signal) %>% outer(., ., comp.op) %>% {. * . * adj.mat} %>% sum()
  })

  pvalue <- (sum(perm.vars <= obs.var) + 1) / (n.permutations + 1)
  pr2 <- 1 - obs.var / median(perm.vars)
  return(list(pvalue=pvalue, pr2=pr2))
}

#' @keywords internal
adjacencyMatrixFromPaiwiseDists <- function(p.dists, trim=0.05, k=NULL) {
  adj.mat <- p.dists %>% pmin(quantile(., 1 - trim)) %>%
    {pmax(0, . - quantile(., trim))} %>% {1 - . / max(.)} %>%
    matrix(ncol=ncol(p.dists))
  diag(adj.mat) <- 0

  if (!is.null(k) && (k < ncol(adj.mat))) { # remove edges for all but k nearest neighbors
    adj.mat %<>% apply(1, function(r) ifelse(r < sort(r, decreasing=TRUE)[k], 0, r)) %>% {(. + t(.)) / 2}
  }

  dimnames(adj.mat) <- dimnames(p.dists)

  return(adj.mat)
}

#' @keywords internal
estimateUMAPOnDistances <- function(p.dists, n.neighbors=15, verbose=FALSE, ...) {
  set.seed(42)
  idx <- (1:ncol(p.dists)) %>% lapply(function(i) head(order(p.dists[i,]), n.neighbors))
  dists <- (1:ncol(p.dists)) %>% lapply(function(i) p.dists[i, idx[[i]]])

  idx <- do.call(rbind, idx)
  dists <- do.call(rbind, dists)

  umap <- uwot::umap(data.frame(x=rep(0, nrow(dists))), nn_method=list(idx=idx, dist=dists), verbose=verbose, ...)

  return(umap)
}

#' @keywords internal
validateDesignFormula <- function(formula, contrast = NULL, override = FALSE, verbose = FALSE) {
  if (is.null(formula)) {
    stop("Design formula must be provided.")
  }
  if (inherits(formula, "formula")) {
    formula_str <- deparse(formula)
    formula_str <- paste(formula_str, collapse = "")
  } else if (is.character(formula)) {
    formula_str <- formula
  } else {
    stop("Design formula must be a string or formula object.")
  }
  if (!grepl("~", formula_str)) {
    stop("Design formula must contain a '~' to separate response and predictors.")
  }

  containsRandomEffects <- grepl("\\([^\\|]*\\|[^\\)]*\\)", formula_str)
  if (containsRandomEffects) {
    rand_eff_vars <- unlist(regmatches(formula_str, gregexpr("(?<=\\|)[^\\)]+", formula_str, perl = TRUE)))
    rand_eff_vars <- trimws(rand_eff_vars)
    warning(sprintf(
      "Random effects (terms with '|') are not supported in this workflow. The variable(s) '%s' will be treated as fixed effects.",
      paste(rand_eff_vars, collapse = ", ")
    ))
    formula_str <- gsub("\\([^\\|]*\\|[^\\)]*\\)", "", formula_str)
    rhs <- gsub("~", "", formula_str)
    rhs <- gsub("\\++", "+", rhs)
    rhs <- gsub("^\\s*\\+|\\+\\s*$", "", rhs)
    rhs <- trimws(rhs)
    rhs_terms <- unlist(strsplit(rhs, "\\+"))
    rhs_terms <- trimws(rhs_terms)
    rhs_terms <- unique(c(rhs_terms, rand_eff_vars))
    rhs_terms <- rhs_terms[rhs_terms != ""]
    formula_str <- paste("~", paste(rhs_terms, collapse = " + "))
  }
  if (verbose) {
    if (exists("containsRandomEffects") && containsRandomEffects) {
      message("Random effect terms provided in the model will be treated as fixed effect terms.")
    }
    message(sprintf("Final design formula: %s", formula_str))
  }

  parsedTerms <- terms(as.formula(formula_str))
  termLabels <- attr(parsedTerms, "term.labels")

  if (is.null(contrast)) {
    # Extract first term from formula
    if (length(termLabels) == 0) {
      stop("No terms found in the design formula to use as contrast.")
    }
    contrast_var <- termLabels[1]
    if (!contrast_var %in% names(self$sample_meta)) {
      stop(sprintf("Contrast variable '%s' not found in sample metadata.", contrast_var))
    }
    levels_contrast <- levels(factor(self$sample_meta[[contrast_var]]))
    if (length(levels_contrast) < 2) {
      stop(sprintf("Contrast variable '%s' must have at least two levels.", contrast_var))
    }
    # Compare last level vs first level
    contrast <- c(contrast_var, levels_contrast[1], levels_contrast[length(levels_contrast)])   
  }

  if (is.null(contrast) || length(contrast) != 3) {
    stop("Contrast must be a vector of length 3: c('variable', 'level1', 'level2').")
  }

  if (override) {
    self$formula <- parsedTerms
    self$contrast <- contrast
  }

  return(parsedTerms)
}

#' @keywords internal
constructModelMatrix <- function(
  sample_meta, formula,
  contrast = NULL,
  override = FALSE,
  keep.intercept = FALSE,
  verbose = FALSE
) {
  if (override && !is.null(formula) && !is.null(contrast)) {
    formula <- validateDesignFormula(formula = formula, contrast = contrast, override = override, verbose = verbose)
  }
  predictorVars <- all.vars(formula)
  sample_meta <- sample_meta[, predictorVars, drop = FALSE]

  valid_vars <- sapply(predictorVars, function(cov) {
        x <- sample_meta[[cov]]
        if (is.factor(x) || is.character(x)) {
            length(unique(x)) > 1
        } else {
            var_x <- var(x, na.rm = TRUE)
            !is.na(var_x) && var_x > 0
        }
    })
    
    if (!all(valid_vars)) {
        removed_vars <- predictorVars[!valid_vars]
        warning(sprintf(
            "Removing covariates with only one factor level or zero variance: %s",
            paste(removed_vars, collapse = ", ")
        ))
        predictorVars <- predictorVars[valid_vars]
        sample_meta <- sample_meta[, predictorVars, drop = FALSE]
        intercept_present <- attr(terms(formula), "intercept") == 1
        formula <- reformulate(predictorVars, intercept = intercept_present)
    }

  #sample_meta <- filterDEMetadata(sample_meta)
  #if (is.null(sample_meta)) {
  #  stop("Redundant covariate columns detected.")
  #}
  predictorVars <- colnames(sample_meta)
  n_samples <- nrow(sample_meta)
  n_covs <- length(predictorVars)
  if (n_samples < n_covs + 1) {
    stop(sprintf("Too few samples (%d) relative to number of covariates (%d).", n_samples, n_covs))
  }
  if (length(predictorVars) > 1) {
        to_remove <- c()
        for (i in seq_along(predictorVars)) {
            for (j in seq_along(predictorVars)) {
                if (i < j) {
                    col1 <- sample_meta[[predictorVars[i]]]
                    col2 <- sample_meta[[predictorVars[j]]]
                    if (is.factor(col1)) col1 <- as.character(col1)
                    if (is.factor(col2)) col2 <- as.character(col2)
                    if (all(col1 == col2, na.rm = TRUE)) {
                        warning(sprintf("Covariates '%s' and '%s' are identical. Removing '%s'.", predictorVars[i], predictorVars[j], predictorVars[j]))
                        to_remove <- c(to_remove, predictorVars[j])
                    } else if (is.numeric(col1) && is.numeric(col2)) {
                        cor_val <- suppressWarnings(cor(col1, col2, use = "pairwise.complete.obs"))
                        if (!is.na(cor_val) && abs(cor_val) > 0.95) {
                            warning(sprintf(
                                "Covariates '%s' and '%s' are highly correlated (cor = %.2f).",
                                predictorVars[i], predictorVars[j], cor_val
                            ))
                        }
                    }
                }
            }
        }
        if (length(to_remove) > 0) {
            to_remove <- unique(to_remove)
            sample_meta <- sample_meta[, !(colnames(sample_meta) %in% to_remove), drop = FALSE]
            predictorVars <- colnames(sample_meta)
            intercept_present <- attr(terms(formula), "intercept") == 1
            formula <- reformulate(predictorVars, intercept = intercept_present)
        }
    }
  for (cov in predictorVars) {
    if (!is.numeric(sample_meta[[cov]]) && !is.factor(sample_meta[[cov]])) {
      if (isTRUE(verbose)) message(sprintf("Converting covariate '%s' to factor.", cov))
      unique_vals <- sort(unique(na.omit(sample_meta[[cov]])))
      sample_meta[[cov]] <- factor(sample_meta[[cov]], levels = unique_vals)
    }
  }
  
  if (!is.null(contrast)) {
    contrast_var <- contrast[1]
    ref.level <- contrast[2]
    target.level <- contrast[3]
  } else {
    contrast_var <- self$contrast[1]
    ref.level <- self$contrast[2]
    target.level <- self$contrast[3]
  }
  if (!is.null(contrast_var) && contrast_var %in% names(sample_meta)) {
    unique_levels <- unique(sample_meta[[contrast_var]])
    if (!ref.level %in% unique_levels || !target.level %in% unique_levels) {
      stop("ref.level or target.level not found in levels of ", contrast_var)
    }
    sample_meta[[contrast_var]] <- factor(
      sample_meta[[contrast_var]],
      levels = c(ref.level, setdiff(unique_levels, ref.level))
    )
  }
  model_mat <- model.matrix(formula, data = sample_meta)
  if (!keep.intercept) {
    if ("(Intercept)" %in% colnames(model_mat)) {
      model_mat <- model_mat[, -1, drop = FALSE]
    }
  }
  # check rank
  qr_decomp <- qr(model_mat)
  if (qr_decomp$rank < ncol(model_mat)) {
    warning(sprintf(
      "Model matrix has linear dependencies: rank %d < number of columns %d. Possible confounding or redundant covariates.",
      qr_decomp$rank, ncol(model_mat)
    ))
  }
  dimnames(model_mat) <- list(rownames(model_mat), colnames(model_mat))
  return(model_mat)
}

#' @keywords internal
getSampleGroups <- function(sample.meta, contrast, sample.id) {
  if (!is.null(contrast)) {
    sample.groups <- setNames(as.character(sample.meta[[contrast[1]]]), sample.meta[,sample.id])
    sample.groups <- sample.groups[sample.groups %in% c(contrast[2], contrast[3])]
    sample.groups <- factor(sample.groups, levels = c(contrast[2], contrast[3]))
  }
  return(sample.groups)
}