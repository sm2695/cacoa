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
buildModelMatrix <- function(
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

# Build a p×K contrast matrix from:
#   - X: model matrix used for fitting (preferably with intercept kept)
#   - formula: the model formula used to build X
#   - contrasts: list of c(var, ref, alt) triplets, e.g. list(c("diagnosis","Control","ASD"))
#   - expand:
#       "marginal" -> only main contrasts for each requested var
#       "simple"   -> main + simple effects at each non-baseline level of partner factors found in interactions
#       "both"     -> "simple" plus, when a partner is 2-level, include the difference-of-simple-effects (the pure interaction)
#
# Notes:
# - Works with treatment coding (the default from model.matrix). If you dropped the intercept, it still works.
# - Supports factor×factor and continuous×factor interactions. For >2-level partners, it generates simple effects per level.
buildContrastMatrix <- function(X, formula, contrasts, expand = c("marginal","simple","both")) {
  stopifnot(is.matrix(X))
  expand <- match.arg(expand)
  p  <- ncol(X); cn <- colnames(X)

  # ---- helpers -------------------------------------------------------
  # Main-effect dummy columns for a factor var (exclude interactions)
  main_level_to_col <- function(var) {
    idx <- grep(paste0("^", var), cn)
    idx <- idx[!grepl(":", cn[idx], fixed = TRUE)]
    lev <- sub(paste0("^", var), "", cn[idx])
    keep <- nzchar(lev)
    setNames(idx[keep], lev[keep])
  }
  # Check if a variable is continuous in X (has a single main column named exactly var)
  cont_col <- function(var) {
    j <- which(cn == var)
    if (length(j) == 1L) j else integer(0)
  }
  # Build the regex for an interaction column (order-agnostic)
  tok  <- function(var, lev) if (missing(lev) || is.null(lev)) var else paste0(var, lev)
  find_inter_col <- function(t1, t2) {
    which(cn == paste0(t1, ":", t2) | cn == paste0(t2, ":", t1))
  }
  add <- function(v, j, w) { if (length(j)==1 && j>0) v[j] <- v[j] + w; v }

  # Parse formula to discover interactions containing each requested var
  tr <- terms(formula)
  term_labels <- attr(tr, "term.labels")
  inter_pairs <- strsplit(term_labels[grepl(":", term_labels, fixed = TRUE)], ":", fixed = TRUE)

  partners_for <- function(var) {
    if (length(inter_pairs) == 0) return(character())
    unique(unlist(lapply(inter_pairs, function(ab) {
      if (var %in% ab) setdiff(ab, var) else character()
    })))
  }

  # Normalize contrasts input into list of lists
  specs <- lapply(contrasts, function(x) {
    stopifnot(is.character(x), length(x) == 3)
    list(var = x[1], ref = x[2], alt = x[3])
  })

  C <- NULL; cnames <- character()

  for (sp in specs) {
    A <- sp$var; ref <- sp$ref; alt <- sp$alt
    v_main <- numeric(p); names(v_main) <- cn

    # ----- main contrast for A (alt vs ref) -----
    mapA <- main_level_to_col(A)
    jA_ref <- unname(mapA[ref]); if (length(jA_ref)==0) jA_ref <- NA_integer_
    jA_alt <- unname(mapA[alt]); if (length(jA_alt)==0) jA_alt <- NA_integer_
    jA_cont <- cont_col(A)

    if (length(jA_cont) == 1L) {
      # A is continuous -> "main contrast" is just slope of A
      v_main <- add(v_main, jA_cont, +1)
      main_name <- paste0("Slope(", A, ")")
    } else {
      # A is factor
      if (!is.na(jA_alt) && !is.na(jA_ref)) {
        v_main <- add(v_main, jA_alt, +1); v_main <- add(v_main, jA_ref, -1)
      } else if (!is.na(jA_alt) && is.na(jA_ref)) {
        v_main <- add(v_main, jA_alt, +1)          # ref is baseline (intercept coding)
      } else if (is.na(jA_alt) && !is.na(jA_ref)) {
        v_main <- add(v_main, jA_ref, -1)          # alt is baseline
      } else {
        stop("For '", A, "': neither '", ref, "' nor '", alt, "' has a main column in X.")
      }
      main_name <- paste0(A, alt, "_vs_", ref)
    }

    # always include the marginal main contrast
    C <- cbind(C, v_main); cnames <- c(cnames, main_name)

    # ---------- expansions through interactions ----------
    if (expand != "marginal") {
      partners <- partners_for(A)
      for (B in partners) {
        # continuous partner?
        jB_cont <- cont_col(B)
        if (length(jB_cont) == 1L) {
          # A × continuous partner is uncommon for "simple effect"; skip by default
          next
        }

        # factor partner: simple effects at each non-baseline level of B
        mapB <- main_level_to_col(B)
        if (length(mapB) == 0) next  # B not factor-coded in X

        # Determine B's non-baseline levels present as columns
        levB <- names(mapB)  # these are the dummy-coded levels (exclude baseline)
        for (b in levB) {
          v_simple <- v_main  # start from main effect of A
          # add interaction column for A_alt : B_b (or slope(A) : B_b if A is continuous)
          if (length(jA_cont) == 1L) {
            jInt <- find_inter_col(tok(A), tok(B, b))
            if (length(jInt) == 1L) v_simple <- add(v_simple, jInt, +1)
            simple_name <- paste0(main_name, " | ", B, "=", sub(paste0("^", B), "", b))
          } else {
            jInt <- find_inter_col(tok(A, alt), tok(B, b))
            if (length(jInt) == 1L) v_simple <- add(v_simple, jInt, +1)
            simple_name <- paste0("Simple(", A, " ", alt, "_vs_", ref, " | ",
                                  B, "=", sub(paste0("^", B), "", b), ")")
          }
          C <- cbind(C, v_simple); cnames <- c(cnames, simple_name)
        }

        # Differences of simple effects (pure interaction) when B has exactly 1 dummy (i.e., 2 levels)
        if (expand == "both" && length(levB) == 1L) {
          b <- levB[1]
          v_diff <- numeric(p); names(v_diff) <- cn
          if (length(jA_cont) == 1L) {
            # slope(A|B=b) - slope(A|B=baseline) = beta_{A:B_b}
            jInt <- find_inter_col(tok(A), tok(B, b))
            if (length(jInt) == 1L) v_diff <- add(v_diff, jInt, +1)
            diff_name <- paste0("DiffSlope(", A, " | ", B, "=", sub(paste0("^", B), "", b), " vs baseline)")
          } else {
            # [A_alt vs ref at B=b] - [A_alt vs ref at baseline] = beta_{A_alt:B_b}
            jInt <- find_inter_col(tok(A, alt), tok(B, b))
            if (length(jInt) == 1L) v_diff <- add(v_diff, jInt, +1)
            diff_name <- paste0("Interaction(", A, " ", alt, " × ", B, "=", sub(paste0("^", B), "", b), ")")
          }
          if (any(v_diff != 0)) { C <- cbind(C, v_diff); cnames <- c(cnames, diff_name) }
        }
      }
    }
  }

  if (is.null(C)) C <- matrix(numeric(0), nrow = p, ncol = 0, dimnames = list(colnames(X), character()))
  rownames(C) <- colnames(X); colnames(C) <- cnames
  C
}