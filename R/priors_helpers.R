#' Fill in sensible default priors for a given p_family,
#' then override with any user-supplied priors.
#'
#' @param priors A named list of user-specified priors, e.g.
#'   list(beta1 = "normal(0,5)", theta = "beta(2,2)")
#' @param p_family One of "linear", "poisson", "logistic", or "gamma"
#' @return A named list of priors in the same format as `priors`,
#'   containing one entry for every required parameter.
#' @keywords internal

fill_defaults <- function(priors = list(), p_family) {
  defaults <- switch(
    p_family,
    "linear" = list(
      mu1 = "normal(0,5)",
      beta1 = "normal(0,5)",
      sigma1 = "cauchy(0,2.5)",
      mu2 = "normal(0,5)",
      beta2  = "normal(0,5)",
      sigma2 = "cauchy(0,2.5)",
      theta = "beta(1,1)"
    ),
    "poisson" = list(
      beta1 = "normal(0,5)",
      beta2 = "normal(0,5)",
      theta = "beta(1,1)"
    ),
    "logistic" = list(
      beta1 = "normal(0, 2.5)",
      beta2 = "normal(0,5)",
      theta = "beta(1,1)"
    ),
    "gamma" = list(
      beta1 = "normal(0,5)",
      beta2 = "normal(0,5)",
      shape = "gamma(2,0.1)",
      theta = "beta(1,1)"
    ),
    stop("`p_family` must be one of 'linear','poisson','logistic','gamma'")
  )
  # merge user overrides into the defaults
  utils::modifyList(defaults, priors)
}


#' Lightly validate user-supplied arguments; further validation occurs during
#' iteration over list elements in process_Variables
#' 
#' @param priors A named list of Stan‐style prior strings, e.g.
#' list(mu1 = "normal(0,5))", sigma1 = "cauchy(0,2.5))
#' @param p_family One of "linear", "poisson", "logistic", or "gamma"
#' @return Invisibly TRUE if all checks pass or only raise warning; errors otherwise.
#' @keywords internal
are_valid_args <- function(priors, p_family) {
  # p_family must be one of the supported families
  supported_families <- c("linear","poisson","logistic","gamma")
  if (!(p_family %in% supported_families)) {
    stop("`p_family` must be one of: ", paste(supported_families, collapse = ", "))
  }
  
  # priors list must be a named list
  if (!is.list(priors) || is.null(names(priors))) {
    stop("`priors` must be a named list of prior strings")
  }
  
  #some blocks were removed because they dont make sense for the implementation
  
  # this block is temporary disabled---needs major changes
  # 5) string must look like dist(…)
  #    allowing both univariate and multivariate families
  #tested_dists <- c("normal", "cauchy", "beta", "gamma", "multi_normal")
  #pattern <- sprintf("^(%s)\\s*\\((.+)\\)$", paste(tested_dists, collapse = "|"))
  #bad_syntax <- names(priors)[!grepl(pattern, priors)]
  #if (length(bad_syntax) > 0) {
  #  stop("Invalid prior syntax for: ", paste(bad_syntax, collapse = ", "),
  #       ". Each must look like dist(args), where dist ∈ {",
  #       paste(tested_dists, collapse = ", "), "}.")
  #}
  
  invisible(TRUE)
}

process_variable <- function(value, key) {
  if (is.character(value)) {
    # dont include strings i.e. priors themselves ex. normal(x,y) 
    # in final definitions block
    return(list(declaration="", definition="")) 
  } else if (is.matrix(value)) { 
    if (nrow(value) == ncol(value) && isSymmetric(value)) {
      dim = nrow(value)
    } else {
      stop("Non-square or unsymmetric matrix found in list 'priors'")
    }
    dec = paste0("cov_matrix[", dim, "] ", key, ";") # declaration
    component_defs = ""
    for (i in 1:dim) {
      for (j in 1:dim) {
        component <- paste0(key, "[", i, ", ", j, "] = ", value[i, j], ";" )
        component_defs <- paste0(component_defs, component)
      }
    }
    return(list(declaration=dec, definition=component_defs))
  } else if (is.numeric(value)) { # non-matrix numeric
    if (length(value) == 1L) { # not a vector; scalar
      def = paste0("real ", key, " = ", value, ";")
      return(declaration="", definition=def)
    } else { # value is a vector
      len = length(value)
      # convert vector to stan list as a string:
      stan_vector = paste0("[", paste(value, collapse=", "), "]'")
      def = paste0("vector[", len, "] ", key, " = ", stan_vector, ";")
      return(list(declaration="", definition=def))
    }
  } else {
    stop("Element of invalid type or form found in list 'priors'")
  }
}