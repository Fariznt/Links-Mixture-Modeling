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


#' Helper function that lightly validates a string as (intended to be) a stan
#' function. Specifically, returns true if str is of the form 
#' ...(...)...{...return...}... where '...' is any arbitrary string AND the
#' string does not contain R function assignment syntax '<-'
is_valid_stan_func <- function(str) {
  txt <- gsub("[\r\n]+", " ", str)
  grepl("\\([^)]*\\).*\\{.*return.*\\}", txt, perl = TRUE) &&
    !grepl("<-", txt, fixed = TRUE)
}

#' Wraps a string element in priors list to indicate that the string is intended 
#' as a stan function.This string will be directly injected into a stan function 
#' block with little validation and is intended for niche or advanced use cases. 
#' Custom functions allow the definition of an arbitrary prior. 
#' @param stan_function_str A string defining a function in stan syntax
#' @return Inputted string wrapped to indicate the string as a stan function 
#' injection
#' @examples
#' lin_fit <- LinksMixtureModeling::fit_model(
#'   ... # other arguments here
#'   priors = list(
#'   my_function = stan_fun("
#'     real my_function(...) { 
#'       ... // function body
#'       return ... ; // return a prior distribution
#'     }
#'   "), # this function is injected into stan's function {} block
#'   theta1 = "my_function(...)" # reference function for other definitions
#'   )
#'   ... # other arguments here
#' )
#' 
stan_func <- function(stan_function_str) {
  if (!is.character(stan_function_str)) {
    stop("`stan_function_str` passed to stan_func() must be a single string.", call. = FALSE)
  } else if (!is_valid_stan_func(stan_function_str)) {
    stop("stan_func() called on invalid stan function")
  }
  structure(stan_function_str, class = c("stan_function_string", "character"))
}

#' Checks if a prior list element was passed in as a stan function; used for
#' validation and prior list processing
#' @param stan_function 
#' @return true if argument was constructed as stan_func("..."), false otherwise
is_func <- function(stan_function) {
  inherits(stan_function, "stan_function_string")
}

#' Initially validate user-supplied arguments; further validation occurs during
#' iteration over list elements in process_Variables. Primarily verifies
#' format for prior string follows <dist>(<comma separated args>) where
#' each arg is either a number or defined in the priors list.
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
  
  for (elt in priors) {
    if (is.character(elt) && !is_func(elt)) { # only looking at prior strings
      trimmed_elt <- trimws(elt) # remove whitespace
      
      # verify basic structure <dist>(<args>)
      pattern <- "^([^(\\s]+)\\s*\\(([^)]*)\\)\\s*$"
      match <- regexec(pattern, trimmed_elt, perl = TRUE)
      parts <- regmatches(trimmed_elt, match)[[1]]
      dist = parts[2]
      hyperparam = parts[3]
      
      if (length(parts) == 0 || hyperparam == "" || dist == "") {
        stop("Prior defining string must be of form '<dist>(<comma-separated args>)'")
      }
      
      # send warning about untested distributions
      if (!(dist %in% c("normal", "cauchy", "beta", "gamma", "multi_normal"))) {
        warning("An untested distribution was defined as prior. Tested ones", 
                " include but are not limited to normal, cauchy, beta, gamma, and",
                " multi_normal.")
      }
      
      # validate args
      if (nzchar(trimws(hyperparam))) { # if args has non-whitespace
        # get individual args
        raw_hyperparam_list <- strsplit(hyperparam, ",", fixed = TRUE)[[1]] 
        hyperparam_list <- trimws(raw_hyperparam_list) # trim white space
        
        if (any(hyperparam_list == "") || grepl(",\\s*$", hyperparam)) {
          stop("Missing values in list of hyperparameters of prior string")
        } 
        
        for (hyperparam in hyperparam_list) {
          # if hyperparameter is not a number and not a key defined in the priors
          if (is.na(suppressWarnings(as.numeric(hyperparam))) 
              && !(hyperparam %in% names(priors))) {
            stop("A variable used in a prior string is not defined.",
                 " If a string defining a prior in priors list is normal(x,y),",
                 " x and y must be their own keys in list 'priors' to be used.")
          }
        }
      }
    }
  }
  invisible(TRUE)
}

#' Processes a given key-value pair from the priors list into syntactically
#' valid stan to be inserted during stan generation. 
#' @param value Value ex. "normal(0,5)"
#' @param key key ex. beta1
#' @return a list that organizes stan code by definition, declaration, and/or 
#' stan function so it is apparent where the stan code should be inserted
process_variable <- function(value, key) {
  if (is.character(value)) { # is string
    if (is_func(value)) { # the string is a function i.e. stan_func("...stan function...")
      warning("Detected stan function definition. Note that stan function 
                injection does not fully validate stan syntax.")
      return(list(declaration="", definition="", stan_func=as.character(value)))
    } else { # the string is the prior definition itself ex. normal(x,y) 
      return(list(declaration="", definition="", stan_func=""))  
    }
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
    return(list(declaration=dec, definition=component_defs, stan_func=""))
  } else if (is.numeric(value)) { # non-matrix numeric
    if (length(value) == 1L) { # not a vector; scalar
      def = paste0("real ", key, " = ", value, ";")
      return(declaration="", definition=def)
    } else { # value is a vector
      len = length(value)
      # convert vector to stan list as a string:
      stan_vector = paste0("[", paste(value, collapse=", "), "]'")
      def = paste0("vector[", len, "] ", key, " = ", stan_vector, ";")
      return(list(declaration="", definition=def, stan_func=""))
    }
  } else {
    stop("Element of invalid type or form found in list 'priors'")
  }
}