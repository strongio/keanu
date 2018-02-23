
# Keanu Model Factory -------------------------------------------------------------------------


#' Create your own modelling function, given a log-likelihood function.
#'
#' This function takes a log-likelihood function, and returns a modelling function, much like
#' \code{stats:lm} or \code{stats:glm}. This modelling function will find coefficients that maximize
#' the likelihood.
#'
#' @param loglik_fun A function whose first argument is the target variable, and whose subsequent
#'   arguments are parameters of the distribution which gives the log-likelihood. (You can indicate
#'   arguments that are *not* parameters with \code{non_parameter_args}.)
#' @param penalty_fun Optional. A function which, when given a dataframe with columns \code{c('parameter',
#'   'term', 'beta')}, yields a penalty term (a positive number to be added to the function to be
#'   minimized by \code{stats::optim}).
#' @param non_parameter_args A character-vector indicating arguments to \code{loglik_fun} that are
#'   not parameters. (This is in addition to the first argument, which is always assumed to be the
#'   target variable, not a parameter.)
#'
#' @return A modelling function. See \code{modeller_template}.
#' @export
modeller_from_loglik <- function(loglik_fun,
                                 penalty_fun = function(...) return(0),
                                 non_parameter_args = character()
) {
  .parameters <- setdiff(names(formals(loglik_fun))[-1], non_parameter_args)
  out_fun <- modeller_template
  environment(out_fun) <- list2env(x = list(loglik_fun = loglik_fun,
                                            penalty_fun = penalty_fun,
                                            .parameters = .parameters),
                                   parent = globalenv())
  return(out_fun)
}

#' A "template" for a modelling function
#'
#' This function will not work by itself. Instead, \code{modeller_from_loglik} provides a modified version of it.
#'
#' @param formulas A list of formulas. The LHS should contain the dependent variable wrapped around a parameter-- e.g., \code{mu(y) ~ x}.
#' @param data A data.frame.
#' @param standardize_x Should data be standardized before passing to optimizer? Can help with convergence. Affects interpretation of coefficients.
#' @param na.action Action when records in data are NA
#' @param loglik_fun_args Further arguments to pass to the loglik function (from \code{modeller_from_loglik}).
#' @param penalty_fun_args Further arguments to pass to the penalty function (from \code{modeller_from_loglik}).
#' @param optim_args A list of arguments to pass to \code{stats:optim}.
#'
#' @return An object of class 'keanu_model', with print and predict methods.
#' @export
modeller_template <- function(formulas, data, standardize_x = FALSE, na.action = na.warn,
                              loglik_fun_args = list(),
                              penalty_fun_args = list(),
                              optim_args = list(method = "BFGS", control = list(trace=as.integer(interactive()),maxit=250)),
                              contrasts = NULL
) {
  the_call <- match.call()
  stopifnot(is.list(formulas))
  is_forms <- purrr::map_lgl(formulas, rlang::is_formula)
  has_both_sides <- purrr::map_int(formulas, length) == 3

  if (!all(is_forms) | !all(has_both_sides)) stop(call. = FALSE, "`formulas` must be a list of formulas with both RHS and LHS.")
  form_info <- parse_formulas_list(formulas = formulas, .parameters=.parameters)

  mc <- prepare_model_components(forms = form_info$formulas, data = data,
                                 standardize_x = standardize_x, na.action = na.action,
                                 contrasts = contrasts)

  out <- optim_from_mc(mc, purrr::map(form_info$links, 'inv_link'), optim_args, loglik_fun, loglik_fun_args, penalty_fun, penalty_fun_args)

  out$call <- the_call
  out$model_components <- mc
  out$parameters <- .parameters
  par_tbl <- roll_par_tbl(out$optim$par)
  par_tbl_split <- split(x = select(par_tbl, -parameter), f = par_tbl$parameter)
  out$coefficients <- purrr::map(.x = par_tbl_split, .f = tibble::deframe)
  out$loglik_fun <- loglik_fun
  out$penalty_fun <- penalty_fun
  out$links <- form_info$links
  out$optim_args <- optim_args

  class(out) <- c("keanu_model", class(out))
  out

}

#' Given the output of `prepare_model_components`, find coefficient-estimates using `stats::optim`.
#'
#' @param mc The output of `prepare_model_components`
#' @param inv_links A list of inverse-link-functions, one corresponding to each `mc$model_mats`
#' @param optim_args A named list of arguments to be passed to `stats::optim`.
#' @param loglik_fun A log-liklihood function to optimize
#' @param loglik_fun_args Additional arguments (aside from target and parameters) to pass to the
#'   `loglik_fun`.
#' @param penalty_fun A penalty function that takes a table with parameter, term, beta, and returns
#'   a numeric penalty that will be factored into the optimization.
#' @param penalty_fun_args Additional arguments (aside from the table) to pass to the `penalty_fun`.
#'
#' @return A list with 'optim' (the results of the call to stats::optim), and 'cov' (the approximate
#'   covariance matrix obtaining from the hessian).
#' @export
optim_from_mc <- function(mc, inv_links, optim_args, loglik_fun, loglik_fun_args, penalty_fun, penalty_fun_args) {
  par_tbl_init <- initialize_par_tbl(mc$model_mats)
  par_init <- unroll_par_tbl(par_tbl_init)

  fun_to_minimize <- function(par, list_of_mm, Y, inv_link_funs) {
    par_tbl <- roll_par_tbl(par)
    df_dist_params <- get_distribution_params(par_tbl, list_of_mm, inv_link_funs = inv_link_funs)
    # loglik:
    ll_args <- c(list(quote(Y)), as.list(df_dist_params), loglik_fun_args)
    lik <- do.call(loglik_fun, ll_args)
    # penalty:
    pen_args <- c(as.list(par_tbl)[1], as.list(par_tbl)[-1], penalty_fun_args)
    pen <- do.call(penalty_fun, pen_args)
    # to minimize:
    sum(c(-lik , pen))
  }

  out <- list()
  optim_args$par <- par_init
  optim_args$fn <- fun_to_minimize
  if ( !is.element('hessian', names(optim_args)) ) optim_args$hessian <- TRUE
  optim_args$list_of_mm <- quote(mc$model_mats)
  optim_args$Y <- quote(mc$response_object)
  optim_args$inv_link_funs <- inv_links
  out$optim <- do.call(stats::optim, optim_args)

  out$cov <- matrix(nrow = length(par_init), ncol = length(par_init), dimnames = list(names(par_init),names(par_init)))
  if (!is.null(out$optim$hessian) &&
      all(!is.na(out$optim$hessian)) &&
      all(!is.nan(out$optim$hessian)) &&
      all(is.finite(out$optim$hessian)) &&
      all(eigen(out$optim$hessian)$values > 0)) {
    fixed_lgl <- logical(length(par)) # TO DO
    replace_idx <- matrix(!fixed_lgl,nrow = length(fixed_lgl), ncol = length(fixed_lgl)) &
      matrix(!fixed_lgl,nrow = length(fixed_lgl), ncol = length(fixed_lgl), byrow = TRUE)
    out$cov[replace_idx] <- solve(out$optim$hessian)

  } else {
    if (optim_args$hessian)
      warning(call. = FALSE,
              "Optimization has probably not converged - Hessian is not positive definite. ")
  }
  return(out)
}


#' Parse links: a helper function for `modeller_template`
#'
#' @param the_link A character naming a link or a named list of functions with elements
#' 'inv_link' and 'link'.
#'
#' @return A named list of functions with elements 'inv_link' and 'link'.
#' @export
parse_link <- function(the_link) {
  if (is.null(the_link)) the_link <- "identity"
  if (is.list(the_link)) {
    if (!all(c("link","inv_link") %in% names(the_link)))
      stop(call. = FALSE, "Link should either be a character or a list with 'link' and 'inv_link'.")
    if (!all(purrr::map_lgl(the_link, purrr::is_function)))
      stop(call. = FALSE, "Link is a list, but one or more elements is not a function.")
    return(the_link)
  } else {
    stopifnot(is.character(the_link))
  }
  out <- switch (the_link,
                 "log" = list(link = log, inv_link = exp),
                 "logit" = list(link = qlogis, inv_link = plogis),
                 "identity" = list(link = identity, inv_link = identity),
                 FALSE
  )
  if (identical(out,FALSE)) stop(call. = FALSE, glue::glue("Did not recognize link function {the_link}."))
  return(out)
}

#' Parse formulas: a helper function for `modeller_template`
#'
#' @param formulas A list of formulas.
#' @param .parameters A character-vector of the expected parameters.
#'
#' @return A list consisting of (1) a list of reformatted formulas ready to be
#'   passed to `prepare_model_components` and (2) a list of link/inverse-link pairs.
#' @export
parse_formulas_list <- function(formulas, .parameters) {
  eval_env <- environment(formulas[[1]])
  for (i in seq_along(formulas)) environment(formulas[[i]]) <- eval_env
  lh_sides <- purrr::map(formulas, function(f) f[[2]])
  if (any(purrr::map_lgl(lh_sides, is.symbol))) {
    stop("The left-hand side of your formula should be the DV wrapped with a parameter of the loglik function (e.g., mu(dv)).",
         call. = FALSE)
  }
  lh_params <- purrr::map_chr(purrr::map(lh_sides, function(s) s[[1]]), deparse)
  unexpected_params <- setdiff(lh_params, .parameters)
  if (length(unexpected_params)>0)
    stop(call. = FALSE, "The following were not expected on the left-hand side of the formula: ",
         paste(collapse = ", ", unexpected_params))
  dv <- unique(purrr::map(lh_sides, function(s) s[[2]]))
  if (length(dv)>1) stop(call. = FALSE, "The DV must be identical across all formulas.")
  dv <- dv[[1]]
  links <- purrr::map(lh_sides, function(f) parse_link(eval(f$link, envir = eval_env)))

  forms_out <- list()
  for (i in seq_along(formulas)) {
    this_form <- formulas[[i]]
    this_form[[2]] <- dv
    forms_out[[ lh_params[[i]] ]] <- this_form
  }
  names(links) <- names(forms_out)

  return( list(formulas = forms_out, links = links) )
}


# Keanu Model Methods -------------------------------------------------------------------------

#' A helper function for summarizing keanu models.
#'
#' Used for `print` and `plot_coefficients`.
#'
#' @param object A `keanu_model`.
#' @param se_multi Standard-errors will be mulitplied by this to obtain confidence-intervals.
#'
#' @return A data_frame with coefficient-estimates and lower, upper confidence intervals.
#' @export
.estimate_summary <- function(object, se_multi = 1.96) {
  df_est <- roll_par_tbl( object$optim$par )
  out <- left_join(x = rename(df_est, estimate=beta),
                   y = rename(roll_par_tbl(sqrt(diag(object$cov))), se=beta),
                   by = c('parameter','term')) %>%
    mutate(lower = estimate - se_multi*se,
           upper = estimate + se_multi*se,
           se=NULL)
  out_split <- split(out, out$parameter)
  inv_links <- map(object$links, 'inv_link')
  out_split <- map2(out_split, inv_links[names(out_split)],
                    function(df, inv_link) {
                      mutate_at(.tbl = df, .vars = c('estimate','lower','upper'), .funs = funs(ilink=inv_link))
                    })
  needs_ilink <- map_lgl(out_split, function(df) any(!near(df$estimate, df$estimate_ilink)))
  out <- map2_df(out_split, needs_ilink, function(df, keep_ilink) {
    if (keep_ilink) {
      df1 <- select(df, -ends_with('_ilink')) %>%
        mutate(parameter = paste0("link(", parameter, ")"))
      df2 <- df %>%
        select(-(estimate:upper)) %>%
        rename_all(.funs = funs(gsub(., pattern = "_ilink", replacement = "")))
      return( bind_rows(df1, df2) )
    } else {
      return( select(df, -ends_with('_ilink')) )
    }
  })
  return(out)
}

#' @export
print.keanu_model <- function(x, ...) {
  print(.estimate_summary(x), n=Inf)
}

#' @export
predict.keanu_model <- function(object, newdata, na.action = na.pass, ...) {
  mc <- prepare_model_components(forms = purrr::map(object$model_components$forms, remove_lhs),
                                 data = newdata,
                                 predvars = object$model_components$predvars,
                                 standardize_x = object$model_components$standardize_x,
                                 na.action = na.action,
                                 contrasts = object$model_components$contrasts,
                                 xlev = object$model_components$xlev,
                                 drop.unused.levels = FALSE)

  out <- get_distribution_params(par_tbl = roll_par_tbl( object$optim$par ),
                                 model_mats = mc$model_mats,
                                 inv_link_funs = map(object$links, 'inv_link'))
  as.matrix(out)
}

#' @export
update.keanu_model <- function(object, forms. = NULL, ...) {
  call <- getCall(object)
  eval_env <- environment(object$model_components$forms[[1]])
  extras <- match.call(expand.dots = FALSE)$...
  if (!is.null(forms.)) {
    call$formulas <- purrr::map2(eval(call$formulas, eval_env), forms., update)
    eval_env <- environment(forms.[[1]])
  }
  if (length(extras) > 0) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) {
      call <- c(as.list(call), extras[!existing])
      call <- as.call(call)
    }
  }
  eval(call, envir =  eval_env)
}

#' @export
augment.keanu_model <- function(x, data) {
  preds <- as.data.frame(predict(x, data))
  colnames(preds) <- paste(sep = "_", ".fitted", colnames(preds))
  data[colnames(preds)] <- preds
  return(data)
}

#' @export
logLik.keanu_model <- function(object, ...) {
  orig_env <- environment(object$model_components$forms[[1]])
  loglik_fun <- environment(eval(object$call[[1]], envir = orig_env))$loglik_fun
  data <- eval(pryr::standardise_call(object$call, env = orig_env)$data, envir = orig_env)
  na_dropped_idx <- object$model_components$na_dropped_idx
  if (is.null(na_dropped_idx)) na_keep_idx <- seq_len(nrow(data))
  else na_keep_idx <- -na_dropped_idx
  data <- data[na_keep_idx,,drop=FALSE]
  args <- as.data.frame(predict(object, newdata=data, envir = orig_env))
  args <- c(list(object$model_components$response_object),args)
  lls <- sum(do.call(loglik_fun, args))
  attr(lls, "df") <- length(object$optim$par)
  return(lls)
}

#' @export
plot_coefficients.keanu_model <- function(x, terms = NULL, se_multi = 1.96, link_only=TRUE) {
  df_est <- .estimate_summary(x, se_multi=se_multi)

  # parse terms:
  if (is.null(terms)) terms <- dplyr::vars(-`(Intercept)`)
  stopifnot(is.character(terms) || inherits(terms, "col_list") || inherits(terms,"quosures") )
  if (is.character(terms)) terms <- paste0("`", terms, "`")
  terms <- dplyr::select_vars_(vars = unique(df_est$term), args = terms)

  # prep table:
  df_est <- filter(.data = df_est, term %in% terms)
  if (link_only) {
    df_est <- df_est %>%
      group_by(parameter2 = inside_parens(df_est$parameter)) %>%
      filter(if_else(rep(any(parameter!=parameter2),n()), parameter!=parameter2, TRUE)) %>%
      ungroup() %>%
      select(-parameter2)
  }
  df_est$parameter <- forcats::fct_inorder(df_est$parameter)

  # prep plot:
  ggplot(df_est, aes(x=term, y=estimate)) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    facet_wrap(~parameter, scales = "free", labeller='label_both') +
    geom_hline(yintercept = 0, linetype='dashed') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


#' @export
univariate_predictions.keanu_model <- function(object) {
  df_mf <- object$model_components$model_frame_merged
  averages <- purrr::map(df_mf, function(x) {
    if (is.numeric(x))
      return(mean(x,na.rm = TRUE))
    else
      return(names(sort(table(x), decreasing = TRUE))[[1]])
  })
  df_avg <- select(df_mf)
  for (var in colnames(df_mf)) {
    df_avg[[var]] <- averages[[var]]
  }
  out <- list()
  for (var in colnames(df_mf)) {
    df_this_var <- df_avg
    df_this_var[[var]] <- df_mf[[var]]
    out[[var]] <- df_this_var[,var,drop=FALSE]
    preds <- as.data.frame(predict(object, df_this_var))
    colnames(preds) <- paste0(".fitted_",colnames(preds))
    out[[var]][colnames(preds)] <- preds
  }
  return(out)
}





