#' keanu package
#'
#' @import dplyr
#' @import tidyr
#'
#' @importFrom tibble data_frame deframe
#' @importFrom purrr map map2 pmap map_dbl map_chr map_lgl
#'             map2_chr map2_dbl map2_lgl transpose flatten_chr flatten_dbl flatten_lgl
#'             walk walk2 map_df map2_df
#' @importFrom stats terms integrate delete.response as.formula model.response model.frame
#'             na.pass na.exclude na.fail model.matrix model.weights
#'             .getXlevels .checkMFClasses reformulate model.extract formula getCall update.formula
#' @importFrom pryr dots
#'
#' @docType package
#' @name keanu
NULL
#> NULL

#' Remove the LHS of a formula
#'
#' @param x A formula
#'
#' @return A formula, without its left-hand-side
#' @export
remove_lhs <- function(x) {
  stopifnot( inherits(x, "formula") )
  if (length(x)==3) return(x[-2])
  else return(x)
}

#' Get mapping from column names to model-terms and vice versa
#'
#' @param formula The formula to be passed to model.frame
#' @param data The data.frame
#' @param ... Arguments to be passed to \code{stats::model.matrix}
#'
#' @return A nested list, where the first level is 'from', the second level is 'to'.
#' @export
get_terms_mapping <- function(formula, data, ...) {

  flip_mapping <- function(mapping_list)  {
    df_mapper <- tibble::data_frame(map_from = names(mapping_list), map_to = mapping_list)
    df_mapper <- tidyr::unnest(df_mapper)
    df_mapper <- group_by(.data = df_mapper, .data$map_to)
    df_mapper <- ungroup(dplyr::do(.data = df_mapper, model_mat_cols = .$map_from))
    return( tibble::deframe(df_mapper) )
  }

  the_dots <- pryr::dots(...)
  extra_args_lgl <- c('model_frame','terms_obj') %in% names(the_dots)
  if (all(extra_args_lgl)) {
    model_frame <- eval(the_dots$model_frame, envir = parent.frame())
    the_dots$model_frame <- NULL
    terms_obj <- eval(the_dots$terms_obj, envir = parent.frame())
    the_dots$terms_obj <- NULL
  } else if (!any(extra_args_lgl)) {
    formula <- remove_lhs(formula)
    model_frame <- model.frame(formula = formula, data, na.action=na.pass)
    terms_obj <- terms(model_frame, data = model_frame)
  } else {
    stop(call. = FALSE, "Both 'model_frame' and 'terms_obj' or neither should be supplied.")
  }
  model_matrix <- do.call(model.matrix,
                          c( list(quote(terms_obj), data = quote(model_frame) ),
                             the_dots) )

  if (ncol(model_frame)==0) {
    # check for intercept?
    from_mm_to_od <- from_od_to_mm <- from_mf_to_od <- list()
    stop(call. = FALSE, "Please report this error to the package maintainer.")
  } else {
    # for each column in the model.matrix, get corresponding model.frame col(s):
    cols_in_mm <- colnames(model_matrix)
    fact_mat <- attr(terms(model_frame),'factors')
    from_mm_to_mf <- purrr::map(attr(model_matrix,'assign'),
                         function(assign_idx) row.names(fact_mat)[1==fact_mat[,assign_idx,drop=TRUE]])
    names(from_mm_to_mf) <- cols_in_mm

    # for each term.label, get corresponding model.frame col(s):
    from_tl_to_mf <-
      purrr::map(.x = seq_len(ncol(fact_mat)), .f = ~row.names(fact_mat)[1==fact_mat[,.x,drop=TRUE]])
    names(from_tl_to_mf) <- colnames(fact_mat)

    # for each column in the model.matrix, get corresponding original-data col(s):
    is_null_lgl <- purrr::map_lgl(from_mm_to_mf, is.null)
    from_mm_to_od <- vector(mode = 'list', length = length(from_mm_to_mf))
    names(from_mm_to_od) <- names(from_mm_to_mf)
    if (any(!is_null_lgl)) {
      from_mm_to_od[!is_null_lgl] <- purrr::map(
        .x = from_mm_to_mf[!is_null_lgl],
        .f = function(vec_of_mm_cols) purrr::flatten_chr(purrr::map(vec_of_mm_cols, ~all.vars(parse(text = .x)[[1]]))))
      # we parse the syntax for each term, getting any corresponding to a variable. then we have to check
      # if that variable is a column in the original-data (e.g., in case a term is
      # `poly(column, degree=my_variable_specifying_degree_that_isnt_a_column)`
      from_mm_to_od[!is_null_lgl] <- purrr::map(from_mm_to_od[!is_null_lgl], ~.x[.x%in%colnames(data)])
    }

    # for each column in the model.frame, get corresponding original-data col(s):
    from_mf_to_od <- purrr::map(colnames(model_frame), ~all.vars(parse(text = .x)[[1]])) %>%
      purrr::map(~.x[.x%in%colnames(data)])
    names(from_mf_to_od) <- colnames(model_frame)
  }

  #
  list(
    original_data = list(
      term_labels = purrr::map(flip_mapping(from_mf_to_od), ~purrr::flatten_chr(flip_mapping(from_tl_to_mf)[.x])),
      model_frame = flip_mapping(from_mf_to_od),
      model_matrix = flip_mapping(from_mm_to_od)
    ),
    term_labels = list(
      original_data = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(from_mf_to_od[.x])),
      model_frame = from_tl_to_mf,
      model_matrix = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(flip_mapping(from_mm_to_mf)[.x]))
    ),
    model_frame = list(
      original_data = from_mf_to_od,
      term_labels = flip_mapping(from_tl_to_mf),
      model_matrix = flip_mapping(from_mm_to_mf)
    ),
    model_mat = list(
      original_data = from_mm_to_od,
      term_labels = purrr::map(from_mm_to_mf, ~purrr::flatten_chr(flip_mapping(from_tl_to_mf)[.x])),
      model_frame = from_mm_to_mf
    )
  )

}

#' Merge a list of formulae
#'
#' @param forms List of formulae/formulas
#' @param data Data.frame
#'
#' @return Single formula, with \code{environment(forms[[1]])}
merge_formulae <- function(forms, data) {
  list_of_term_labels <- purrr::map(purrr::map(forms, terms, data=data), attr, 'term.labels')
  term_labels_unique <- unique(purrr::flatten_chr(list_of_term_labels))
  if (length(term_labels_unique)>0) form_out <- stats::reformulate(term_labels_unique)
  else form_out <- ~1
  lazyeval::f_lhs(form_out) <- lazyeval::f_lhs(forms[[1]])
  environment(form_out) <- environment(forms[[1]])
  form_out
}

#' Handle missing values with a warning
#'
#' Modified version from \code{modelr}, where \code{warning(immediate.=TRUE)}.
#'
#' @param object
#'
#' @export
na.warn <- function (object) {
  missing <- sum(!stats::complete.cases(object))
  if (missing > 0) {
    warning("Dropping ", missing, " rows with missing values",
            call. = FALSE, immediate. = TRUE)
  }
  stats::na.exclude(object)
}

#' Prepare model-components
#'
#' @param forms A list of formulas.
#' @param data A data.frame
#' @param predvars The 'predvars' attribute of the terms object, from a previous call to this
#'   function.
#' @param standardize_x Either a logical specifying whether to center/scale numeric predictors in
#'   the model.frame, or a list with names 'center' and 'scale', each of which in turn are named
#'   numeric vectors specifying centering/scaling. Defaults
#' @param na.action A function specifying what to do with NAs. Defaults to \code{na.warn}.
#' @param contrasts Passed to \code{model.matrix}.
#' @param xlev Passed to \code{model.frame}.
#' @param drop.unused.levels Passed to \code{model.frame}.
#'
#' @return A list of components useful for building models.
#' @export
prepare_model_components <- function(forms, data,
                                     predvars = NULL,
                                     standardize_x = FALSE,
                                     na.action = na.warn,
                                     contrasts = NULL,
                                     xlev = NULL,
                                     drop.unused.levels = FALSE) {

  formula_merged <- merge_formulae(forms, data)
  terms_merged <- terms(formula_merged)
  if (!is.null(predvars)) {
    if (is.character(predvars)) predvars <- parse(text = predvars)[[1]]
    if (length(formula_merged)==2)
      predvars <- predvars[-2] # if no response in formula, remove response from predvars
    attr(terms_merged,'predvars') <- predvars
  }
  model_frame_merged <-
    model.frame(terms_merged, data = data, na.action = na.action,
                drop.unused.levels = drop.unused.levels, xlev = xlev)
  xlevels <- .getXlevels(attr(model_frame_merged, "terms"), model_frame_merged)
  if (length(xlevels)==0) xlevels <- NULL

  # separate out response object: --
  response_idx <- attr(terms_merged,'response')
  terms_full <- terms(model_frame_merged)
  if (response_idx!=0) {
    response_object <- model_frame_merged[[response_idx]]
    model_frame_merged <- model_frame_merged[,-response_idx,drop=FALSE]
    attr(model_frame_merged,'terms') <- delete.response(terms_full)
  } else {
    response_object <- NULL
  }
  if (is.null(predvars)) predvars <- attr(terms_full,'predvars')

  ## standardize model-frame: --
  mf_is_numeric_lgl <- purrr::map_lgl(model_frame_merged, is.numeric) & !purrr::map_lgl(model_frame_merged, is.matrix)
  standardize_x_arg <- standardize_x
  standardize_x <- list(
    center= purrr::map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], mean, na.rm=TRUE),
    scale= purrr::map_dbl(model_frame_merged[,mf_is_numeric_lgl,drop=FALSE], sd, na.rm=TRUE)
  )
  if (!is.logical(standardize_x_arg)) {
    stopifnot(is.list(standardize_x_arg))
    stopifnot(c('center','scale') %in% names(standardize_x_arg))
    for (nm in c('center','scale')) {
      # check for factor names:
      if ( any( names(standardize_x_arg[[nm]]) %in% names(which(!mf_is_numeric_lgl)) ) )
        stop(call. = FALSE, "The following variables cannot be standardized because they aren't numeric:",
             paste0(deparse(names(which(!mf_is_numeric_lgl & does_something_lgl))),collapse="\n"),
             "\n(Hint: for factors, use contrasts for centering instead.)")

      # add to defaults:
      standardize_x[[nm]][names(standardize_x_arg[[nm]])] <- standardize_x_arg[[nm]]
    }
  } else {
    if (!standardize_x_arg) {
      standardize_x <- FALSE
    }
  }
  model_frame_std <- model_frame_merged
  if (!identical(standardize_x,FALSE)) {
    model_frame_std[names(standardize_x$center)] <- purrr::map(.x = names(standardize_x$center),
                                                               .f= ~model_frame_std[[.x]]-standardize_x$center[[.x]])
    model_frame_std[names(standardize_x$scale)] <- purrr::map(.x = names(standardize_x$scale),
                                                              .f= ~model_frame_std[[.x]]/standardize_x$scale[[.x]])
  }

  ## smarter contrasts: --
  if (length(contrasts)==0) contrasts <- NULL
  if (is.null(contrasts)) {
    matrix_contrasts_lgl <- purrr::map_lgl(purrr::map(model_frame_std, attr,'contrasts'), is.matrix)
    if (any(matrix_contrasts_lgl) & !is.null(xlevels))  # any variables with matrix-contrasts doesn't need an xlev (it'd override it)
      xlevels <- xlevels[intersect(names(xlevels),names(which(!matrix_contrasts_lgl)))]
    contrasts <- purrr::compact(purrr::map(model_frame_std, attr,'contrasts'))
  }

  ## get model-mat(s): --
  terms_obj <- terms(model_frame_std)
  model_matrix_merged_std <-
    model.matrix(terms_obj, data = model_frame_std, contrasts.arg = contrasts)
  term_mapping <- get_terms_mapping(data = data,
                                    contrasts.arg = contrasts,
                                    model_frame = model_frame_std,
                                    terms_obj = terms_obj)
  list_of_model_mats <- purrr::map(forms, function(this_form) {
    this_terms <- terms(this_form)
    needed_cols <- unique(purrr::flatten_chr(term_mapping$term_labels$model_matrix[attr(this_terms,'term.labels')]))
    if (attr(this_terms,'intercept')==1) needed_cols <- c("(Intercept)", needed_cols)
    model_matrix_merged_std[,needed_cols,drop=FALSE]
  })

  # return--
  out <- list(list_of_model_mats= list_of_model_mats,
              response_object = response_object,
              na_dropped_idx= attr(model_frame_merged,"na.action"),
              standardize_x= standardize_x,
              predvars = predvars,
              xlevels = xlevels,
              contrasts = contrasts,
              forms = forms)
  class(out) <- c("model_components",class(out))
  out

}


#' Take a data.frame with parameter, term, beta, unroll into a flat vector
#'
#' @param df_pars A data.frame with parameter, term, beta
#'
#' @return A named numeric vector, named with parameter__sep__term
#' @export
unroll_df_pars <- function(df_pars) {
  with(df_pars, structure(beta, names=paste0(parameter, "__sep__", term)))
}

#' Roll up a flat vector into a data.frame with parameter, term, beta
#'
#' @param unrolled_par A named numeric vector, named with parameter__sep__term
#'
#' @return A data.frame with parameter, term, beta
#' @export
roll_df_pars <- function(unrolled_par) {
  tidyr::separate(col = "name", into = c("parameter", "term"), sep = "__sep__", data = tibble::enframe(unrolled_par, value = 'beta'))
}

#' Given coefficients (in a data.frame) and a list of model-mats, get predicted parameters for the
#' distribution
#'
#' @param df_pars A data.frame with parameter, term, beta
#' @param list_of_model_mats A list of model-matrices
#' @param inv_link_funs A names list of inverse-link functions (names must match parameters)
#'
#' @return A data.frame where each row corresponds to the rows of the model-mats, and each column is
#'   a distribution-parameter.
#' @export
get_distribution_params <- function(df_pars, list_of_model_mats, inv_link_funs = NULL) {
  list_of_pars <- purrr::map(.x = split(df_pars, df_pars$parameter),
                             .f = function(df) tibble::deframe(df[,c('term','beta'),drop=FALSE]) )
  list_of_model_mats <- list_of_model_mats[names(list_of_pars)]
  if (is.null(inv_link_funs))
    inv_link_funs <- purrr::map(list_of_pars, ~identity)
  else
    inv_link_funs <- inv_link_funs[names(list_of_pars)]

  list_of_pars_per_obs <- purrr::map(
    .x = names(list_of_pars),
    .f = function(param_name) {
      inv_link_funs[[param_name]]( as.matrix(list_of_model_mats[[param_name]]) %*% matrix( list_of_pars[[param_name]] ) )
    })
  names(list_of_pars_per_obs) <- names(list_of_pars)
  tibble::as_data_frame(purrr::map(list_of_pars_per_obs, as.numeric))
}

#' Initialize a dataframe with parameter, term, beta
#'
#' @param list_of_model_mats A list of model-matrices
#'
#' @return A data.frame with parameter, term, beta
#' @export
initialize_df_pars <- function(list_of_model_mats) {
  nms <- purrr::map(list_of_model_mats, colnames)
  df_nested <- tibble::enframe(name = 'parameter', value = 'term', nms)
  df_unnested <- tidyr::unnest(df_nested)
  df_unnested$beta <- 0
  df_unnested
}

#' Given the results of `prepare_model_components`, find coefficients which maximize the
#' log-likelihood.
#'
#' @param model_components The results of `prepare_model_components`
#' @param log_likelihood_fun A function which takes two arguments: (1) a data.frame where the
#'   columns are parameters of the distribution being predicted, and rows correspond to
#'   observations, and (2) the target variable being predicted.
#' @param inits Not yet supported.
#' @param inv_link_funs A named list of inverse-link functions: i.e., functions that transform
#'   real-valued predictions into (potentially constrained) parameter-estimates. The list must be
#'   named for each parameter of the distribution. Defaults to the identity function for all.
#' @param optim_args A named list that will be passed as arguments to \code{stats::optim}.
#'
#' @return An object of class 'keanu_model', with a predict method.
#' @export
optim_from_mc <- function(
  model_components,
  log_likelihood_fun,
  inits = NULL,
  inv_link_funs = NULL,
  optim_args = list(method = "BFGS", control = list(trace=as.integer(interactive()),maxit=250)) ) {

  df_pars <- initialize_df_pars(model_components$list_of_model_mats)
  if (!is.null(inits))
    stop("Please report this error to the package maintainer.")
  par <- unroll_df_pars(df_pars)

  neg_loglik_fun <- function(par, list_of_mm, Y) {
    df_pars <- roll_df_pars(par)
    df_dist_params <- get_distribution_params(df_pars, list_of_mm, inv_link_funs = inv_link_funs)
    - sum( log_likelihood_fun(df_dist_params, Y) )
  }

  out <- list()
  optim_args$par <- par
  optim_args$fn <- neg_loglik_fun
  optim_args$hessian <- TRUE
  optim_args$list_of_mm <- quote(model_components$list_of_model_mats)
  optim_args$Y <- quote(model_components$response_object)
  out$optim <- do.call(optim,optim_args)

  ## collect results, add prior
  out$cov <- matrix(nrow = length(par), ncol = length(par),
                    dimnames = list(names(par),names(par)))
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
    warning(call. = FALSE,
            "Optimisation has probably not converged - Hessian is not positive definite. ")
  }
  df_est <- roll_df_pars( out$optim$par )
  out$res_trans <-
    left_join(x = rename(.data = df_est , estimate = beta),
              y = rename(.data = roll_df_pars( sqrt(diag(out$cov)) ), se = beta) ,
              by = c('parameter','term'))
  out$res_trans <- mutate(.data = out$res_trans,
                          upper = .data$estimate + 1.96*se,
                          lower = .data$estimate - 1.96*se)

  out$model_components <- model_components
  out$inv_link_funs <- inv_link_funs

  class(out) <- c("keanu_model", class(out))
  out

}

#' @export
predict.keanu_model <- function(object, newdata, na.action = na.pass) {
  mc <- prepare_model_components(forms = purrr::map(object$model_components$forms, remove_lhs),
                                 data = newdata,
                                 predvars = object$model_components$predvars,
                                 standardize_x = object$model_components$standardize_x,
                                 na.action = na.action,
                                 contrasts = object$model_components$contrasts,
                                 xlev = object$model_components$xlev,
                                 drop.unused.levels = FALSE)
  out <- get_distribution_params(df_pars = select(object$res_trans, parameter, term, beta = estimate),
                                 list_of_model_mats = mc$list_of_model_mats,
                                 inv_link_funs = object$inv_link_funs)
  as.matrix(out)
}

