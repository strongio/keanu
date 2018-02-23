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

#' Parse a character into an expression
#'
#' @param x A character
#'
#' @return An expression
parse_char <- function(x) {
  stopifnot(is.character(x))
  parse(text = x)[[1]]
}

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
    if (attr(terms_obj,'intercept')==1) mm_entry <- list(`(Intercept)`=NULL)
    else mm_entry <- NULL
    out <- list(original_data = list(term_labels=NULL, model_frame=NULL, model_matrix=NULL),
         model_frame = list(term_labels=NULL, original_data=NULL, model_matrix=NULL),
         term_labels = list(original_data=NULL, model_frame=NULL, model_matrix=NULL),
         model_matrix = list(term_labels=mm_entry, model_frame=mm_entry, original_data=mm_entry))
    return(out)
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

  from_od_to_mf <- flip_mapping(from_mf_to_od)
  from_mf_to_tl <- flip_mapping(from_tl_to_mf)
  from_mf_to_mm <- flip_mapping(from_mm_to_mf)
  from_mf_to_tl <- flip_mapping(from_tl_to_mf)

  #
  list(
    original_data = list(
      term_labels = purrr::map(from_od_to_mf, ~purrr::flatten_chr(from_mf_to_tl[.x])),
      model_frame = from_od_to_mf,
      model_matrix = flip_mapping(from_mm_to_od)
    ),
    term_labels = list(
      original_data = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(from_mf_to_od[.x])),
      model_frame = from_tl_to_mf,
      model_matrix = purrr::map(from_tl_to_mf, ~purrr::flatten_chr(from_mf_to_mm[.x]))
    ),
    model_frame = list(
      original_data = from_mf_to_od,
      term_labels = from_mf_to_tl,
      model_matrix = from_mf_to_mm
    ),
    model_mat = list(
      original_data = from_mm_to_od,
      term_labels = purrr::map(from_mm_to_mf, ~purrr::flatten_chr(from_mf_to_tl[.x])),
      model_frame = from_mm_to_mf
    )
  )

}

response_to_rhs <- function(formula) {
  resp <- formula[[2]]
  # when parentheses are used on the LHS, interpreted differently than on RHS.
  # need to wrap in a temp-function that'll be stripped off later
  if (!is.symbol(resp)) {
    if (resp[[1]]=="(")
      resp[[1]] <- quote(.response)
    else if (as.character(resp[[1]])%in%c("+","-"))
      resp <- call('.response',resp)
  }
  out <- as.formula(unclass(quo(!!resp)))
  environment(out) <- environment(formula)
  out
}

#' Response term-indicator
#'
#' This is used to indicate a term on the RHS of a formula should be moved back to the LHS.
#'
#' @param x A variable
#'
#' @return The same variable
#' @export
.response <- function(x) {
  return(x)
}

#' Merge a list of formulae
#'
#' @param forms List of formulae/formulas
#' @param data Data.frame
#' @param include_response_on_rhs Defaults to FALSE. Should response terms be moved to the right-hand side, so that the model-frame for the merged formula includes them in (e.g.) NA-based row-deletion? If
#'
#' @return Single formula, with \code{environment(forms[[1]])}
merge_formulae <- function(forms, data, include_response_on_rhs = FALSE) {

  ## check for response:
  has_response_lgl <- purrr::map_int(forms, length)==3L
  if (include_response_on_rhs) {
    if (any(has_response_lgl))
      forms$.response <- merge_formulae(forms = purrr::map(forms[which(has_response_lgl)], response_to_rhs) ,
                                        data=data)
  } else {
    if (any(has_response_lgl[-1]))
      warning(call. = FALSE, "`forms` list includes LHS terms after first element; these will be removed.")
  }

  ## get all term-labels, convert to formula:
  list_of_term_labels <- purrr::map(purrr::map(forms, terms, data=data), attr, 'term.labels')
  term_labels_unique <- unique(purrr::flatten_chr(list_of_term_labels))
  if (length(term_labels_unique)>0) form_out <- stats::reformulate(term_labels_unique)
  else form_out <- ~1

  ## if response was moved to rhs, add attribute for removing/selecting it:
  if (include_response_on_rhs & any(has_response_lgl)) {
    attr(form_out,'response_idx_in_terms') <- which(attr(terms(form_out, data=data), 'term.labels') %in% list_of_term_labels$.response)
    attr(form_out,'response_remover') <- paste0("~ .",paste0("-", list_of_term_labels$.response, collapse=""))
    attr(form_out,'response_remover') <- as.formula(attr(form_out,'response_remover'), env = environment(forms[[1]]))
  } else {
    ## otherwise, first formula's lhs is used
    lazyeval::f_lhs(form_out) <- lazyeval::f_lhs(forms[[1]])
  }
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

  formula_with_response <- merge_formulae(forms, data, include_response_on_rhs = TRUE)
  response_idx <- attr(formula_with_response,'response_idx_in_terms')
  terms_with_response <- terms(formula_with_response)
  if (!is.null(predvars)) {
    if (is.character(predvars)) predvars <- parse(text = predvars)[[1]]
    # if they supplied predvars from a previous call, these wont have response-terms
    # if those response-terms were in the forms, we want to add them back
    if (!is.null(response_idx))
      predvars[length(predvars)+seq_along(response_idx)] <- attr(terms_with_response,'variables')[response_idx+1]
    attr(terms_with_response,'predvars') <- predvars
  }
  model_frame_with_response <-
    model.frame(terms_with_response, data = data, na.action = na.action,
                drop.unused.levels = drop.unused.levels, xlev = xlev)
  if (is.null(predvars)) {
    # we just computed predvars. let's remember them for future calls, but
    # don't remember response-terms, because they might not be present in
    # future calls (e.g., if youre predicting you dont need the response)
    predvars <- attr(attr(model_frame_with_response,'terms'), 'predvars')
    if (!is.null(response_idx)) predvars <- predvars[-(response_idx+1)]
  }

  # separate out response object: --
  if (!is.null(response_idx)) {
    formula_merged <- update(formula_with_response, attr(formula_with_response,"response_remover"))
    tm <- get_terms_mapping(data=data, model_frame = model_frame_with_response, terms_obj = terms_with_response)
    response_cols_in_mf <- purrr::flatten_chr(tm$term_labels$model_frame[attr(terms_with_response,'term.labels')[response_idx]])
    response_object <- model_frame_with_response[,response_cols_in_mf,drop=TRUE]
    if (!is.null(dim(response_object)))
      colnames(response_object) <- gsub(pattern = "\\.response", replacement = "", x = colnames(response_object))
    reponse_mf_idx <- which(colnames(model_frame_with_response) %in% response_cols_in_mf)
    model_frame_merged <- model_frame_with_response[,-reponse_mf_idx,drop=FALSE]
    attr(model_frame_merged,'terms') <- terms(formula_merged, data = data)
    attr(attr(model_frame_merged,'terms'),'predvars') <- predvars
    attr(attr(model_frame_merged,'terms'),'dataClasses') <- attr(terms(model_frame_with_response),'dataClasses')[-response_idx]
  } else {
    model_frame_merged <- model_frame_with_response
    formula_merged <- formula_with_response
    response_object <- NULL
  }
  xlevels <- .getXlevels(attr(model_frame_merged, "terms"), model_frame_merged)
  if (length(xlevels)==0) xlevels <- NULL

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
  list_of_terms <- purrr::map(forms, terms, data = data)
  model_mats <- purrr::map(list_of_terms, function(this_terms) {
    this_term_map <- term_mapping$term_labels$model_matrix
    if (length(this_term_map)>0) needed_cols <- unique(purrr::flatten_chr(this_term_map[attr(this_terms,'term.labels')]))
    else needed_cols <- c()
    if (attr(this_terms,'intercept')==1) needed_cols <- c("(Intercept)", needed_cols)
    model_matrix_merged_std[,needed_cols,drop=FALSE]
  })

  # return--
  out <- list(model_mats= model_mats,
              formula_merged = formula_merged,
              model_frame_merged = model_frame_merged,
              response_object = response_object,
              na_dropped_idx= attr(model_frame_merged,"na.action"),
              standardize_x= standardize_x,
              predvars = predvars,
              xlevels = xlevels,
              contrasts = contrasts,
              forms = purrr::map(list_of_terms, formula))
  class(out) <- c("model_components",class(out))
  out

}

#' Plot coefficients from a model
#'
#' @param x A model object
#' @param terms Terms to include in plot, defaulting to everything but the intercept. This can be a
#'   character-vector, or, more conveniently, a call to \code{dplyr::vars}, which provides a
#'   flexible semantics for including/excluding arbitrary terms.
#' @param se_multi Standard-errors will be mulitplied by this to obtain confidence-intervals.
#' @param ... Arguments to be passed to subsequent methods.
#'
#' @return A ggplot2 object.
#' @export
plot_coefficients <- function(x, terms, se_multi, ...) {
  UseMethod("plot_coefficients")
}

#' @export
plot_coefficients.default <- function(x, terms = NULL, se_multi = 1.96, ...) {

  df_est <- broom::tidy(x)

  if (!all(c('term', 'estimate', 'std.error') %in% colnames(df_est) ))
    stop(call. = FALSE,
         "No `plot_coefficients` method defined for object with classes:\n",
         paste(class(x), collapse = ", ") )

  if (is.null(terms)) if (is.null(terms)) terms <- dplyr::vars(-`(Intercept)`)
  stopifnot(is.character(terms) || inherits(terms, "col_list") || inherits(terms,"quosures") )
  if (is.character(terms)) terms <- paste0("`", terms, "`")
  terms <- dplyr::select_vars_(vars = unique(df_est$term), args = terms)

  df_est <- filter(.data = df_est, term %in% terms)

  ggplot(df_est, aes(x=term, y=estimate)) +
    geom_pointrange(aes(ymin = estimate - se_multi*std.error, ymax = estimate + se_multi*std.error)) +
    geom_hline(yintercept = 0, linetype='dashed') +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

#' Return whatever is inside the parantheses of a string.
#'
#' @param x A character-vector.
#'
#' @return Whatever is inside the parentheses. If no parantheses, return original string.
inside_parens <- function(x) {
  re <- "\\(([^()]+)\\)"
  result <- map(x, ~gsub(re, "\\1", str_extract_all(.x, re)[[1]]))
  map2_chr(x, result, function(orig, match) {if (length(match)>0) match else orig})
}

#' Take a data.frame with parameter, term, beta, unroll into a flat vector
#'
#' @param par_tbl A data.frame with parameter, term, beta
#'
#' @return A named numeric vector, named with parameter__sep__term
#' @export
unroll_par_tbl <- function(par_tbl) {
  with(par_tbl, structure(beta, names=paste0(parameter, "__sep__", term)))
}

#' Roll up a flat vector into a data.frame with parameter, term, beta
#'
#' @param unrolled_par A named numeric vector, named with parameter__sep__term
#'
#' @return A data.frame with parameter, term, beta
#' @export
roll_par_tbl <- function(unrolled_par) {
  tidyr::separate(col = "name", into = c("parameter", "term"), sep = "__sep__", data = tibble::enframe(unrolled_par, value = 'beta'))
}

#' Given coefficients (in a data.frame) and a list of model-mats, get predicted parameters for the
#' distribution
#'
#' @param par_tbl A data.frame with parameter, term, beta
#' @param model_mats A list of model-matrices
#' @param inv_link_funs A named list of inverse-link functions (names must match parameters)
#'
#' @return A data.frame where each row corresponds to the rows of the model-mats, and each column is
#'   a distribution-parameter.
#' @export
get_distribution_params <- function(par_tbl, model_mats, inv_link_funs = NULL) {
  list_of_pars <- purrr::map(.x = split(par_tbl, par_tbl$parameter),
                             .f = function(df) tibble::deframe(df[,c('term','beta'),drop=FALSE]) )
  model_mats <- model_mats[names(list_of_pars)]
  if (is.null(inv_link_funs))
    inv_link_funs <- purrr::map(list_of_pars, ~identity)
  else
    inv_link_funs <- inv_link_funs[names(list_of_pars)]

  list_of_pars_per_obs <- purrr::map(
    .x = names(list_of_pars),
    .f = function(param_name) {
      inv_link_funs[[param_name]]( as.matrix(model_mats[[param_name]]) %*% matrix( list_of_pars[[param_name]] ) )
    })
  names(list_of_pars_per_obs) <- names(list_of_pars)
  tibble::as_data_frame(purrr::map(list_of_pars_per_obs, as.numeric))
}

#' Initialize a dataframe with parameter, term, beta
#'
#' @param model_mats A list of model-matrices
#'
#' @return A data.frame with parameter, term, beta
#' @export
initialize_par_tbl <- function(model_mats) {
  nms <- purrr::map(model_mats, colnames)
  df_nested <- tibble::enframe(name = 'parameter', value = 'term', nms)
  df_unnested <- tidyr::unnest(df_nested)
  df_unnested$beta <- runif(n = nrow(df_unnested), min = -2, max = 2)
  df_unnested
}


#' Produce univariate predictions
#'
#' @param object A model object
#'
#' @return A list of data.frames, one for each column in the model.frame, each consisting of that
#'   column and the corresponding predictions.
#' @export
univariate_predictions <- function(object, ...) {
  UseMethod("univariate_predictions")
}


