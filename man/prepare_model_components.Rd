% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/keanu-package.R
\name{prepare_model_components}
\alias{prepare_model_components}
\title{Prepare model-components}
\usage{
prepare_model_components(forms, data, predvars = NULL,
  standardize_x = FALSE, na.action = na.warn, contrasts = NULL,
  xlev = NULL, drop.unused.levels = FALSE)
}
\arguments{
\item{forms}{A list of formulas.}

\item{data}{A data.frame}

\item{predvars}{The 'predvars' attribute of the terms object, from a previous call to this
function.}

\item{standardize_x}{Either a logical specifying whether to center/scale numeric predictors in
the model.frame, or a list with names 'center' and 'scale', each of which in turn are named
numeric vectors specifying centering/scaling. Defaults}

\item{na.action}{A function specifying what to do with NAs. Defaults to \code{na.warn}.}

\item{contrasts}{Passed to \code{model.matrix}.}

\item{xlev}{Passed to \code{model.frame}.}

\item{drop.unused.levels}{Passed to \code{model.frame}.}
}
\value{
A list of components useful for building models.
}
\description{
Prepare model-components
}
