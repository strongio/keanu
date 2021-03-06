% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modeller_from_likelihood.R
\name{parse_formulas_list}
\alias{parse_formulas_list}
\title{Parse formulas: a helper function for `modeller_template`}
\usage{
parse_formulas_list(formulas, .parameters)
}
\arguments{
\item{formulas}{A list of formulas.}

\item{.parameters}{A character-vector of the expected parameters.}
}
\value{
A list consisting of (1) a list of reformatted formulas ready to be
  passed to `prepare_model_components` and (2) a list of link/inverse-link pairs.
}
\description{
Parse formulas: a helper function for `modeller_template`
}
