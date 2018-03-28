# keanu

A package for easier custom model fitting in R.

For a complete introduction, read [our blog post](https://www.strong.io/blog/keanu-enter-the-model-matrix).

## Installation

	devtools::install_github('strongio/keanu')

## Example

Build your likelihood function.

	lognormal_loglik_fun <- function(Y, mu, sigma) dlnorm(x = Y, meanlog = mu, sdlog = sigma, log = TRUE)

Generate a new, custom model interface.

	lognormal_model <- modeller_from_loglik(loglik_fun = lognormal_loglik_fun)

Estimate the parameters of your model with your new interface!

	keanu_model <- lognormal_model(
	  formulas =  list( mu(outcome) ~ .,
	                    sigma(outcome, link='log') ~ . ),
	  data = df[c('outcome', predictors)]
	)