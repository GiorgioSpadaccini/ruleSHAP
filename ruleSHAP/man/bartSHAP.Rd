% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SHAP.R
\name{bartSHAP}
\alias{bartSHAP}
\title{Compute Shapley values for BART model}
\usage{
bartSHAP(model, data, n_test = nrow(data), ping = 1000)
}
\arguments{
\item{model}{An object of class \code{dbartsSampler} as obtained by fitting a BART with the \code{dbarts} function}

\item{data}{A data frame containing the data used to compute the shapley values. Should not contain the outcome \eqn{y}}

\item{n_test}{While using the whole dataset to compute the shapley values, one might only
want to compute the Shapley values of the first few observations. Such number can be regulated with \code{n_test}}

\item{ping}{An update on the progress of the fit will be given every \code{ping} iterations.
To receive no updates at all, set this to a higher value than the number of MCMC samples of the model.}
}
\value{
A dataframe containing the shapley values of each of the first \code{n_test} points of the dataset,
per predictor, with confidence intervals (\code{CIinf} being the lower bound and \code{CIsup} being the upper bound)
}
\description{
This function computes the marginal Shapley values for a BART model as fitted with dbarts.
}
