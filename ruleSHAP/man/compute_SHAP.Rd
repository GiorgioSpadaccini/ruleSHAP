% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SHAP.R
\name{compute_SHAP}
\alias{compute_SHAP}
\title{Compute Shapley values for RuleSHAP model}
\usage{
compute_SHAP(
  model,
  data,
  intercept = T,
  alpha = 0.05,
  interactions = F,
  combine = NULL,
  block_size = 5000
)
}
\arguments{
\item{model}{An object of class \code{dbartsSampler} as obtained by fitting a BART with the \code{dbarts} function}

\item{data}{A data frame containing the data used to compute the Shapley values. Should not contain the outcome \eqn{y}}

\item{intercept}{Logical value denoting whether the fitted modeluses an intercept}

\item{alpha}{The significance level of the confidence intervals.}

\item{interactions}{A logical value denoting whether marginal interaction Shapley values should be computed,
on top of the overall values.}

\item{combine}{A list containing all combinations of Shapley values that should be additively
combined to explore the joint effect of two variables. For instance, if the joint effect
of variables \code{x.1} and \code{x.2} should be added to the dataframe, the list will
contain a vector \code{c('x.1','x.2')}.}

\item{block_size}{The computation of these Shapley values involves computations on large matrices which,
for computational reasons, are split into small submatrices. This integer denotes the size of a submatrix
in terms of how many rows it should contain.}
}
\value{
A dataframe containing the shapley values of each point of the dataset,
per predictor, with confidence intervals (\code{CIinf} being the lower bound and \code{CIsup} being the upper bound)
}
\description{
This function computes the marginal (interaction) Shapley values for a RuleSHAP model as fitted with the
\code{ruleSHAP} function. The computation includes confidence intervals.
}
