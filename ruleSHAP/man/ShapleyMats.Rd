% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SHAP.R
\name{ShapleyMats}
\alias{ShapleyMats}
\title{Compute a matrix with Shapley values of each rule}
\usage{
ShapleyMats(data, data_test = data, Rs, Rs_test = Rs, id_mat, interactions = F)
}
\arguments{
\item{data}{a dataframe containing the data used to estimate the expectations in Shapley values}

\item{data_test}{a dataframe containing the points to compute the Shapley values of. By default,
this coincides with the datapoints used to estimate Shapley values.}

\item{Rs}{a list of as many matrices as there are rules to compute the Shapley values for.
The \eqn{j}-th element is a matrix corresponding to the \eqn{j}-th rule. Each of its columns
corresponds to the 0-1 encoding of a subrule of the \eqn{j}-th rule, as observed in the data
provided with input parameter \code{data}. Can be computed with the \link{RuleMats} function.}

\item{Rs_test}{same as Rs, but computed for the (possibly different) observations provided
from the \code{data_test} parameter.}

\item{interactions}{A logical parameter determining whether interaction Shapley values should
also be computed}

\item{RulePredMat}{A matrix with as many rows as there are rules and as many columns as there
are predictors. The \eqn{j}-th column has ones on the entries corresponding to rules that
involve the \eqn{j}-th predictor, while all remaining entries are zeroes. Can be computed
with the \link{RuleMats} function.}
}
\value{
marginal A matrix with \eqn{n\cdot p} rows and as many columns as
there are terms (both linear and rules). It is obtained by vertically
stacking matrices of \eqn{n} rows. Each submatrix focuses on the shapley values
of a different predictor: the \eqn{(i,k)}-th entry of the \eqn{j}-th of
such submatrices represents the contribution of the \eqn{k}-th term to
the Shapley value of the \eqn{i}-th datapoint for the \eqn{j}-th predictor.
The first \eqn{p} terms are the linear terms, and the remaining columns
refer to the rules.

interaction A matrix with \eqn{n\cdot p^2} rows and as many columns as
there are terms (both linear and rules). It is obtained by vertically
stacking matrices of \eqn{n \cdot p} rows. Each submatrix is in turn split
into \eqn{p} subsubmatrices which focuses on the interaction shapley values
of a different pair of predictors: the \eqn{(i,k)}-th entry of the \eqn{j}-th
subsubmatrix of the \eqn{j'}-th submatrix represents the contribution of the
\eqn{k}-th term to the Shapley value of the \eqn{i}-th datapoint for the
interaction between the \eqn{j}-th and the \eqn{j'}-th predictor.
The first \eqn{p} terms are the linear terms, and the remaining columns
refer to the rules.
}
\description{
This function computes the Shapley values of each rule of a RuleSHAP model and
arranges them in a matrix that can be used to compute Shapley values of the
model as a whole.
}
