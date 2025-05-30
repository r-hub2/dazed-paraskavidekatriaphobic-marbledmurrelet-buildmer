% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility.r
\name{re2mgcv}
\alias{re2mgcv}
\title{Convert lme4 random-effect terms to mgcv 're' smooths}
\usage{
re2mgcv(formula, data, drop = TRUE)
}
\arguments{
\item{formula}{The lme4 formula.}

\item{data}{The data.}

\item{drop}{Logical indicating whether constant, non-intercept columns should be dropped. Default \code{TRUE}. A warning is issued if a column needed to be dropped. Note that repeated intercept columns are silently merged without a warning.}
}
\description{
Convert lme4 random-effect terms to mgcv 're' smooths
}
\examples{
library(buildmer)
re <- re2mgcv(temp ~ angle + (1|replicate) + (1|recipe),lme4::cake)
model <- buildgam(re$formula,re$data)
# note: the below does NOT work, as the dependent variable is looked up in the data by name!
\dontshow{if (FALSE)}
re <- re2mgcv(log(Reaction) ~ Days + (Days|Subject),lme4::sleepstudy)
}
