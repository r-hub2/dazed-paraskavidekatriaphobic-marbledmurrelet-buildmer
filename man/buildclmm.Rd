% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buildmer.r
\name{buildclmm}
\alias{buildclmm}
\title{Use \code{buildmer} to fit cumulative link mixed models using \code{clmm} from package \code{ordinal}}
\usage{
buildclmm(formula, data = NULL, buildmerControl = buildmerControl())
}
\arguments{
\item{formula}{A formula specifying both fixed and random effects using \code{lme4} syntax}

\item{data}{See the general documentation under \code{\link{buildmer-package}}}

\item{buildmerControl}{Control arguments for buildmer --- see the general documentation under \code{\link{buildmerControl}}}
}
\description{
Use \code{buildmer} to fit cumulative link mixed models using \code{clmm} from package \code{ordinal}
}
\examples{
if (requireNamespace('ordinal')) {
model <- buildclmm(SURENESS ~ PROD + (1|RESP),data=ordinal::soup,
buildmerControl=list(args=list(link='probit',threshold='equidistant')))
}
}
\seealso{
\code{\link{buildmer-package}}
}
