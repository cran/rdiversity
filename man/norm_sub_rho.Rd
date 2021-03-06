% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diversity-measures.R
\name{norm_sub_rho}
\alias{norm_sub_rho}
\title{Normalised subcommunity rho diversity}
\usage{
norm_sub_rho(meta, qs)
}
\arguments{
\item{meta}{object of class \code{metacommunity}}

\item{qs}{\code{vector} of mode \code{numeric} containing \emph{q} values}
}
\value{
\code{norm_sub_rho} returns a standard output of class \code{rdiv}
}
\description{
Calculates similarity-sensitive normalised subcommunity rho diversity (the
representativeness of subcommunity \emph{j}). This measure may be calculated
for a series of orders, represented as a vector of \code{qs}.
}
\examples{
pop <- data.frame(a = c(1,3), b = c(1,1))
row.names(pop) <- paste0("sp", 1:2)
pop <- pop/sum(pop)
meta <- metacommunity(pop)

# Calculate normalised subcommunity rho diversity
norm_sub_rho(meta, 0:2)

}
\references{
R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
S. Mitchell, and L. Matthews. 2016. How to partition diversity.
arXiv 1404.6520v3:1–9.
}
