% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diversity-components.R
\name{raw_gamma}
\alias{raw_gamma}
\title{Gamma (low level diversity component)}
\usage{
raw_gamma(meta)
}
\arguments{
\item{meta}{object of class \code{metacommunity}}
}
\value{
\code{raw_gamma} returns an object of class \code{powermean}
}
\description{
Calculates the low-level diversity component necessary for calculating gamma
diversity.
}
\details{
Values generated from \code{raw_gamma()} may be input into \code{subdiv()} and
\code{metadiv()} to calculate subcommunity and metacommunity gamma diversity.
}
\examples{
pop <- cbind.data.frame(A = c(1,1), B = c(2,0), C = c(3,1))
row.names(pop) <- paste0("sp", 1:2)
pop <- pop/sum(pop)
meta <- metacommunity(pop)

# Calculate gamma component
g <- raw_gamma(meta)
subdiv(g, 1)
metadiv(g, 1)

}
\references{
R. Reeve, T. Leinster, C. Cobbold, J. Thompson, N. Brummitt,
S. Mitchell, and L. Matthews. 2016. How to partition diversity.
arXiv 1404.6520v3:1–9.
}
