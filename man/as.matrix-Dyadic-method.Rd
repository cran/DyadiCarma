% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meth_as_matrix.R
\name{as.matrix,Dyadic-method}
\alias{as.matrix,Dyadic-method}
\alias{as.matrix}
\title{Matrix representation of dyadic objects}
\usage{
\S4method{as.matrix}{Dyadic}(x)
}
\arguments{
\item{x}{\code{Dyadic}-object.}
}
\value{
The result is a \code{width*(2^height-1) x width*(2^height-1)} matrix.
}
\description{
Extracting the matrix representation
of a \code{Dyadic}-object.
}
\details{
The dyadic structure contains information about the type of matrix and its width and height.
}
\section{References}{

Kos, M., Podgórski, K., & Wu, H. (2025). Dyadic Factorization and Efficient Inversion of Sparse Positive Definite Matrices. arXiv. https://arxiv.org/abs/2505.08144
}

\examples{
#--------------------------------------------------------#
#------- Matrix representation of dyadic objects --------#
#--------------------------------------------------------#

N <- 4
k <- 3

# Construct four types of dyadic matrices with made of 1's
V <- construct(N, k, type = "vert") # vertical
H <- construct(N, k, type = "horiz") # horizontal
S <- construct(N, k, type = "symm") # symmetric
AS <- construct(N, k, type = "asymm", distr="norm") # asymmetric

# Convert the dyadic matrices to matrix format
mat_V <- as.matrix(V)
mat_H <- as.matrix(H)
mat_S <- as.matrix(S)
mat_AS <- as.matrix(AS)
}
\seealso{
\code{\link{Dyadic-class}} for the definition of the \code{Dyadic}-class;
\code{\link{dyadFac}} for the dyadic decomposition of dyadic matrices;
}
