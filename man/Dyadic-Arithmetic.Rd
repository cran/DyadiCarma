% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meth_arithmetic.R
\name{Dyadic Arithmetic}
\alias{Dyadic Arithmetic}
\alias{+,Dyadic,Dyadic-method}
\alias{-,Dyadic,Dyadic-method}
\alias{*,Dyadic,numeric-method}
\alias{*,numeric,Dyadic-method}
\alias{-,Dyadic,ANY-method}
\title{Arithmetic methods for Dyadic objects}
\value{
A Dyadic object representing the corresponding result of the
arithmetic operation.
}
\description{
Implements arithmetic operations for Dyadic objects, including negation,
addition, subtraction, and scalar multiplication.
}
\section{Methods}{

\describe{
  \item{Unary `-'}{Negates a Dyadic object.}
  \item{`+'}{Adds two Dyadic objects.}
  \item{`-'}{Subtracts one Dyadic object from another.}
  \item{`*'}{Multiplies a Dyadic object by a scalar or vice versa.}
}
}

\section{References}{

Kos, M., Podgórski, K., & Wu, H. (2025). Dyadic Factorization and Efficient Inversion of Sparse Positive Definite Matrices. arXiv. https://arxiv.org/abs/2505.08144
}

\examples{
#--------------------------------------------------------#
#-------- Arithmetic methods for dyadic objects ---------#
#--------------------------------------------------------#

N <- 4
k <- 3

# Construct four types of dyadic matrices with made of 1's
V <- construct(N, k, type = "vert") # vertical
H <- construct(N, k, type = "horiz") # horizontal
S <- construct(N, k, type = "symm") # symmetric
AS <- construct(N, k, type = "asymm") # asymmetric

# Negation of dyadic objects (matrices)
NegV <- -V
NegV@type
all(as.matrix(NegV) == -as.matrix(V)) # Should be TRUE

# Addition of dyadic objects (matrices)
HpV <- H + V # horizontal + vertical = asymmetric
HpV@type

# Subtraction of dyadic objects (matrices)
SmAS <- S - AS # symmetric - asymmetric = asymmetric
SmAS@type

# Scalar multiplication of dyadic objects (matrices)
DoubleV <- 2 * V # Scalar multiplication does not change the type
VDouble <- V * 2 # Scalar multiplication does not change the type
DoubleV@type
VDouble@type
all(as.matrix(DoubleV) == 2 * as.matrix(V)) # Should be TRUE
all(as.matrix(VDouble) == as.matrix(DoubleV)) # Should be TRUE

# Linear combination
linearComb <- -S + 3 * H - 6 * AS + V # linear combination of dyadic matrices
linearComb@type # "asymm"
}
\seealso{
\code{\link{Dyadic-class}} for the definition of the \code{Dyadic}-class;
\code{\link[DyadiCarma]{as.matrix}} for extracting the matrix representation
of a \code{Dyadic}-object
}
