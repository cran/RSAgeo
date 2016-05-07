\name{RSA}
\alias{RSA}

\title{
 A Resampling-based Stochastic Approximation Method for Analysis of Large Geostatitical Data
}
\description{
Performs parameter estimation using a resampling-based Stochastic Approximation (RSA) method.
It is a stochatic approximation method. At every iteration, only a subset of the data is drawn and used to update the estimation of the parameters.
The data are assumed to have a powered exponential correlation structure.
}
\usage{
RSA(Data, N_subset, Stepscale, Total_Iteration, Warm)
}

\arguments{
  \item{Data}{
	a Matrix. The first two column of the matrix consists of x1 and x2, where x1 and x2 provide the 2D location for each obsevation. The third column of the matrix gives the y value, which is the response value. Extra covariates can be included in Column 4 and beyond.
}
  \item{N_subset}{
        an integer. It specifies the size of the subset drawn from the data.
}
  \item{Stepscale}{
    a positive number. It controls the gain factor. It specifies the number of iterations when the gain factor begins to shrink. For example, one can be set it equal to 2 times the burn-in steps.
}
  \item{Total_Iteration}{
    an interger. It gives the total number of iterations for stochatic approximation.
}
  \item{Warm}{
    an integer. It is the number of burn-in iterations. 
}
}

\value{
  \item{beta }{gives the coefficient estimates of the mean effect. It is a vector of lenght equal to the number of coefficients plus 1.}
  \item{phi }{gives the shape estimate in the powered exponential correlation matrix.}
  \item{sigmasq }{gives the estimate of error variance.}
  \item{tausq }{gives the estimate of nugget variance.}
}
\references{
Liang, F., Cheng, Y., Song, Q., Park, J., and Yang, P. (2013) A Resampling-based Stochastic Approximation Method for Analysis of Large Geostatistical Data. J. Amer. Statist. Assoc., 108, 325-339.
}
\author{
  Yichen Cheng, Faming Liang
}

\examples{
data(gdata)
## This is an illustration of the function. 
## In practice, Total_Iteration is recommended to be set to 2500 or higher.
## N_subset is recommemded to be set to 300 or higher.
RSA(gdata,N_subset=50,Stepscale=40,Total_Iteration=100,Warm=20)
}
