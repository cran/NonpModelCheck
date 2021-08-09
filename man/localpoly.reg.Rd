\name{localpoly.reg}
\alias{localpoly.reg}
\title{Local Polynomial Regression Fitting
}
\description{
Computes the smoothed response or its derivatives in a nonparametric regression using local polynomial fitting.}
\usage{
localpoly.reg(X, Y, points = NULL, bandwidth = "CV", 
gridsize = 30, degree.pol = 0, kernel.type = "epanech", 
deriv = 0)}
\arguments{
  \item{X}{matrix with observations, rows corresponding to data points and columns corresponding to covariates.}
  \item{Y}{
vector of observed responses.}
  \item{points}{points at which to get smoothed values. If NULL, estimation is done on the observations of X.}
  \item{bandwidth}{
bandwidth, vector, matrix of bandwidths, "CV", "GCV", "CV2", "GCV2" or "Adp". 
When X is univariate(vector): options for the bandwidth are: "CV" for leave-one-out cross validation with criterion of minimum MSE; "GCV" for Generalized Cross Validation; "Adp" for adaptive bandwidth (see details); a positive vector of same length as X representing a bandwidth that changes with the location of estimation.
When X is multivariate(matrix): options for the bandwidth are: "CV" for leave-one-out cross validation with criterion of minimum MSE (search is done in a grid margilally for each covariate); "GCV" for Generalized Cross Validation (search is done in a grid marginally for each covariate); "CV2" for leave-one-out cross validation (search is done in any combination of grids for each covariate); "GCV2" for GCV for each covariate(search is done in any combination of grids for each covariate); or a vector or matrix of the same size as 'points';  
See \emph{Details}.
}
  \item{gridsize}{
number of possible bandwidths to be searched in cross-validation. If cross-validation is not performed, it is ignored.
}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit. In the univariate case there is no restriction; in the multivariate case, the degree can be 0,1 or 2.}
  \item{kernel.type}{kernel type, options are "box", "trun.normal", "gaussian", "epanech",\cr 
  	"biweight", "triweight" and "triangular". "trun.normal" is a gaussian kernel truncated between -3 and 3.
  	}
  \item{deriv}{order of the derivative of the regression function to be estimated.}
}
\details{
Computes smoothed values using local polynomial fitting with the specified kernel type. If multidimensional, a multiplicative(product) kernel is used as weight.

In cross validation, for multivariate X and bandwidth options "CV" and "GCV", the procedure searches individually for each covariate, the bandwidth that produces the smallest MSE from a grid of \emph{gridsize} possible bandwidths evenly distributed between the minimum and maximum/2 distance of any of the points of that covariate. In other words, one separate cross-validation is performed in each dimension of X.

In cross validation, for multivariate X and bandwidth options "CV2" and "GCV2", a d-dimensional (number of covariates) grid is created, where each dimension of the grid is a vector of \emph{gridsize} possible bandwidths evenly distributed between the minimum and maximum/2 distance of any of the points in each covariate. Then a search is done crossing all possible combinations of values of each dimension of the grid, where the resulting vector of bandwidths correspond to those which yeild minimum MSE.

Adaptive bandwidth, for univariate X only, is obtained by a similar procedure to the one proposed by Fan and Gijbels (1995). The interval is split into [1.5*n/(10*log(n))] intervals, a leave-one-out cross validation is performed in each interval to obtain a local bandwidth. These bandwidths are then smoothed to obtain the bandwidth for each point in X.
}
\value{
\item{X}{the same input matrix}
\item{Y}{the same input response vector}
\item{points}{points at which smoothed values were computed}
\item{bandwidth}{bandwidth used for the polynomial fit}
\item{predicted}{vector with the predicted(smoothed) values}
}
\references{
Fan J. and Gijbels I. (1995). Data-driven bandwidth selection in local polynomial fitting: Variable Bandwidth and Spatial Adaptation. JRSS-B. Vol 57(2), 371-394.

Wand M. P. and Jones M. C. (1995). Kernel Smoothing. Chapman and Hall.
}
\author{
Adriano Zanin Zambom <adriano.zambom@gmail.com>
}

\seealso{
\code{\link{npvarselec},   \link{npmodelcheck}}
}
\examples{
X = rnorm(100)
Y = X^3 + rnorm(100)

localpoly.reg(X, Y, degree.pol = 0, kernel.type = "box",
 bandwidth = "CV")
localpoly.reg(X, Y, degree.pol = 1, kernel.type = "box",
 bandwidth = "CV")
##--
X = runif(100,-3,3)
Y = sin(1/2*pi*X) + rnorm(100,0,.5)

localpoly.reg(X, Y, degree.pol = 0, kernel.type = "gaussian",
  bandwidth = "CV")
localpoly.reg(X, Y, degree.pol = 1, kernel.type = "gaussian", 
  bandwidth = "CV")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.