\name{npmodelcheck}
\alias{npmodelcheck}
\title{Hypothesis Testing for Covariate or Group effect in Nonparametric Regression 
}
\description{
Tests the significance of a covariate or a group of covariates in a nonparametric regression based on residuals from a local polynomial fit of the remaining covariates using high dimensional one-way ANOVA.
}
\usage{
npmodelcheck(X, Y, ind_test, p = 7, degree.pol = 0, kernel.type =
 "epanech", bandwidth = "CV", gridsize = 30, dim.red = c(1, 10))
}
\arguments{
  \item{X}{
matrix with observations, rows corresponding to data points and columns correspond to covariates.
}
  \item{Y}{
vector of observed responses.}
  \item{ind_test}{
index or vector with indices of covariates to be tested.}
  \item{p}{
size of the window W_i. See Details.}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit.}
  \item{kernel.type}{
kernel type, options are "box", "trun.normal", "gaussian", "epanech",\cr 
"biweight", "triweight" and "triangular". "trun.normal" is a gaussian kernel truncated between -3 and 3.
}
  \item{bandwidth}{
bandwidth, vector or matrix of bandwidths for the local polynomial fit. If a vector of bandwidths, it must correspond to each covariate of X_{-(ind_test)}, that is, the covariates not being tested. If "CV", leave-one-out cross validation with criterion of minimum MSE is performed to select a unique bandwidth that will be used for all dimensions of X_{-(ind_test)}; if "GCV", Generalized Cross Validation is performed to select a unique bandwidth that will be used for all dimensions of X_{-(ind_test)}; if "CV2" leave-one-out cross validation for each covariate of X_{-(ind_test)}; and if "GCV2", GCV for each covariate of X_{-(ind_test)}. It can be a matrix of bandwidths (not to be confused with bandwidth matrix H), where each row is a vector of the same dimension of the columns of X_{-(ind_test)}, representing a bandwidth that changes with the location of estimation for multidimensional X. See \link{localpoly.reg}.
}
  \item{gridsize}{
number of possible bandwidths to be searched in cross-validation. If left as \emph{default} 0,  gridsize is taken to be 5+as.integer(100/d^3). If cross-validation is not performed, it is ignored.
}
  \item{dim.red}{
vector with first element indicating 1 for Sliced Inverse Regression (SIR) and 2 for Supervised Principal Components (SPC); the second element of the vector should be number of slices (if SIR), or number of principal components (if SPC). If 0, no dimension reduction is performed. See Details.}
}
\details{
	To test the significance of a single covariate, say X_j, assume that its observations X_{ij}, i = 1,...n, define the factor levels of a one-way ANOVA. To construct the ANOVA, each of these factor levels is augmented by including residuals from nearby covariate values. Specifically, cell "i" is augmented by the values of the residuals corresponding to observations X_{ij} for "i" in W_i (W_i defines the neighborhood, and has size "p"). These residuals are obtained from a local polynomial fit of the remaining covariates X_{-(j)}. Then, the test for the significance of X_j is the test for no factor effects in the high-dimensional one-way ANOVA. See references for further details.
	
	When testing the significance of a group of covariates, the window W_i is defined using the fist supervised principal component (SPC) of the covariates in that group; and the local polynomial fit uses the remaining covariates X_{-(ind_test)}.

Dimension reduction (SIR or SPC) is applied on the remaining covariates (X_{-(ind_test)}), which are used on the local polynomial fit. This reduction is used to moderate the effect of the curse of dimensionality when fitting nonparametric regression for several covariates. For SPC, the supervision is done in the following way: only covariates with p-values (from univariate "npmodelcheck" test with Y) < 0.3 can be selected to compose the principal components. If no covariate has p-value < 0.3, then the most significant covariate will be the only component. For SIR, the size of the effective dimension reduction space is selected automatically through sequential testing (see references for details).
}
\value{
\item{bandwidth}{bandwidth used for the local polynomial fit}
\item{predicted}{vector with the predicted values with the remaining covariates}
\item{p-value}{p-value of the test}
}
\references{
Zambom, A. Z. and Akritas, M. G. (2014). a) Nonparametric Lack-of-fit Testing and Consistent Variable Selection. Statistica Sinica, v. 24, pp. 1837-1858.

Zambom, A. Z. and Akritas, M. G. (2015). b) Signicance Testing and Group Variable Selection. Journal of Multivariate Analysis, v. 133, pp. 51-60.

Li, K. C. (1991). Sliced Inverse Regression for Dimension Reduction. Journal of the American Statistical Association, 86, 316-327.

Bair E., Hastie T., Paul D. and Tibshirani R. (2006). Prediction by supervised principal components. Journal of the American Statistical Association, 101, 119-137.

Zambom, A. Z. and Akritas, M. G. (2017) NonpModelCheck: An R Package for Nonparametric Lack-of-Fit Testing and Variable Selection, Journal of Statistical Software, 77(10), 1-28. \cr
doi:10.18637/jss.v077.i10
}
\author{
Adriano Zanin Zambom <adriano.zambom@gmail.com>
}

\seealso{
\code{\link{localpoly.reg},   \link{npvarselec}}
}
\examples{
X = matrix(1,100,5)

X[,1] = rnorm(100)
X[,2] = rnorm(100)
X[,3] = rnorm(100)
X[,4] = rnorm(100)
X[,5] = rnorm(100)
Y = X[,3]^3 + rnorm(100)

npmodelcheck(X, Y, 2, p = 9, degree.pol = 0, kernel.type = "trun.normal", 
bandwidth = c(0.85, 0.09, 2.5, 2.2),  dim.red = 0) # can use bandwidth = "CV"

# not run: can also try
#npmodelcheck(X, Y, 3, p = 7, degree.pol = 0, kernel.type = "trun.normal", 
#bandwidth = "CV",  dim.red = c(2,2)) 

#npmodelcheck(X, Y, c(1,2), p = 11, degree.pol = 0, kernel.type = "box", 
#bandwidth = c(0.5, 0.5, 0.5),  dim.red = c(1,10)) 

#npmodelcheck(X, Y, c(3,4), p = 5, degree.pol = 0, kernel.type = "box", 
#bandwidth = c(2.8, 2.8, 2.8),  dim.red = c(1,20))

#npmodelcheck(rnorm(100), rnorm(100), 1, p = 5, degree.pol = 1, 
#kernel.type = "box", bandwidth = .5,  dim.red = c(1,20))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.