\name{npvarselec}
\alias{npvarselec}
\title{Variable selection for covariates in nonparametric regression
}
\description{Performs variable selection using hypothesis tests of covariates in high-dimensional \cr 
	one-way ANOVA for a completely nonparametric regression model.}
\usage{
npvarselec(X, Y, method = "backward", p = 7, degree.pol = 0, 
    kernel.type = "epanech", bandwidth = "CV", gridsize = 10, 
    dim.red = c(1, 10))
}
\arguments{
  \item{X}{
matrix with observations, rows corresponding to data points and columns correspond to covariates.
}
  \item{Y}{
vector of observed responses.}
  \item{method}{
type of algorithm to run variable selection, options are "backward", "forward" and "forward2".}
  \item{p}{
size of the window W_i. See \link{npmodelcheck} for details.}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit.}
  \item{kernel.type}{
kernel type, options are "box", "trun.normal", "gaussian", "epanech",\cr 
"biweight", "triweight" and "triangular". "trun.normal" is a gaussian kernel truncated between -3 and 3.
}
  \item{bandwidth}{
bandwidth for the local polynomial fit at each step of the elimination (or selection). Options are: "CV" for leave-one-out cross validation with criterion of minimum MSE to select a unique bandwidth that will be used for all dimensions; "GCV" for Generalized Cross Validation to select a unique bandwidth that will be used for all dimensions; "CV2" for leave-one-out cross validation for each covariate; and "GCV2" for GCV for each covariate. See \link{localpoly.reg}.
}
  \item{gridsize}{
number of possible bandwidths to be searched in cross-validation. \emph{Default} is set to 10. If cross-validation is not performed, it is ignored.
}
  \item{dim.red}{
vector with first element indicating 1 for Sliced Inverse Regression (SIR) and 2 for Supervised Principal Components (SPC); the second element of the vector should be number of slices (if SIR), or number of principal components (if SPC). If 0, no dimension reduction is performed. This is used to moderate the curse of dimensionality in the local polynomial estimation at each step of the elimination (or selection). See \link{npmodelcheck} for details.}
}
\details{
Backward elimination is done by removing, at each step, the least significant covariate in the model if its p-value, obtained from the the test \link{npmodelcheck}, is not significant according to False Discovery Rate (FDR) corrections (Benjamini and Yekutieli, 2001). The precudere continues until all covariates left have significant p-values based on FDR.

Forward selection is done by adding to the model, at each step, the covariate with the smallest p-value (when tested with all covariates that are already in the model), if when added, every covariate in the model is significant according to FDR corrections.

Forward2 selection as follows: at each step, denote by Z = (Z_1, ..., Z_q) the covariates in the model and by W = (W_1, ..., W_r) the covariates not in the model (note that (Z,W) = X). Let p_j, j = 1,...r, be the maximum of the set of q+1 p-values obtained from testing each the covariates (Z1,...,Z_q,W_j). Add to the model the covariate corresponding to the smallest p_j as long as, when added, all the p-values of the covariates in the model are significant according to FDR corrections.

See also details of \link{npmodelcheck}.
}
\value{
\item{selected}{selected covariates}
\item{p_values}{p-values of the tests of the selected covariates}
}
\references{
Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. arXiv 1205.6761.

Benjamini, Y. and Yekutieli, D. (2001) The control of false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188.
}
\author{
Adriano Zanin Zambom <adriano.zambom@gmail.com>
}

\seealso{
\code{\link{npmodelcheck},   \link{localpoly.reg},   \link{group.npvarselec}}
}
\examples{
d = 10	
X = matrix(1,90,d)

for (i in 1:d)
   X[,i] = rnorm(90)
Y = X[,3]^3 + X[,6]^2 + sin(1/2*pi*X[,9]) + rnorm(90)

npvarselec(X, Y, method = "forward", p = 9, degree.pol = 0, 
kernel.type = "trun.normal", bandwidth = "CV")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
