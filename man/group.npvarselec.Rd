\name{group.npvarselec}
\alias{group.npvarselec}
\title{Group variable selection for nonparametric regression
}
\description{
Performs group variable selection in a completely nonparametric regression model using hypothesis testing for high-dimensional one-way ANOVA and False Discovery Rate (FDR) corrections.
}
\usage{
group.npvarselec(X, Y, groups, method = "backward", p = 7, fitSPC 
    = TRUE, degree.pol = 0, kernel.type = "epanech", bandwidth = 
    "CV", gridsize = 10, dim.red = c(1, 10))
}
\arguments{
  \item{X}{
matrix with observations, rows corresponding to data points and columns correspond to covariates.
}
  \item{Y}{
vector of observed responses.}
  \item{groups}{
a variable of type "list" containing, in each item, a vector of indices of the covariates in each group.}
  \item{method}{
type of algorithm to run variable selection, options are "backward", "forward" and "forward2".}
  \item{p}{
size of the window W_i. See \link{npmodelcheck} for details.}
  \item{fitSPC}{
a logical indicating whether to use the first supervised principal component (SPC) of each group in the local polynomial fitting. See Details.
}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit.}
  \item{kernel.type}{
kernel type, options are "box", "trun.normal", "gaussian", "epanech",\cr 
"biweight", "triweight" and "triangular". "trun.normal" is a gaussian 
 kernel truncated between -3 and 3.
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
The selection procedure is based on the nonparametric test \link{npmodelcheck}, which  for testing the significance of group i, uses the residuals of the local polynomial regression of all the other covariates that are not in that group. 
When "fitSPC" is TRUE, the residuals for each test are computed based on the estimated regression curve m(S_{(-i)}), where S has d columns, each containing the first SPC of the corresponding group, and S_{(-i)} is the matrix S without column i.

Backward elimination is done by removing, at each step, the least significant group in the model if its p-value, obtained from the test \link{npmodelcheck}, is not significant according to False Discovery Rate (FDR) corrections (Benjamini and Yekutieli, 2001). The final model contains only groups that have significant p-values (based on FDR).

Forward selection is done by adding to the model, at each step, the group with the smallest p-value (when tested with all covariates that are already in the model), if when added, every group in the model is significant according to FDR corrections.

Forward2 selection is as follows: at each step, denote by Z = (Z_1, ..., Z_q) the groups in the model and by W = (W_1, ..., W_r) the groups not in the model (note that (Z,W) = X). Let p_j, j = 1,...r, be the maximum of the set of q+1 p-values obtained from testing each group (Z1,...,Z_q,W_j). Add to the model the group corresponding to the smallest p_j as long as, when added, all the p-values of the groups in the model are significant according to FDR corrections.


See also details of \link{npmodelcheck} and \link{localpoly.reg}.
}
\value{
\item{selected}{groups selected}
\item{p_values}{p-values of the tests of the selected groups}
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
\code{\link{npmodelcheck},   \link{localpoly.reg},   \link{npvarselec}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.