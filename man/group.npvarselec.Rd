\name{group.npvarselec}
\alias{group.npvarselec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Group variable selection for nonparametric regression
}
\description{
Performs group variable selection in a completely nonparametric regression model using hypothesis testing for high-dimensional one-way ANOVA and False Discovery Rate (FDR) corrections.
}
\usage{
group.npvarselec(X, Y, groups, method = "backward", p = 7, fitSPC = TRUE,
    degree.pol = 0, kernel.type = "epanech", bandwidth = 0, gridsize = 10, 
    dim.red = c(1, 10))
}
%- maybe also 'usage' for other objects documented here.
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
size of the window Wi. See \link{npmodelcheck} for details.}
  \item{fitSPC}{
a logical indicating whether to use the first supervised principal component (SPC) of each group in the local polynomial fitting. See Details.
}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit.}
  \item{kernel.type}{
kernel type, options are "box", "trun.normal", "gaussian", "epanech", "biweight", "triweight" and "triangular". "trun.normal" is a gaussian kernel truncated between -3 and 3.}
  \item{bandwidth}{
bandwidth for the local polynomial fit at each step of the elimination (or selection). Options are: 0 for leave-one-out cross validation with criterion of minimum MSE to select a unique bandwidth that will be used for all dimensions; -1 for Generalized Cross Validation to select a unique bandwidth that will be used for all dimensions; -2 for leave-one-out cross validation for each covariate; and -3 for GCV for each covariate. See \link{localpoly.reg}.
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
Zambom, A. Z. and Akritas, M. G. (2012). a) Nonparametric Model Checking and Variable Selection. arXiv 1205.6761.

Zambom, A. Z. and Akritas, M. G. (2012). b) Signicance Testing and Group Variable
Selection. arXiv 1205.6843.

Benjamini, Y. and Yekutieli, D. (2001) The control of false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188.

}
\author{
Adriano Zanin Zambom <adriano.zambom@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{npmodelcheck},   \link{localpoly.reg},   \link{npvarselec}}
}
\examples{
groups = vector("list",7)

groups[[1]] = c(3,8,5,12,14)
groups[[2]] = c(6,7,9,10)
groups[[3]] = 13
groups[[4]] = c(1,2,4,11)
groups[[5]] = c(15,16, 20)
groups[[6]] = 17
groups[[7]] = c(18,19)

X = matrix(1,100,20)
for (i in 1:20)
X[,i] = rnorm(100)

Y = X[,13]^3 + X[,7] + X[,15]^2 + X[,16] + rnorm(100)

group.npvarselec(X,Y,groups)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
