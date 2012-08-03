\name{plot3d.localpoly.reg}
\alias{plot3d.localpoly.reg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{3d plot from a local polynomial fit
}
\description{
Create a 3d plot from a local polynomial fit of two covariates and a response variable.}
\usage{
plot3d.localpoly.reg(X,Y, bandwidth = 0, gridsize = 30, degree.pol = 0,
   kernel.type = "epanech", gridsurface = 30, xlab=expression(X_1),
   ylab=expression(X_2), zlab=expression(Y), theta = 30, phi = 30, 
   expand = 0.5, col = "lightblue", ltheta = 120, shade = 0.75, 
   ticktype = "detailed", pch = 16,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
n by 2 matrix with observations, rows corresponding to data points and columns correspond to covariates.
}
  \item{Y}{
vector of observed responses.}
  \item{bandwidth}{
bandwidth, vector or matrix. If 0, leave-one-out cross validation with criterion of minimum MSE is performed to select a unique bandwidth that will be used for all dimensions of X; if -1, Generalized Cross Validation is performed to select a unique bandwidth that will be used for all dimensions of X; if -2 leave-one-out cross validation for each covariate; and if -3, GCV for each covariate. It may be a vector for each dimension of the X; or a matrix of bandwidths (not to be confused with bandwidth matrix H), where each row is a vector of size 2, representing a bandwidth that changes with the location of estimation for the grid. See \link{localpoly.reg}.
}
  \item{gridsize}{
number of possible bandwidths to be searched in cross-validation. If left as \emph{default} 0,  gridsize is taken to be 5+as.integer(100/d^3). If cross-validation is not performed, it is ignored.
}
  \item{degree.pol}{
degree of the polynomial to be used in the local fit.}
  \item{kernel.type}{
kernel type, options are "box", "trun.normal", "gaussian", "epanech", "biweight", "triweight" and "triangular". "trun.normal" is a gaussian kernel truncated between -3 and 3.}
  \item{gridsurface}{number of points on each axis at which to estimate the local polynomial surface.}
  \item{xlab}{parameter for persp}
  \item{ylab}{parameter for persp}
  \item{zlab}{parameter for persp}
  \item{theta}{parameter for persp}
  \item{phi}{parameter for persp}
  \item{expand}{parameter for persp}
  \item{col}{parameter for persp}
  \item{ltheta}{parameter for persp}
  \item{shade}{parameter for persp}
  \item{ticktype}{parameter for persp}
  \item{pch}{parameter for persp}
  \item{...}{further parameters for plotting persp}
}
\details{
	Uses function "persp" to plot the estimated surface of a local polynomial fit in a nonparametric model with two covariates. The surface is estimated at points of a grid with size "gridsurface", which are evenly distributed between the minimum and maximum of the observed predictors. It also adds the observed points to the plot.
}
\value{
\item{X}{the same input matrix}
\item{Y}{the same input response vector}
\item{points}{points at which to get smoothed values}
\item{bandwidth}{bandwidth used for the polynomial fit}
\item{predicted}{matrix with the predicted values at grid points}
}
\author{
Adriano Zanin Zambom <adriano.zambom@gmail.com>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{localpoly.reg}}
}
\examples{
X = matrix(0,50,2)
X[,1] = runif(50,-2,2)
X[,2] = runif(50,-2,2)
Y = 4*sin(pi*X[,1]) + X[,2] + rnorm(50)

plot3d.localpoly.reg(X,Y, bandwidth=-2, gridsize = 15, degree.pol = 0, gridsurface=20)	


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
