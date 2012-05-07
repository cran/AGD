\name{wp.twin}
\alias{wp.twin}
\title{Superposes two worm plots}
\description{
Superposes two worm plots from GAMLSS fitted objects. This is 
a diagnostic tool for comparing two solutions.
          }
\usage{
wp.twin(obj1, obj2 = NULL, xvar = NULL, xvar.column = 2, n.inter = 16, 
  show.given = FALSE, ylim.worm = 0.5, line = FALSE, 
  cex = 1, col1 = "black", col2 = "orange", warnings = FALSE, ...) 
}

\arguments{
  \item{obj1}{a GAMLSS fitted object}
  \item{obj2}{an optional second GAMLSS fitted object}
  \item{xvar}{ the explanatory variable against which the worm plots will be plotted }
  \item{xvar.column}{ the number referring to the column of \code{obj1$mu.x} and \code{obj2$mu.x}. 
  If \code{xvar=NULL} then the explanatory variable is set to 
  \code{xvar=obj1$mu.x[,xvar.column]} respectively \code{xvar=obj2$mu.x[,xvar.column]}.
  The default is \code{xvar.column=2}, which selects the variable following the intercept 
  (which is typically age in most applications).}
  \item{n.inter}{the number of intervals in which the explanatory variable \code{xvar} will be cut. The default is 16.}
  \item{show.given}{whether to show the x-variable intervals in the top of the graph, default is \code{show.given=FALSE} }
  \item{ylim.worm}{for multiple plots, this values is the y-variable limit, default value is \code{ylim.worm=0.5}}
  \item{line}{whether to plot the polynomial line in the worm plot, default value is \code{line=FALSE}}
  \item{cex}{ the cex plotting parameter with default \code{cex=1}}
  \item{col1}{ the color for the points of \code{obj1}. The default \code{col="black"} }
  \item{col2}{ the color for the points of \code{obj2}. The default \code{col="orange"} }
  \item{warnings}{ a logical indicating whether warnings should be produced. The default \code{warnings=FALSE} }
   \item{\dots}{for extra arguments, \code{overlap}, \code{xlim.worm} or \code{pch}}
 }
\details{
 This function is a customized version of the \code{wp()} function found in the \code{gamlss} package. Function
 \code{wp.twin()} allows overplotting of two worm plots, each in its own color. 
 The points of \code{obj1} are plotted first, the points of \code{obj2} are superposed.
 This twin worm plot provide a visual assessment of  the differences between the solutions. 
 Extra arguments can be specified (e.g. \code{xvar}) that are passed down to the \code{wp()} function of \code{gamlss} 
 if specified.
 The worm plot is a detrended normal QQ-plot that highlight departures from normality.

 Argument \code{xvar} takes priority over \code{xvar.column}. The \code{xvar} variable is cut into \code{n.iter}
 intervals with an equal number observations and detrended normal QQ (i.e. worm) plots for each interval are plotted. 
 This is a way of highlighting failures of the model within different ranges of the 
 explanatory variable. 
 
 If \code{line=TRUE} and \code{n.inter>1}, the fitted coefficients from fitting cubic polynomials to the residuals (within each x-variable interval) can be obtain by e.g. 
 \code{coeffs<-wp.twin(model1,xvar=x,n.iner=9)}.  van Buuren \emph{et al.} (2001) used these residuals to identify regions (intervals) of the 
 explanatory variable within which the model does not fit adequately the data (called "model violation")  
 
 
}
\value{
  For multiple plots the \code{xvar} intervals and the coefficients of the fitted cubic 
  polynomials to the residuals (within each \code{xvar} interval) are returned.   
}
\references{
Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
           
van Buuren and Fredriks M. (2001) Worm plot: simple diagnostic device for modelling growth reference curves. 
            \emph{Statistics in Medicine}, \bold{20}, 1259--1277.
            \url{http://www.stefvanbuuren.nl/publications/Worm plot - Stat Med 2001.pdf}

van Buuren and Fredriks M. (2007) Worm plot to diagnose fit in quantile regression. 
            \emph{Statistical Modelling}, \bold{7}, 4, 363--376.
            \url{http://www.stefvanbuuren.nl/publications/WP quantile regression - Stat Mod 2007.pdf}
}
            
\author{Stef van Buuren, using R code of Mikis Stasinopoulos and Bob Rigby}

\seealso{  \code{\link{wp}}}

\examples{
data(abdom)
a <- gamlss(y~cs(x,df=1),sigma.fo=~cs(x,0),family=LO,data=abdom)
b <- gamlss(y~cs(x,df=3),sigma.fo=~cs(x,1),family=LO,data=abdom)
coeff1 <- wp.twin(a,b,line=TRUE)
coeff1
rm(a,b,coeff1)
}

\keyword{smooth}

