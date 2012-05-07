\name{extractLMS}
\alias{extractLMS}
\title{Extracts LMS values from a gamlss object.}
\description{
  Extract LMS values from a gamlss object for solutions that transform the
  age axis according to the M-curve.
}
\usage{
    extractLMS(fit, data, sex = "M", grid = "classic",
                decimals = c(4,4,4), flatAge=NULL)
}

\arguments{
  \item{fit}{A gamlss object containing the final fit on transformed age, \code{t.age}.}
  \item{data}{A data frame containing the original data, with both \code{age} and 
  \code{t.age}}
  \item{sex}{A character vector indicating whether the fit applied to males
   \code{sex="M"} or females \code{sex="F"}. The default is \code{sex="M"}.}
  \item{grid}{A character vector indicating the desired age grid. See \code{ageGrid()} for 
  possible options. The default is a \code{grid="classic"}, a grid of 59 age points.}
  \item{decimals}{A numerical vector of length 3 indicating the number of significant 
  digits for rounding of the L, M and S curves, respectively.}
  \item{flatAge}{A scalar indicating the age beyond which the L, M and S values should
  be constant. The default (NULL) is not to flatten the curves.}
}
\details{
It is crucial that \code{t.age} in \code{data} correspond to exactly the 
same age transformation as used to fit the \code{gamlss} object.
Age grid values beyond the range of \code{data$age} produce \code{NA} in the
L, M and S values.
Parameter \code{flatAge} should be one of the values of the age grid.
}


\value{
  A data frame with rows corresponding to time points, and with the following columns:
  \code{sex},\code{x},\code{L},\code{M},\code{S}.
  }

\author{Stef van Buuren, 2010}

\examples{

#
library(AGD)
boys <- boys7482

# calculate initial M curve
data <- na.omit(boys[,1:2])
f0154  <- gamlss(hgt~cs(age,df=15,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(age,df=4,c.spar=c(-1.5,2.5)),
                 data=data,family=NO,
                 control=gamlss.control(n.cyc=3))                      

# calculate transformed age
t.age <- fitted(lm(data$age~fitted(f0154)))
t.age <- t.age - min(t.age)
data.t <- data.frame(data,t.age=t.age)

# calculate final solution
f0106r <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=3))

# extract the LMS reference table in the 'classic' age grid
nl4.hgt.boys <- extractLMS(fit = f0106r, data=data.t, grid="compact", 
                dec = c(0,2,5))
nl4.hgt.boys


# flatten the reference beyond age 20Y (not very useful in this data)
nl4.hgt.boys.flat <- extractLMS(fit = f0106r, data=data.t, flatAge=20)
nl4.hgt.boys.flat

# use log age transformation
data.t <- data.frame(data, t.age = log(data$age))
f0106rlog <- gamlss(hgt~cs(t.age,df=10,c.spar=c(-1.5,2.5)),
                 sigma.formula=~cs(t.age,df=6,c.spar=c(-1.5,2.5)),
                 data=data.t,family=NO,
                 control=gamlss.control(n.cyc=1))

nl4.hgt.boys.log <- extractLMS(fit = f0106rlog, data=data.t)
nl4.hgt.boys.log
}


\keyword{distribution} 

