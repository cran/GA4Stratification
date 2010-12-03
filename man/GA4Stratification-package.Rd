\name{GA4Stratification-package}
\alias{GA4Stratification-package}
\docType{package}
\title{A genetic algorithm approach to determine stratum boundaries and sample sizes of each stratum in stratified sampling}
\description{This is a package for Genetic Algorithm approach for the determination of the stratum boundaries and sample sizes in each stratum in stratified sampling developed by Keskinturk and Er in 2007.}
\details{
\tabular{ll}{
Package: \tab GA4Stratification\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-12-03\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr}
}
\author{Sebnem Er, Timur Keskinturk, Charlie Daly 

Maintainer: Sebnem Er <er.sebnem@gmail.com>
}
\note{
This study is part of a project supported by the Scientific and Technological Research Council of Turkey (TUBITAK).
}
\references{
http://ideas.repec.org/a/eee/csdana/v52y2007i1p53-67.html
}
\keyword{ stratification }
\seealso{
\code{\link{iso2004}}
\code{\link{P75}}
\code{\link{chi1}}
\code{\link{chi5}}
\code{\link{chi10}}
\code{\link{chi15}}
\code{\link{normal100_10}}
\code{\link{beta10_3}}

\code{\link{GA4StratificationP1}}
\code{\link{GA4StratificationP1fit}}
\code{\link{GA4StratificationP1fitt}}
\code{\link{GA4StratificationP1m}}
\code{\link{GA4StratificationP1x}}

\code{\link{GA4StratificationP2}}
\code{\link{GA4StratificationP2fit}}
\code{\link{GA4StratificationP2fitt}}
\code{\link{GA4StratificationP2m}}
\code{\link{GA4StratificationP2x}}

\code{\link{GA4StratificationP3}}
\code{\link{GA4StratificationP3fit}}
\code{\link{GA4StratificationP3fitt}}
\code{\link{GA4StratificationP3m}}
\code{\link{GA4StratificationP3x}}

\code{\link{GA4StratificationP4}}
\code{\link{GA4StratificationP4fit}}
\code{\link{GA4StratificationP4fitt}}
\code{\link{GA4StratificationP4m}}
\code{\link{GA4StratificationP4x}}

\code{\link{GA4StratificationSelection}}
\code{\link{randomnumGenerator}}
}
\examples{
data(chi1)
data(chi10)
data(chi15)
data(chi5)
data(iso2004)
data(normal100_10)
data(beta10_3)
data(P75)
GA4Stratification(chi1,2,80,10,35,0.15,"Equal")
GA4Stratification(chi10,3,81,10,35,0.15,"Proportional")
GA4Stratification(chi15,4,80,10,35,0.15,"Neyman")
GA4Stratification(chi5,5,80,10,35,0.15,"GA")
GA4Stratification(iso2004,6,84,10,35,0.15,"Equal")
GA4Stratification(normal100_10,2,80,10,35,0.15,"Proportional")
GA4Stratification(beta10_3,3,81,10,35,0.15,"Neyman")
GA4Stratification(P75,4,80,10,35,0.15,"GA")

}
