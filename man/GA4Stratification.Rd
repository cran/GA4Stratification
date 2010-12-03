\name{GA4Stratification}
\alias{GA4Stratification}
\title{A genetic algorithm approach to determine stratum boundaries and sample sizes of each stratum in stratified sampling}
\description{This is a package for Genetic Algorithm approach for the determination of the stratum boundaries and sample sizes in each stratum in stratified sampling.}
\usage{
GA4Stratification(dataName, numberOfStrata, sampleSize, iteration, GAgenerationSize, mutationRate, sampleAllocation)
}
\arguments{
  \item{dataName}{A string: The name of any type of numeric data containing the values of a univariate stratification variable.}
  \item{numberOfStrata}{An integer: The number of strata.}
  \item{sampleSize}{An integer: The total sample size.}
  \item{iteration}{An integer: The number of iterations in the Genetic Algorithm process.}
  \item{GAgenerationSize}{An integer: The number of generations in the Genetic Algorithm process.}
  \item{mutationRate}{A numeric: The mutation rate in the Genetic Algorithm process. Mutation rate must be in between 0 and 1, inclusive. Small levels of mutation rate is preferable in Genetic Algorithm approach.}
  \item{sampleAllocation}{A string: \code{sampleAllocation} defines the type of the sample allocation method, which could be either "Equal", "Proportional", "Neyman" or "GA".}
}
\references{
http://ideas.repec.org/a/eee/csdana/v52y2007i1p53-67.html

http://www.sciencedirect.com/science/article/B6V8V-4NHM520-1/2/a21e0295aa1616ff56da1ddf2c0ba1ac
}
\author{Sebnem Er, Timur Keskinturk, Charlie Daly

Maintainer: Sebnem Er <er.sebnem@gmail.com>}
\note{
This study is part of a project supported by the Scientific and Technological Research Council of Turkey (TUBITAK).
}

\examples{
data(iso2004)
GA4Stratification(iso2004,2,80,10,35,0.15,"Equal")
## The function is currently defined as
function(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate,sampleAllocation)
{
   if(sampleAllocation == "Equal")
   {
	GA4StratificationP1(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)
	
   } else if(sampleAllocation=="Proportional") {
	
	GA4StratificationP2(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)

   } else if(sampleAllocation == "Neyman") {
	
	GA4StratificationP3(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)
   
   } else if (sampleAllocation == "GA") {

	GA4StratificationP4(dataName,numberOfStrata,sampleSize,iteration,GAgenerationSize,mutationRate)

   }

  }
}
\keyword{ stratification}