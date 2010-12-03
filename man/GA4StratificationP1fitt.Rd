\name{GA4StratificationP1fitt}
\alias{GA4StratificationP1fitt}
\title{The genetic algorithm fitness function to calculate the variance of the estimate of each chromosome in the generation with an Equal Sample Allocation Scheme}
\description{This is the fitness function in Genetic Algorithm that  calculates the variance of the estimate according to the boundaries obtained with GA and sample sizes obtained with Equal Allocation for the final chromosome.}
\usage{
GA4StratificationP1fitt(bestGeneration, dataName, numberOfStrata, sampleSize, bestValue, cumTotal, sumSquares)
}
\arguments{
  \item{bestGeneration}{The generation that has the smallest fitness value is the best generation that will be delivered to the next step.}
  \item{dataName}{A string: The name of any type of numeric data containing the values of a univariate stratification variable.}
  \item{numberOfStrata}{An integer: The number of strata.}
  \item{sampleSize}{An integer: The total sample size.}
  \item{bestValue}{A numeric: The best fitness value that is the minimum variance of the estimate for the best generation.}
  \item{cumTotal}{An array: The cumulative total of the data elements from i=1 to i=N}
  \item{sumSquares}{An array: The cumulative total of the squares of the data elements from i=1 to i=N}
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
\seealso{
\code{\link{GA4StratificationP1}}
\code{\link{GA4StratificationP1fit}}
\code{\link{GA4StratificationP1m}}
\code{\link{GA4StratificationP1x}}

\code{\link{GA4StratificationSelection}}
\code{\link{randomnumGenerator}}
}

\examples{
## The function is currently defined as
function(bestGeneration,dataName,numberOfStrata,sampleSize,bestValue,cumTotal,sumSquares)
{
   c=nrow(dataName)
   nocrom=length(bestGeneration)/c
   bestGeneration=array(bestGeneration,dim=c(length(bestGeneration)/c,c))
   fitp1=array(0,dim=c(1,nocrom))
   fit=array(0,dim=c(nocrom,1))
   N=means=s=n=vars=mas=NN=k=p=t=array(0,dim=c(nocrom,numberOfStrata))

   dd=array(0,dim=c(nocrom,1))

   for ( i in 1:nocrom )
   {

   	mas[i,]=which(bestGeneration[i,]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5

    	for ( j in 2:numberOfStrata )
		{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
			N
	      	means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
			means
      		s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
			s
		}   

    	for ( j in 1:numberOfStrata )
		{
      	n[i,j]=(sampleSize/numberOfStrata)
        	}
  }

   return(array(c(N,n,bestValue),dim=c(numberOfStrata,3)))
}
}
\keyword{ stratification}
