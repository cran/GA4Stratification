\name{GA4StratificationP3fit}
\alias{GA4StratificationP3fit}
\title{The genetic algorithm (GA) fitness function that calculates the variance of the estimate of each chromosome in the GA generation with an Neyman Sample Allocation Scheme}
\description{This is the fitness function in GA that calculates the variance of the estimate according to the boundaries obtained with GA and sample sizes obtained with Neyman Allocation.}
\usage{
GA4StratificationP3fit(randomGeneration, dataName, numberOfStrata, sampleSize, cumTotal, sumSquares, c, dd, nocrom, fitp1, fit, N, means, s, n, vars, mas, NN, k, p, t)
}
\arguments{
  \item{randomGeneration}{The generation which a fitness function will be applied and a fitness value will be calculated. This is initially a random generation and after each iteration it is the mutated, crossovered and selected generation.}
  \item{dataName}{A string: The name of any type of numeric data containing the values of a univariate stratification variable.}
  \item{numberOfStrata}{An integer: The number of strata.}
  \item{sampleSize}{An integer: The total sample size.}
  \item{cumTotal}{An array: The cumulative total of the data elements from i=1 to i=N}
  \item{sumSquares}{An array: The cumulative total of the squares of the data elements from i=1 to i=N}
  \item{c}{An integer: The length of the data.}
  \item{dd}{An array (nocrom X 1): The minimum of the difference in between the size of each stratum and the sample size to be drawn from that stratum.}
  \item{nocrom}{An integer: The number of the chromosomes in the generation.}
  \item{fitp1}{An array (1 X nocrom): The fitness value for each chromosome in the generation.}
  \item{fit}{An array (nocrom X 1): The fitness value for each chromosome in the generation.}
  \item{N}{An array (nocrom X nofstrata): The number of the elements in each stratum for each chromosome in the generation.}
  \item{means}{An array (nocrom X nofstrata): The mean of the elements in each stratum for each chromosome in the generation.}
  \item{s}{An array (nocrom X nofstrata): The standart deviation of the elements in each stratum for each chromosome in the generation.}
  \item{n}{An array (nocrom X nofstrata): The number of the sample size to be drawn from each stratum for each chromosome in the generation.}
  \item{vars}{An array (nocrom X nofstrata): The variance of the estimate in each stratum for each chromosome in the generation.}
  \item{mas}{An array (nocrom X nofstrata): The indice of each stratum for each chromosome in the generation. }
  \item{NN}{An array (nocrom X nofstrata): The cumulative sum of the number of the elements in each stratum.}
  \item{k}{An array (nocrom X nofstrata): The difference between the number of the elements and the sample sizes in each stratum.}
  \item{t}{An array (nocrom X nofstrata): The maximum of the \code{k}}
  \item{p}{An array (nocrom X nofstrata): The indice of the element where \code{k} is equal to \code{t}}
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
\code{\link{GA4StratificationP3}}
\code{\link{GA4StratificationP3fitt}}
\code{\link{GA4StratificationP3m}}
\code{\link{GA4StratificationP3x}}

\code{\link{GA4StratificationSelection}}
\code{\link{randomnumGenerator}}
}

\examples{
## The function is currently defined as
function(randomGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,c,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
{

   for ( i in 1:nocrom )
   {
   	mas[i,]=which(randomGeneration[i,]==1,arr.ind=TRUE)
    	N[i,1]=min(mas[i,])
    	means[i,1]=cumTotal[mas[i,1],]/N[i,1]
	if(N[i,1]==1)
	{s[i,1]=0
	} else
    	s[i,1]=((N[i,1]/(N[i,1]-1))*(sumSquares[N[i,1]]/N[i,1]-means[i,1]^2))^.5

    	for ( j in 2:numberOfStrata )
		{    
    		N[i,j]=mas[i,j]-mas[i,(j-1)]
			
	      	means[i,j]=(cumTotal[mas[i,j],]-cumTotal[mas[i,j-1],])/N[i,j]
				if(N[i,j]==1)
				{s[i,j]=0
				} else
      		s[i,j]=((N[i,j]/(N[i,j]-1))*((sumSquares[mas[i,j]]-sumSquares[mas[i,j-1]])/N[i,j]-means[i,j]^2))^.5
			
		}   

    	for ( j in 1:numberOfStrata )
		{
      		n[i,j]=max(1,round(sampleSize*N[i,j]*s[i,j]/sum(N[i,]*s[i,])))
        		n[i,j]=min(n[i,j],N[i,j])
			
		}
		if ( (sampleSize-sum(n[i,]))>0 )
		{
      		k[i,]=N[i,]-n[i,]
        		t[i,]=max(k[i,])
        		p[i,]=which(k[i,]==t[i,],arr.ind=TRUE)[1]
        		n[i,p[i,1]]=min(n[i,p[i,1]]+sampleSize-sum(n[i,]),N[i,p[i,1]])
		}

		for ( j in 1:numberOfStrata )
		{
        		vars[i,j]=((N[i,j]-n[i,j])*s[i,j]^2*N[i,j]^2)/(c^2*n[i,j]*N[i,j])
		}

    	dd[i,]=min((N[i,]-n[i,]))
    	NN[i,]=cumsum(N[i,])
    	kl=0

   		fit[i,]=sum(vars[i,])

    	if ( dd[i]<0 ) 
		{
			fit[i]= 9999999999999999
    	} else if (!all(N[i,]!=1)) 
		{
			fit[i]= 999999999999999999
    	} else if (!all(N[i,]!=0)) 
		{
			fit[i]= 999999999999999999
		} else
		{
			fit[i]=fit[i]
		}
   
    	for ( j in 1:(numberOfStrata-1) )
		{
        	kl=kl+dataName[(NN[i,j]+1),]-dataName[NN[i,j],]
		}

    	p3fit=array(-fit,dim=c(nocrom,1))

   }
   return(p3fit)
}

}
\keyword{ stratification}
