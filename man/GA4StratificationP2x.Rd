\name{GA4StratificationP2x}
\alias{GA4StratificationP2x}
\title{The crossover operation in genetic algorithm for the determination of stratum boundaries and sample sizes of each stratum in stratified sampling}
\description{This is the crossover operation in the Genetic Algorithm approach for the determination of the stratum boundaries and sample sizes in each stratum in stratified sampling with an Proportional Sample Allocation Scheme.}
\usage{
GA4StratificationP2x(crossoverGeneration,bestGeneration,dataName,numberOfStrata,sampleSize,fitnessValueGeneration,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)

}
\arguments{
  \item{crossoverGeneration}{The generation that will be crossovered and transfered to the next generation.}
  \item{bestGeneration}{The generation that has the best fitness value after the selection process.}
  \item{dataName}{A string: The name of any type of numeric data containing the values of a univariate stratification variable.}
  \item{numberOfStrata}{An integer: The number of strata.}
  \item{sampleSize}{An integer: The total sample size.}
  \item{fitnessValueGeneration}{The fitness value -the variance of the estimate- of the best generation after selection}
  \item{cumTotal}{An array: The cumulative total of the data elements from i=1 to i=N}
  \item{sumSquares}{An array: The cumulative total of the squares of the data elements from i=1 to i=N}
  \item{lengthData}{An integer: The size of the data.}
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
\code{\link{GA4StratificationP2}}
\code{\link{GA4StratificationP2fit}}
\code{\link{GA4StratificationP2fitt}}
\code{\link{GA4StratificationP2m}}

\code{\link{GA4StratificationSelection}}
\code{\link{randomnumGenerator}}
}

\examples{
## The function is currently defined as
function(crossoverGeneration,bestGeneration,dataName,numberOfStrata,sampleSize,fitnessValueGeneration,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
 {
   fitnessValueParents=fitnessValueGeneration
   parents=cbind(crossoverGeneration,fitnessValueParents)
   crossoverGenerationp=crossoverGeneration
   rowCrossoverGenerationp=nocrom

	tableData=as.data.frame(table(dataName))
 	randomnumRange=cumsum(tableData[,2])
	lengthRandomnum=length(randomnumRange)


   randomNumbers=array(0,dim=c(rowCrossoverGenerationp,2))
   
   for (i in 1:rowCrossoverGenerationp)
   {
	randomNumbers[i,]=randomnumGenerator((1:rowCrossoverGenerationp),(rowCrossoverGenerationp+1),2)
   }
   mother=father=NULL
   for (i in 1:rowCrossoverGenerationp)
   {
	  mother=randomNumbers[i,1]
	  father=randomNumbers[i,2]
   
    	  crossoverPoint=sample((randomnumRange[1:lengthRandomnum-1]),1)
        while ( sum(crossoverGenerationp[mother,c(1:crossoverPoint)])!=sum(crossoverGenerationp[father,c(1:crossoverPoint)]) )
        {
       	  crossoverPoint=sample((randomnumRange[1:lengthRandomnum-1]),1)
	  }
    	  crossoverGeneration[i,c(1:crossoverPoint)]=crossoverGenerationp[mother,c(1:crossoverPoint)]
        crossoverGeneration[i,c((crossoverPoint+1):lengthData)]=crossoverGenerationp[father,c((crossoverPoint+1):lengthData)]
   }
      s=GA4StratificationP2fit(crossoverGeneration,dataName,numberOfStrata,sampleSize,cumTotal,sumSquares,lengthData,dd,nocrom,fitp1,fit,N,means,s,n,vars,mas,NN,k,p,t)
      crossoverGenerationx=cbind(crossoverGeneration,s)
      GA4StratificationP2x=rbind(parents, crossoverGenerationx)
      GA4StratificationP2x=GA4StratificationP2x[order(GA4StratificationP2x[,(lengthData+1)]),]
      GA4StratificationP2x=GA4StratificationP2x[c((rowCrossoverGenerationp+1):(rowCrossoverGenerationp*2)),c(1:lengthData)]
	return(GA4StratificationP2x)
 }

}
\keyword{ stratification}
