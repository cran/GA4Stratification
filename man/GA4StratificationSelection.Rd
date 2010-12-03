\name{GA4StratificationSelection}
\alias{GA4StratificationSelection}
\title{The selection operation in genetic algorithm for the determination of stratum boundaries and sample sizes of each stratum in stratified sampling with Equal sample allocation scheme}
\description{This is the selection operation in the Genetic Algorithm approach for the determination of the stratum boundaries and sample sizes in each stratum in stratified sampling with an Equal Sample Allocation Scheme.}
\usage{
GA4StratificationSelection(selectionGeneration, selectionGenerationFitness)
}
\arguments{
  \item{selectionGeneration}{The generation that selection operator will be applied.}
  \item{selectionGenerationFitness}{The fitness value of the generation that selection operator will be applied.}
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
\code{\link{GA4StratificationP1fitt}}
\code{\link{GA4StratificationP1m}}
\code{\link{GA4StratificationP1x}}

\code{\link{randomnumGenerator}}
}

\examples{
## The function is currently defined as
function(selectionGeneration,selectionGenerationFitness)
   {
   	rowSelectionGeneration=nrow(selectionGeneration)
	colSelectionGeneration=ncol(selectionGeneration)
   	selectionStrata=array(0,dim=c(rowSelectionGeneration,(colSelectionGeneration+1)))
   	newSelectionGeneration=cbind(selectionGeneration, -selectionGenerationFitness)
   	sortedSelectionGeneration=newSelectionGeneration[order(newSelectionGeneration[,colSelectionGeneration+1]),]
   	sortedSelectionGenerationFitness=selectionGenerationFitness[order(-selectionGenerationFitness),]
   	wheelOld = sortedSelectionGenerationFitness / sum(selectionGenerationFitness)
	wheel=1:rowSelectionGeneration

	for(i in rowSelectionGeneration:1)
	{
		wheel[rowSelectionGeneration+1-i]=wheelOld[i]
	}
	wheel=cumsum(wheel)

   	for ( i in  1:rowSelectionGeneration )
   	{
		r = runif(1,0,1)
		for ( j in  1:rowSelectionGeneration )
		{
			if(r < wheel[j])
			{
				selectionStrata[i,] = sortedSelectionGeneration[j,]
				break;
			}
		}
	}
	randomGeneration=selectionStrata[,1:colSelectionGeneration]
	fitnessValueGeneration=selectionStrata[,(colSelectionGeneration+1)]
	return(randomGeneration)
   }
}
\keyword{ stratification}
