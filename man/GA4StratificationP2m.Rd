\name{GA4StratificationP2m}
\alias{GA4StratificationP2m}
\title{The mutation operation in genetic algorithm for the determination of stratum boundaries and sample sizes of each stratum in stratified sampling}
\description{This is the mutation operation in the Genetic Algorithm approach for the determination of the stratum boundaries and sample sizes in each stratum in stratified sampling with an Proportional Sample Allocation Scheme.}
\usage{
GA4StratificationP2m(mutationGeneration, mutationRate, it)
}
\arguments{
  \item{mutationGeneration}{The generation that will be crossovered and transfered to the next generation.}
  \item{mutationRate}{A numeric: The mutation rate in the Genetic Algorithm process. Mutation rate must be in between 0 and 1, inclusive. Small levels of mutation rate is preferable in Genetic Algorithm approach.}
  \item{it}{An integer: The number of the iteration}
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
\code{\link{GA4StratificationP2x}}

\code{\link{GA4StratificationSelection}}
\code{\link{randomnumGenerator}}
}

\examples{
## The function is currently defined as
function(mutationGeneration,mutationRate,it)
   {

   	rowOfMutationGeneration=nrow(mutationGeneration)
   	for ( i in 1:rowOfMutationGeneration )
   	{
   	   for ( k in 1:5 )
  	   {
    		if ( runif(1,0,1) < mutationRate )
 			{
				if ( runif(1,0,1) < mutationRate )
				{
					if ( it<50 )
			 		{
	               		ones=which(mutationGeneration[i,]==1,arr.ind=TRUE)
      	         		zeros=which(mutationGeneration[i,]==0,arr.ind=TRUE)
            	   		mutationPoint=ones[sample((length(ones)-1),1)]
	               		mutationPoint1=zeros[sample(length(zeros),1)]
	               		mutationGeneration[i,mutationPoint]=0
	               		mutationGeneration[i,mutationPoint1]=1
            	 	} else
			 		{ 
	               		ones=which(mutationGeneration[i,]==1,arr.ind=TRUE)
            	   		mutationPoint=ones[sample((length(ones)-1),1)]
                   		if ( runif(1,0,1)<0.51 )
		             	{
							if ( mutationGeneration[i,(mutationPoint+1)]==0 )
							{
								mutationGeneration[i,mutationPoint]=0
		    		            mutationGeneration[i,(mutationPoint+1)]=1
              		       	}
				 	     } else if ( mutationPoint>1)
						 {
							if ( mutationGeneration[i,(mutationPoint-1)]==0 ) 
							{
							       mutationGeneration[i,mutationPoint]=0;
							       mutationGeneration[i,mutationPoint-1]=1;     
							}
						 }
					}
    			}
			}
  		}
   	}
   return(mutationGeneration)
   }

}
\keyword{ stratification}
