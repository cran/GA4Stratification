\name{randomnumGenerator}
\alias{randomnumGenerator}
\title{Random number generator}
\description{This function generates random variables that don't repeat in a given range}
\usage{
randomnumGenerator(randomnumRange, lengthRandomnum, howManyRands)
}
\arguments{
  \item{randomnumRange}{This is the range of the random numbers

randomnumRange can be something given as: randomnumRange=1:487

or

randomnumRange=c(6, 12, 36, 67, 87, 146, 267)

 therefore length of the randonnumRange

 in the first case is 487

 in the second case is 7, defining the index of the randomnumRange to be swapped}

  \item{lengthRandomnum}{This is the length of the random number range}
  \item{howManyRands}{Number of random values}

}
\value{An array of different random numbers}
\references{Swap method of Charlie Daly}
\author{Sebnem Er}
\examples{
randomnumGenerator(c(1,3,5,7,9,15,21),7,3)

function(randomnumRange,lengthRandomnum,howManyRands)
   {

  	for (i in 1:(howManyRands))
   	{
		integer=sample(lengthRandomnum-1,1)
		tmp=randomnumRange[integer]
		randomnumRange[integer]=randomnumRange[i]
		randomnumRange[i]=tmp
   	}
	return(randomnumRange[1:howManyRands])
   
   }
}
\keyword{ different}
\keyword{ randomNumbers}
