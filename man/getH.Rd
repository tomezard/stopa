\name{getH}
\alias{getH}

\title{Generate habitat transition matrix.}
\description{
\code{getH} is used to obtain the habitat transition matrix for a given sequence of environmental states.}

\usage{getH(sts)}

\arguments{
  \item{sts}{A vector of habitat states for successive years that correspond to succssive transition matrices.}
}

\details{
The state at time t is the column number and the state at time t+1 is the row number.  For example, the transition probability from state one to state two is in cell [1,2]
}

\value{A habitat transition matrix.}

\references{
Tuljapurkar, S., Horvitz, C. C., and J. B. Pascarella. (2003). The many growth rates and elasticities of populations in random environments. Am. Nat. 162: 489-502.

Haridas, C. V., and S. Tuljapurkar. (2005). Elasticities in variable environments: properties and implications. Am. Nat. 166: 481-495.
}

\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}

\seealso{\code{\link{projectPop}}, \code{\link{ardisia}}, \code{\link{ardisia.envtl}}}

\examples{
##generate a random (iid) seqeuence of 20 habitat states between 1 and 4
habStates <- ceiling(runif(20)*4)
##obtain the habitat transition matrix for this sequence of environmental states
getH(habStates)

}
			
