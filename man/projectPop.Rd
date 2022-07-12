\name{projectPop}
\alias{projectPop}

\title{Project a population forward.}
\description{
Project a population forward in time based on a sequence of environmental states with corresponding demographic rates.}

\usage{
projectPop(demogT, envtlT, sts=NULL, tS = 100000, tT = 10000)
}

\arguments{
  \item{demogT}{A \code{\link{list}} of population transition matrices, containing demographic rates for each environmental state.}
  \item{envtlT}{The environment transition matrix, containing the transition rates between environmental states.}
  \item{sts}{If there are as many environmental states as years, then elasticities of within- and between-state variance can be calculated.  This vector contains the codes to associate the environmental states with the individual years.}
  \item{tS}{the length of the environmental sequence to be modelled.}
  \item{tT}{the length of the transient to be discarded.}
}
\details{
Based on a given environmental transition matrix, the focal population is projected through a stochastic environment. The environmental state at time (t+1) is probabilistically dependent upon the environmental state at time t and defined by the Markov chain.  This process was iterated to project the population for \code{tS} time steps with a transient of length \code{tT} discarded from either end.  Using the sequence of environmental states that correspond to a set of demographic rates, each transition matrix is obtained in the following manner:

(1) obtain the habitat state at time t from the Markov chain;
(2) randomly select a year when the population was in that habitat state;
(3) fill A(t) using the observed demographic rates for that year.

These processes generate a transition matrix at each time, which was stored with associated vectors of population structure and reproductive value for perturbation analysis in stochPerturb.  Note that simulated, not asymptotic, sequences of population structure and reproductive value are stored, i.e. the population structure does not converge to the stable-age distribution.

}

\value{
A list containing the salient points from the population projection:
  \item{growth}{The one-step population growth rate (see Tuljapurkar et al. 2003).}
  \item{popStruct}{A matrix of population structure vectors (dominant left eigenvector).}
  \item{repVal}{A matrix of reproducitve value vectors (dominant right eigenvector).}
  \item{rates}{The population transition matrices.}
  \item{states}{The habitat state at each time step.}
}

\references{
If you use this functions, it'd be really great if you'd cite either:

Ezard, T. H. G., J.-M. Gaillard, M. J. Crawley, and T. Coulson. 2008. Habitat dependence and correlations between elasticities of long-term growth rates. Am. Nat., 172: 424-430.

or

Ezard, T. H. G., and T. Coulson. in resubmission. How sensitive are elasticities of long-run stochastic growth to how environmental variability is modelled? 


The theory was developed in the following manuscripts:

Tuljapurkar, S., Horvitz, C. C., and J. B. Pascarella. (2003). The many growth rates and elasticities of populations in random environments. Am. Nat. 162: 489-502.

Haridas, C. V., and S. Tuljapurkar. (2005). Elasticities in variable environments: properties and implications. Am. Nat. 166: 481-495.
}

\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk}}
\seealso{\code{\link{stochPerturb}}, \code{\link{ardisia}}, \code{\link{ardisia.envtl}}}

\examples{
data(ardisia)
data(ardisia.envtl)

##really short to save time as an example
####Tuljapurkar et al. (2003) used 100,000 with a transient of 500
ard <- projectPop(ardisia,ardisia.envtl[[1]],tS=1000,tT=100)

##plot a section to illustrate how the states change
plot(ard$states[333:444],type='b',xlab="Time",ylab="Environmental State")

}

