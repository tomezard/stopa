\name{stochPerturb}
\alias{stochPerturb}
\alias{eigen.analysis}
\title{Calculate elasticities of long-run stochastic growth.}
\description{
Performs an elasticity analysis of long-run stochastic growth \code{lambdaS} against various descriptors of demographic rate distributions.
}
\usage{
stochPerturb(simltn)
}
\arguments{
  \item{simltn}{A simulation of a population obtained using \code{\link{projectPop}}.}
}
\details{
MUST follow projectPop, which stores all the required vectors and growth rates.  The length of sequence is adopted from the population simulation (\code{simltn}).  Even if no within-state variation is entered (i.e. exactly one matrix for each habitat state), all quantities are calculated.

Unlike deterministic analysis (when there is only one number to perturb), in stochastic perturbation analysis we need to perturb distributions of demographic rates.  Different perturbations - for example, holding the variation, mean or the variance of the rates fixed - yield different insights into the system.  Full details are given in the manuscripts cited below.  Briefly, the "Stochastic Elasticity" of lambdaS to the vital rates is the proportional change in lambdaS for a proportional change in that vital rate.  The coefficient of variation in the vital rates does not change but both mean and variance of those vital rates do.  This quantity (eS) can be disentangled (Tuljapurkar et al., 2003): eS = eSmu + eSsigma, where eSmu holds the mean constant and eSsigma holds the variance constant.

eSsigma can be further disentangled (Morris et al. 2006): eSsigma = eSigmaB + eSigmaW, where eSigmaB is the elasticity of long-run stochastic growth to between-state variability and eSigmaW is the elasticity of lambdaS to within-state variability.  The elasticity of lambdaS to within-state variability of the vital rates is the proportional change in lambdaS for a proportional change in variability within each environmental state, holding the state-specific mean rates and between-state rate variance constant. Conversely, the elasticity of lambdaS to variability between environmental states is the proportional change in lambdaS for a proportional change in the rate variance between each state, holding the within-state variance and overall rate means constant.

The habitat-stage elasticities are the proportional change in lambdaS for a proportional change in the vital rates if and only if the environment is in habitat state beta (Horvitz et al. 2005). This quantity is, in part, dependent on the frequency of the habitat state in the Markov chain (Horvitz et al. 2005), and hence we might want to rescale these quantities as contributions to habitat-stage elasticity (Horvitz et al. 2005). Habitat-stage elasticities are not output as contributions to habitat-stage elasticity, but can be readily converted (see the example below).

eigen.analysis is from the now missing demogR package by James Holland Jones at Stanford.
}

\value{
A list containing the following quantities
  \item{eD}{eD: deterministic elasticities of asymptotic growth}
  \item{eS}{ES: perturb the rate in every state by same percent, see Tuljapurkar et al. (2003).}
  \item{eSmu}{ESmu: perturb the mean only, see Tuljapurkar et al. (2003).}
  \item{eSsigma}{ESsigma: pertub the variability only, see Tuljapurkar et al. (2003).}
  \item{eSigmaW}{EsigmaW: elasticity to within-phase variability, see Morris et al. (2006).}
  \item{eSigmaB}{EsigmaB: elasticity to between-phase variability, see Morris et al. (2006).}
  \item{esBeta}{Ebeta: contribution to habitat-stage elasticity, see Horvitz et al. (2005).}
}

\author{Thomas H.G. Ezard \email{thomas.ezard@imperial.ac.uk} adapted from MATLAB code written by Horvitz, C. C., Haridas, C.V. and S. Tuljapurkar.}

\references{
If you use this functions, it'd be really great if you'd cite either:

Ezard, T. H. G., J.-M. Gaillard, M. J. Crawley, and T. Coulson. 2008. Habitat dependence and correlations between elasticities of long-term growth rates. Am. Nat., 172: 424-430.

or

Ezard, T.H.G. & Coulson, T. (2010). How sensitive are elasticities of long-run stochastic growth to how environmental variability is modelled? Ecol. Model. 221, 191.

The theory was developed in the following manuscripts:

Tuljapurkar, S., Horvitz, C. C., and J. B. Pascarella. (2003). The many growth rates and elasticities of populations in random environments. Am. Nat. 162: 489-502.

Haridas, C. V., and S. Tuljapurkar. (2005). Elasticities in variable environments: properties and implications. Am. Nat. 166: 481-495.

Horvitz, C. C., S. Tuljapurkar, and J. B. Pascarella. (2005). Plant-animal interactions in random environments: habitat stage elasticity, seed predators, and hurricanes. Ecology 86:3312-3322.

Morris, W. F., S. Tuljapurkar, C. V. Haridas, E. S. Menges, C. C. Horvitz, and C. A. Pfister. (2006). Sensitivity of the population growth rate to demographic variability within and between phases of the disturbance cycle. Ecol. Lett. 9:1331-1341.
}

\seealso{\code{\link{projectPop}}, \code{\link{ardisia}}, \code{\link{ardisia.envtl}}}

\examples{
data(ardisia)
data(ardisia.envtl)

##short to save time as an example, 
##Tuljapurkar et al. (2003) used 100,000 with a transient of 500.
ard <- projectPop(ardisia,ardisia.envtl[[1]],tS=10000,tT=1000)
ardEs <- stochPerturb(ard)
ardEs

##we might want tip express habitat-stage elasticities as contributions
##     as contributions (Horvitz et al. 2005)
lapply(ardEs$eSBeta, function(x) x/sum(x))

}
