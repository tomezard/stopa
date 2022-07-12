\name{ardisia}\alias{ardisia}\docType{data}\title{Population transition matrices for Ardisia escallonioides}
\description{The population transition matrices for Ardisia escallonioides, published by Tuljapurkar et al. (2003) for each environmental state, which is delimited using the percentage of open canopy (from Horvitz and Pascarella 1998 with modifications in the last matrix [<5% open]).

The environmental states are 65% open, 55% open, 45% open, 35% open, 25% open, 15% open and <5% open. 

The eight life-history stages are seed, seedling, juvenile, pre-reproductive, small reproductive, medium reproductive, large reproductive and very large reproductive.}\usage{data(ardisia)}
\format{
A list containing the population transition matrices for the each of the seven environmental states (cf "habitat states" in Tuljapurkar et al. [2003]) identified by Pascarella & Horvitz 1998 in decreasing order of canopy 'openness'.}\source{
Pascarella, J.B., and C.C. Horvitz (1998) Hurricane disturbance and the population dynamics of a tropical understory shrub: megamatrix elasticity analysis. Ecology 79: 547-563.

Tuljapurkar, S., Horvitz, C. C., and J. B. Pascarella (2003) The many growth rates and elasticities of populations in random environments. Am. Nat. 162: 489-502
}

\seealso{\code{\link{projectPop}}, \code{\link{stochPerturb}}, \code{\link{ardisia.envtl}}}

\examples{
data(ardisia)
str(ardisia)

}\keyword{datasets}