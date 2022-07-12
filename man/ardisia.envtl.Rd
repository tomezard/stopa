\name{ardisia.envtl}\alias{ardisia.envtl}\docType{data}\title{Environment transition matrices for Ardisia escalloindes}
\description{Environment transition matrices for Ardisia escalloindes based on three hurricane probabilities. The transitions are between the seven environmental states identified by Pascarella & Horvitz (1998) and modified by Tuljapurkar et al. (2003).  See \code{\link{ardisia}}.}


\format{A list containing the environmental transition matrices under scenarios of increased hurricane probability (P[hur] = 0.81), the historical hurricane probability (P[hur] = 0.081) and a decreased hurricane probability (P[hur] = 0.0081) in that order.}\source{
Pascarella, J.B., and C.C. Horvitz (1998) Hurricane disturbance and the population dynamics of a tropical understory shrub: megamatrix elasticity analysis. Ecology 79: 547-563.

Tuljapurkar, S., Horvitz, C. C., and J. B. Pascarella (2003) The many growth rates and elasticities of populations in random environments. Am. Nat. 162: 489-502}

\seealso{\code{\link{projectPop}}, \code{\link{stochPerturb}}, \code{\link{ardisia}}}
\examples{
data(ardisia.envtl)
str(ardisia.envtl)
}\keyword{datasets}