eigen.analysis <- function (A) 
{
    ev <- eigen(A)
    lmax <- which(Re(ev$values) == max(Re(ev$values)))
    lambda <- Re(ev$values[lmax])
    W <- ev$vectors
    w <- abs(Re(W[, lmax]))
    V <- Conj(solve(W))
    v <- abs(Re(V[lmax, ]))
    s <- v %o% w
    s[A == 0] <- 0
    class(s) <- "leslie.matrix"
    e <- s * A/lambda
    rho <- lambda/abs(Re(ev$values[2]))
    eigen.analysis <- list(lambda1 = lambda, rho = rho, sensitivities = s, 
        elasticities = e, stable.age = w/sum(w), repro.value = v/v[1])
    eigen.analysis
}

