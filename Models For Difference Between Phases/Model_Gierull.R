dFuncGierullEq7 <- function(x, param){
  
  if (!is.numeric(param[3]) || param[3] %% 1 != 0) {
    stop("Please, set the number of looks to an integer value")
  }
  
  param[1]  <- abs(param[1])
  beta <- param[1] * cos(x - param[2])
  #
  soma1 <- 0.5 * beta * exp(param[3] * log(1 - param[1]^2)
                            + lgamma(param[3] + 0.5)
                            - (param[3] + 0.5) * log(1 - beta^2)
                            - lgamma(param[3])) / sqrt(pi)
  soma2 <- 0.5  * exp(param[3] * log(1 - param[1]^2)
                      + lgamma(param[3] - 0.5)
                      - param[3] * log(1 - beta^2)
                      - lgamma(param[3])) / pi^1.5
  soma3 <- beta * asin(beta) * exp(param[3] * log(1 - param[1]^2)
                                   + lgamma(param[3] - 0.5)
                                   + log(param[3] - 0.5)
                                   - (param[3] + 0.5)
                                   * log(1 - beta^2)
                                   - lgamma(param[3])) / pi^1.5
  y  <- 0.5 - param[3]
  vet <- rep(0, param[3] - 1)
  vet[1] <- 1
  for(i in 2: (param[3] - 1)){
    vet[i] <- y + i - 1
  }
  gammaratioln <- rep(0, param[3] - 1)
  for(i in 1: (param[3] - 1)){
    gammaratioln[i] <- log(abs(vet[i]))
  }
  soma4 <- 0
  for(i in 1: (param[3] - 1)){
    soma4 <- soma4 +  (-1)^(i+1)*(-1)^(i-1) * exp(param[3] * log(1 - param[1]^2)
                                                  + lgamma(param[3] - i) - lgamma(param[3])
                                                  + sum(gammaratioln[1: i])
                                                  + log(1 + (2 * i - 1) * beta^2)
                                                  - (i + 1) * log(1 - beta^2)) * 0.25 / pi
  }
  f <- soma1 + soma2 + soma3  + soma4
  return(f)
}

