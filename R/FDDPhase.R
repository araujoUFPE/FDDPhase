require(roxygen2)

################################################################################
## FDDPhase: A package for statistical analysis of phase data under           ##
##          univariate physical and empirical models.                         ##
################################################################################

############################## dFuncLeeEq5 #####################################
#' @title Probability density function
#' @name dFuncLeeEq5
#'
#' @description Article "Intensity and phase statistics of multilook polarimetric
#'              and interferometric SAR imagery", Lee et al. (1994), Eq. 18.
#'
#' @param x  Equivalent to \code{sample=psi}, is the array describing the support of the
#'           phase at which the pdf  is computed.
#'
#' @param r  Complex correction coefficient (coherence), \code{0<r<=1}.
#'
#' @param theta Is a potential mean between, \code{-pi<theta<=pi}.
#'
#' @param L  Number of looks or target, \code{L>=1} positive real.
#'
#' @details Multilook probability density function for phase difference
#'          involving the Hypergeometric Gauss Function (GHF).
#'
#' @return A graphic
#'
#' @importFrom hypergeo hypergeo
#'
#' @author Araújoo, J. E.
#'
#' @examples
#' sample <- seq(-pi,pi, 0.001)
#' param <- c(r = 0.7, theta = 0, L = 4)
#' data<-data.frame(x=sample,y=dFuncLeeEq5(sample, param))
#' ggplot(data,aes(x,y)) +
#'   geom_line() +
#'   xlim(-pi, pi) +
#'   # coord_polar() +
#'   xlab("sample") +
#'   ylab("density") +
#'   theme_pander()
#' @export
dFuncLeeEq5 <- function(x, param){
  r2L <- (1-(param[1])^2)^param[3]
  rcost <- param[1]*cos(x-param[2])

  # Test the condition for the indicator function on both the support and the
  # parameter space
  ifelse(x > -pi & x <= pi & param[1] > 0 & param[1] <= 1 & param[2] > -pi & param[2] <= pi & param[3] >= 1,
         # yes
         return(Re(
           (gamma(param[3]+0.5)/(2*sqrt(pi)*gamma(param[3]))) *
             (r2L/(1-rcost^2)^(param[3]+0.5)) * rcost +
             r2L/(2*pi) * hypergeo(param[3], 1, 0.5, rcost^2))),
         # no
         return(0))
}

############################# dFuncGierullEq7 ##################################

#' @title Probability density function
#' @name dFuncGierullEq7
#'
#' @description This function calculates the pdf of the interferometric phase.
#'              Article "Closed-form expressions for InSAR sample statistics and
#'              its application to non-gaussian data" Gierull (2020), Eq. 13.
#'
#' @param x Equivalent to \code{sample=psi}, is the array describing the support of the
#'          phase at which the pdf  is computed.
#'
#' @param r Denotes the magnitude of the complex correlation (coherence),
#'          \code{0<r<=1}.
#' @param theta  Is a potential mean between, \code{-pi<theta<=pi}.
#'
#' @param L Number of looks or targeted, \code{L>1} integer.
#'
#' @details The correct phase pdf is derived using identity tables for GHF.
#'
#' @return A graphic
#'
#' @author Araújo, J. E.
#'
#' @examples
#' sample = seq(-pi,pi, 0.001)
#' param <- c(r = 0.9, theta = 0, L = 4)
#' data<-data.frame(x=sample,y=dFuncGierullEq7(sample, param))
#' ggplot(data,aes(x,y)) +
#'   geom_line() +
#'   xlim(-pi, pi) +
#'   # coord_polar() +
#'   xlab("sample") +
#'   ylab("density") +
#'   theme_pander()
#' @export
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

############################### Cauchy Truncated ###############################
#' @title Probability density function
#' @name dTruncCauchy
#'
#' @description The truncated Cauchy is a Cauchy distribution bounded between
#'              low a nd high (the pdf is 0 outside these bounds and renormalized).
#'              Samples from this distribution are differentiable with respect to
#'              loc and scale, but not with respect to the bounds low and high
#'
#' @param x  Equivalent to \code{sample=psi}, is the array describing the support of the
#'           phase at which the pdf.
#' @param mu  loc: Floating point tensor; the modes of the corresponding
#'           non-truncated Cauchy distribution(s).
#' @param sigma  scale: Floating point tensor; the scales of the distribution(s).
#'              Must contain only positive values.
#'
#' @details The truncated Cauchy distribution has the same density, with a
#'          different proportionality constant, over a finite interval. It's a
#'          well-behaved distribution, with a mean, variance, and any other
#'          moment you might care about.
#'
#' @return A graphic
#'
#' @importFrom stats dcauchy
#'
#' @importFrom stats pcauchy
#'
#' @author Araújo, J. E.
#'
#' @examples
#' sample = seq(-pi,pi, 0.001)
#' param <- c(mu = 0.5, sigma = 0.5)
#' data<-data.frame(x=sample,y=dTruncCauchy(sample, param))
#' ggplot(data,aes(x,y)) +
#'  geom_line() +
#'  xlim(-pi, pi) +
#'  # coord_polar() +
#'  xlab("sample") +
#'  ylab("density") +
#'  theme_pander()
#'
##Empirical model
##Distribution Cauchy truncated
#' @export
#'
dTruncCauchy <- function(x, param){
  # Transform the value of x to the interval (-π, π), making the function periodic.
  x <- (x - 2 * pi * floor((x + pi) / (2 * pi))) #is the mod
  c <- pcauchy((pi - param[1])/param[2]) - pcauchy((-pi - param[1])/param[2])
  g <- (1/param[2]) * dcauchy((x-param[1])/param[2]) * (1/c)
  return(g)
}

############################### Normal Truncated ###############################
#' @title Probability density function
#' @name dTruncNorm
#'
#' @description The truncated normal is a normal distribution bounded between
#'              low and high (the pdf is zero outside these bounds and
#'              renormalized). Samples from this distribution are differentiable
#'              with respect to loc, scale as well as the bounds,low and high,
#'              i.e., this implementation is fully reparameterizeable
#'
#' @param psi  Equivalent to \code{sample=psi}, is the array describing the support of the
#'             phase at which the pdf.
#' @param mu  loc: Floating point tensor; the means of the distribution(s).
#' @param sigma scale: loating point tensor; the stddevs of the distribution(s).
#'              Must contain only positive values.
#'
#' @details In probability and statistics, the truncated normal distribution is
#'          the probability distribution derived from that of a normally
#'          distributed random variable by bounding the random variable from
#'          either below or above (or both). The truncated normal distribution
#'          has wide applications in statistics and econometrics.
#'
#' @importFrom stats dnorm
#'
#' @importFrom stats pnorm
#'
#' @return A graphic
#'
#' @author Araújo, J. E.
#'
#' @examples
#' sample = seq(-pi,pi, 0.001)
#' param <- c(mu = 0.5, sigma = 0.5)
#' data<-data.frame(x=sample,y=dTruncNorm(sample, param))
#' ggplot(data,aes(x,y)) +
#'  geom_line() +
#'  xlim(-pi, pi) +
#'  # coord_polar() +
#'  xlab("sample") +
#'  ylab("density") +
#'  theme_pander()
#'
##Empirical model
##Distribution Normal Truncated
#' @export
dTruncNorm <- function(x, param){
  # Transform the value of x to the interval (-π, π), making the function periodic.
  x <- (x - 2 * pi * floor((x + pi) / (2 * pi)))
  c <- pnorm((pi - param[1])/param[2]) - pnorm((-pi - param[1])/param[2])
  f <- (1/param[2]) * dnorm((x-param[1])/param[2]) * (1/c)
  return(f)
}

#################### Maximum Likelihood Estimation (MLE) ######################

################################## dFuncLeeEq4 ################################

#' @title Maximization with the maxLik package
#'
#' @name estim.LeeEq5
#'
#' @description If it is necessary to maximize the parameters of the dFuncLeeEq4,
#'              density, we indicate the maxLik package.
#'
#' maxLik::maxLik(LogLikelihoodLeeEq5, start=c(r, theta,L))
#'
#' @details The main function of the maxLik package is to perform the Maximum
#'          Likelihood estimate.
#'
#' @examples ## Estimate the parameters of Lee et al. (1994), Eq. 18 Density.
#' sample <- c( 0.21850138, 0.74192327, 1.42014964, 0.03348661, 0.65550808,
#' 0.46397016, 0.03027485)
#' param <- c(r = 0.7, theta = 0, L = 5)
#' mle.Lee <- estim.LeeEq5(sample, param)
#' summary(mle.Lee)
#' @export
estim.LeeEq5 <- function(sample, param){
  LogLikelihoodLeeEq5 <- function(param){
    return(
      sum(log(dFuncLeeEq5(sample, param))))
  }
  return(maxLik(LogLikelihoodLeeEq5, start = param, fixed = "L"))
}

################################ dFuncGierullEq7 ###############################

#' @title Maximization with the maxLik package
#'
#' @name estim.GierullEq7
#'
#' @description If it is necessary to maximize the parameters of the dFuncGierullEq7,
#'              density, we indicate the maxLik package.
#'
#' @details The main function of the maxLik package is to perform the Maximum
#'          Likelihood estimate.
#'
#' maxLik::maxLik(LogLikelihoodGierullEq7, start=c(r, theta, L))
#'
#' @examples ## Estimate the parameters of Gierull (2020), Eq. 13. Density.
#' sample <- c( 0.21850138, 0.74192327, 1.42014964, 0.03348661, 0.65550808,
#' 0.46397016, 0.03027485)
#' param <- c(r = 0.7, theta = 0, L = 5)
#' mle.Gierull <- estim.GierullEq7(sample, param)
#' summary(mle.Gierull)
#'
#' @export
estim.GierullEq7 <- function(sample, param){
  LogLikelihoodGierullEq7 <- function(param){
    return(
      sum(log(dFuncGierullEq7(sample, param))))
  }
  return(maxLik(LogLikelihoodGierullEq7, start = param, fixed = "L"))
}

################################ dTruncCauchy ##################################

#' @title Maximization with the maxLik package
#'
#' @name estim.TruncCauchy
#'
#' @description If it is necessary to maximize the parameters of the dTruncCauchy,
#'              density, we indicate the maxLik package.
#'
#' @details The main function of the maxLik package is to perform the Maximum
#'          Likelihood estimate.
#'
#' maxLik::maxLik(LogLikelihoodTCauchy, start=c(mu, sigma))
#'
#' @examples ## Estimate the parameters of Truncated Cauchy Density
#' sample <- c( 0.21850138, 0.74192327, 1.42014964, 0.03348661, 0.65550808,
#' 0.46397016, 0.03027485)
#' param <- c(mu = 0.5, sigma = 1)
#' mle.TruncCauchy <- estim.TruncCauchy(sample, param)
#' summary(mle.TruncCauchy)
#'
#' @export
estim.TruncCauchy <- function(sample, param){
  LogLikelihoodTCauchy <- function(param){
    return(
      sum(log(dTruncCauchy(sample, param))))
  }
  return(maxLik(LogLikelihoodTCauchy, start = param))
}

################################ dTruncNorm ####################################

#' @title Maximization with the maxLik package
#'
#' @name estim.TruncNorm
#'
#' @description If it is necessary to maximize the parameters of the dTruncNorm,
#'              density, we indicate the maxLik package.
#'
#' @details The main function of the maxLik package is to perform the Maximum
#'          Likelihood estimate.
#'
#' maxLik::maxLik(LogLikelihoodTNorm, start=c(mu, sigma))
#'
#' @examples ## Estimate the parameters of Truncated Normal Density
#' sample <- c( 0.21850138, 0.74192327, 1.42014964, 0.03348661, 0.65550808,
#' 0.46397016, 0.03027485)
#' param <- c(mu = 0.5, sigma = 1)
#' mle.TruncNorm <- estim.TruncNorm(sample, param)
#' summary(mle.TruncNorm)
#' @export
estim.TruncNorm <- function(sample, param){
  LogLikelihoodTNorm <- function(param){
    return(
      sum(log(dTruncNorm(sample, param))))
  }
  return(maxLik(LogLikelihoodTNorm, start = param))
}
