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
