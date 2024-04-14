dTruncNorm <- function(x, param){
  # Transform the value of x to the interval (-π, π), making the function periodic.
  x <- (x - 2 * pi * floor((x + pi) / (2 * pi)))
  c <- pnorm((pi - param[1])/param[2]) - pnorm((-pi - param[1])/param[2])
  f <- (1/param[2]) * dnorm((x-param[1])/param[2]) * (1/c)
  return(f)
}