dTruncCauchy <- function(x, param){
  # Transform the value of x to the interval (-π, π), making the function periodic.
  x <- (x - 2 * pi * floor((x + pi) / (2 * pi))) #is the mod
  c <- pcauchy((pi - param[1])/param[2]) - pcauchy((-pi - param[1])/param[2])
  g <- (1/param[2]) * dcauchy((x-param[1])/param[2]) * (1/c)
  return(g)
}