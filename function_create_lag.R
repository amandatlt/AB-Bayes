shift_by <- function(x, n) {
  n <- round(n)
  direction <- ifelse(n<0, "lag", "lead")
  n <- sqrt(n^2)
  v <- rep(NA,times = n)
  if (direction == "lag") {
    x <- x[1:(length(x) - n)]
    x <- c(v, x)
  } else {
    x <- x[(n+1):length(x)]
    x <- c(x,v)
  }
  return(x)
}