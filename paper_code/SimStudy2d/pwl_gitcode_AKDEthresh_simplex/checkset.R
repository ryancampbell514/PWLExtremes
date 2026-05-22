euclid.dist = function (x1, x2) 
{
  return(sqrt(sum((x1 - x2)^2)))
}


checkset.2d = function (x, y, w, r0w, k = 1) 
{
  x11 <- x[1]
  x12 <- x[2]
  y11 <- y[1]
  y12 <- y[2]
  if (is.null(dim(w))) {
    w <- cbind(w, 1 - w)
  }
  i1 <- which.min(apply(w, 1, euclid.dist, x2 = c(x11, y11)/(x11 + 
                                                               y11)))
  i2 <- which.min(apply(w, 1, euclid.dist, x2 = c(x12, y11)/(x12 + 
                                                               y11)))
  i3 <- which.min(apply(w, 1, euclid.dist, x2 = c(x11, y12)/(x11 + 
                                                               y12)))
  i4 <- which.min(apply(w, 1, euclid.dist, x2 = c(x12, y12)/(x12 + 
                                                               y12)))
  logic <- all((x11 + y11) > k * r0w[i1], (x12 + y11) > k * 
                 r0w[i2], (x11 + y12) > k * r0w[i3], (x12 + y12) > k * 
                 r0w[i4])
  return(logic)
}
