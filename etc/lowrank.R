require(rjulia)
jDo("using LowRankApprox");
qr_julia = function(L, ptol = 1e-10){
  r2j(L,"L");
  jDo("Q,R,P = pqr(L); R = R[:,invperm(P)]");
  Q = j2r("Q"); R = j2r("R");
  return(list(Q = Q, R = R))
}
svd_julia = function(L, ptol = 1e-10){
  r2j(L,"L");
  jDo("u,d,v = psvd(L); R = R[:,invperm(P)]");
  u = j2r("u"); d = j2r("d"); v = j2r("v");
  return(list(u = u, d = d, v = v))
}
svd_R = function(L, rank = 100, q = 2){
  m = dim(L)[2]
  Y = L %*% matrix(stats::rnorm(m*m), m, m); rank = min(rank, dim(L))
  for(i in 1:q){
    Y <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
    Z <- crossprod(L,Y)
    Z <- qr.Q( qr(Z, complete = FALSE) , complete = FALSE )
    Y <- L %*% Z
  }#End for
  remove(Z)
  Q <- qr.Q( qr(Y, complete = FALSE) , complete = FALSE )
  remove(Y)
  B <- crossprod(Q,L)
  s = svd(B, nv = rank, nu = rank)
  s$Q = Q
  return(s)
}