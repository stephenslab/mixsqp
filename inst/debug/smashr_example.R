load("smashr.RData")
n    <- nrow(L)
w    <- rep(1/n,n)
out1 <- mixkwdual(L)
out2 <- mixsqp(L,x0 = x0,control=list(numiter.em = 0,maxiter.sqp = 10,eps = 0))

x <- x0
x <- x/sum(x)
u <- drop(L %*% x)
g <- drop(-t(L) %*% (w / u) + 1)
H <- crossprod((sqrt(w)/u) * L)
y <- solve(H,-g)
