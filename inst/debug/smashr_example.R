load("smashr.RData")
out <- mixsqp(L,x0 = x0,control=list(numiter.em = 0,maxiter.sqp=16))

