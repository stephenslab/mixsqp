load("smashr.RData")

out <- mixsqp(L,x0 = c(1,1),
              control = list(numiter.em = 0,maxiter.sqp = 20))
