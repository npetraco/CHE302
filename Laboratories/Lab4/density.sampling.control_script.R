library(che302r)
library(parallel)
library(doSNOW)
library(rgl)

#Pick an orbital:
n <- 1
l <- 0
m <- (0)

# Create squared orbital functions for sampling
orb.sq.funcs  <- splined.sqorb.functs(n, l, m, r.max = 70, num.knots = 1000)
rax           <- orb.sq.funcs$r.axis
R.sq          <- orb.sq.funcs$R.sq.func
plot(rax, rax^2 * R.sq(rax), typ="l")

#Sample the orbital density:
detectCores() # Number of processors available

#density.samples <- sample.orbital.density.ais(n, l, m, printQ = T, orb.func.info = orb.sq.funcs)
#density.samples <- sample.orbital.density.mcmc(n, l, m, printQ = T, orb.func.info = orb.sq.funcs)
density.samples <- sample.orbital.density.ais.parallel(n, l, m, printQ = T, orb.func.info = orb.sq.funcs, num.processes = 4)
#density.samples <- sample.orbital.density.mcmc.parallel(n, l, m, printQ = T, orb.func.info = orb.sq.funcs, num.processes = 4)
#density.samples <- sample.orbital.density(n, l, m, printQ = T, orb.func.info = orb.sq.funcs, algorithm = "mcmc", type = "parallel", num.processes = 4)

# Sample in cartesian coordinates
x               <- density.samples[,1]
y               <- density.samples[,2]
z               <- density.samples[,3]

#Plot the sampled orbital density:
#2D
plot(x,y)
plot(x,z)
plot(y,z)

#3D
plot3d(x,y,z,type="s",radius=0.1,xlab="x",ylab="y",zlab="z")

# 3D slices
plot3d(x,y,rep(0,length(z)),type="s",radius=0.5,xlab="x",ylab="y",zlab="z")
plot3d(x,rep(0,length(y)),z,type="s",radius=0.1,xlab="x",ylab="y",zlab="z")
plot3d(rep(0,length(x)), y, z,type="s",radius=0.5,xlab="x",ylab="y",zlab="z")
