#Pick an orbital:
n <- 1
l <- 0
m <- (0)

# Create squared orbital functions for sampling
orb.sq.funcs  <- splined.sqorb.functs(n, l, m, r.max = 320, num.knots = 1000)
rax           <- orb.sq.funcs$r.axis
R.sq          <- orb.sq.funcs$R.sq.func
plot(rax, rax^2 * R.sq(rax), typ="l")

#Sample the orbital density:
detectCores() # Number of processors available

#density.samples <- sample.orbital.density.ais(n, l, m, printQ = T, orb.func.info = orb.sq.funcs)
#density.samples <- sample.orbital.density.mcmc(n, l, m, printQ = T, orb.func.info = orb.sq.funcs)
density.samples <- sample.orbital.density.ais.parallel(n, l, m, printQ = T, orb.func.info = orb.sq.funcs, num.processes = 2)
#density.samples <- sample.orbital.density.mcmc.parallel(n, l, m, printQ = T, orb.func.info = orb.sq.funcs, num.processes = 4)

# Sample in cartesian coordinates
x               <- density.samples[,1]
y               <- density.samples[,2]
z               <- density.samples[,3]

# Sample in spherical coordinates
r               <- density.samples[,4]
theta           <- density.samples[,5]
phi             <- density.samples[,6]

# Check sample to make sure we got all the nodes:
rax           <- orb.sq.funcs$r.axis
R.sq          <- orb.sq.funcs$R.sq.func
plot(rax, rax^2 * R.sq(rax), typ="l", lwd=1, col="green")
hist(r, bre=200, probability = T)

thax          <- orb.sq.funcs$theta.axis
Th.sq         <- orb.sq.funcs$Theta.sq.func
plot(thax, sin(thax) * Th.sq(thax), typ="l", lwd=1, col="green")
hist(theta, bre=200, probability = T)

phax          <- orb.sq.funcs$phi.axis
Ph.sq         <- orb.sq.funcs$Phi.sq.func
plot(phax, Ph.sq(phax), typ="l", lwd=1, col="green")
hist(phi, bre=200, probability = T)


#Plot the sampled orbital density:
#2D 
plot(x,y)
plot(x,z)
plot(y,z)

#3D
plot3d(x,y,z,type="s",radius=1.5,xlab="x",ylab="y",zlab="z")

# 3D slices
plot3d(x,y,rep(0,length(z)),type="s",radius=0.5,xlab="x",ylab="y",zlab="z")
plot3d(x,rep(0,length(y)),z,type="s",radius=0.5,xlab="x",ylab="y",zlab="z")
plot3d(rep(0,length(x)), y, z,type="s",radius=0.5,xlab="x",ylab="y",zlab="z")
