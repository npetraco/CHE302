#Pick an orbital:
n <- 3
l <- (2)
m <- (-1)

# Create squared orbital functions for sampling
orb.sq.funcs  <- splined.sqorb.functs(n, l, m, r.max = 50, num.knots = 1000)
rax           <- orb.sq.funcs$r.axis
R.sq          <- orb.sq.funcs$R.sq.func
plot(rax, rax^2 * R.sq(rax), typ="l")
plot(rax, log(rax^2 * R.sq(rax)), typ="l")


#Sample the orbital density:
density.samples <- sample.orbital.density.parallel(n, l, m, orb.func.info = orb.sq.funcs, num.processes = 4)
x               <- density.samples[,1]
y               <- density.samples[,2]
z               <- density.samples[,3]

#Plot the sampled orbital density:
#2D 
plot(x,y)
plot(x,z)
plot(y,z)
#3D
plot3d(x,y,z,type="s",radius=0.2,xlab="x",ylab="y",zlab="z")

