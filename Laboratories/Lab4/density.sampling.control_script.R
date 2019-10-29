#Pick an orbital:
n <- 1
l <- (0)
m <- (0)

# Create squared orbital functions for sampling:
# Make sure the tail of your radial wave function is long enough!
r.max         <- 20
orb.sq.funcs  <- splined.sqorb.functs(n, l, m, r.max = r.max, num.knots = 1000)

#Sample the orbital density:
density.samples <- sample.density3a(n, l, m, orb.func.info = orb.sq.funcs, num.samples=50000, num.thin = 10, num.burnin = 200)
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
