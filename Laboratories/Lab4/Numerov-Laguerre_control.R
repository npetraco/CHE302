#source(file.choose())

#Domain (x-axis):
dr    <- 0.1                            # Increment on r-axis
r.min <- 1e-15                          # Start of r-axis.
r.max <- 20                             # End of r-axis.
r     <- seq(from=r.min,to=r.max,by=dr) # r-axis

#Numerov's procedure to solve the (almost) Radial Schrodinger Equation:
n            <- 1      # n = 1,2,3,...
l            <- 0      # l = 0, ... , n-1 (e.g. s,p,d,f,g,.....)
Guess.Energy <- (-0.5) # Initial guess for energy of the state
max.iter     <- 50                  
state        <- (n-l-1)                  
Fr.info      <- numerov.procedure(r,dr,"radial",state,Guess.Energy,max.iter,0.00)

#Transform F(r) back into R(r), the proper radial wave function
psi.info     <- Fr.info
psi.info[,2] <- (r^-1)*Fr.info[,2]
plot(r, r^2 * psi.info[,2],typ="l")

#Approximately Normalize Psi:
Npsi.info<-approx.normalize(psi.info, include.jacobian = TRUE)

#Plot the normalized wave function:
plot(r, r^2*Npsi.info[,2], xlim=c(min(r),max(r)),ylab="N psi(x)",typ="l")
