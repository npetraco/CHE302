library(rgl)
library(orthopolynom)
library(parallel)
library(doMC)

#Convert spherical coordinates to cartesian coordinates
spherical.to.cartesian<-function(rl,tha,pha) {
  x<- rl * sin(tha) * cos(pha)
  y<- rl * sin(tha) * sin(pha)
  z<- rl * cos(tha)
  
  return(c(x,y,z))
}


#Sample the wave function's density and plot corresponding r,theta, phi (x,y,z) values
#Recipie adapted from D. Cromer, J Chem Educ. 45(10) 626-631 1968
sample.density2<-function(nqn, lqn, mqn, spherical.grid, num.samples) {
  
  r<-spherical.grid[,1]
  theta<-spherical.grid[,2]
  phi<-spherical.grid[,3]
  
  #den.vals<-Re(psi.func(r,theta,phi)^2) #The Imaginary parts should be 0, but still need to extract Re
  #den.vals<-orbital_density(nqn, lqn, mqn, spherical.grid) #Use the boost functions
  den.vals <- r^2 * sin(theta) * Re((Radial(n, l, r) * LegendreP(l, m, theta) * PhiF(m, phi))^2)
  den.max<-max(den.vals)                  #Approx max of the density
  pratio.vals<-den.vals/den.max           #Ratio for monte carlo selection
  
  count<-0
  num.iter.keep<-num.samples #For now and speed sample the number of samples. This is a big change from Cromers paper. Much faster but may be wrong. 
  #Prune the grid around the most probable values:
  keep.grid<-NULL
  while(count<num.samples) {
    q.iter<-runif(length(pratio.vals))
    idxs<-which(pratio.vals>=q.iter)  
    if(length(idxs>0)) {
      idx.pick<-sample(idxs,num.iter.keep)
      keep.grid<-rbind(keep.grid,spherical.grid[idx.pick,])
      count<-(count+length(idx.pick))
      #print(count)
    }
  }
  orbital<-t(sapply(1:num.samples,function(x){spherical.to.cartesian(keep.grid[x,1],keep.grid[x,2],keep.grid[x,3])}))
  
  return(orbital)
  
}


#
#Radial wavefunction using the orthopoly package
#

#Generates the symbolic expression for the desired Laguerre polynomial
laguerreLgen<-function(nn,ll){
  
  poly.nom <- glaguerre.polynomials(nn, ll, normalized=F)
  poly.nom <- list(poly.nom[[length(poly.nom)]]) #to gen polynomial values orthopolynom func expexts a list
  
  return(poly.nom)
  
}

#Generate a numerical representation of the desired radial wavefunction
Radial<-function(nn,ll,rr){
  
  lag.poly <- laguerreLgen(nn-ll-1, 2*ll+1)
  #print(lag.poly)
  
  Radial.vals <- exp(-rr/nn) * rr^ll * polynomial.values(lag.poly,rr)[[1]]
  
  return(Radial.vals)
  
}

#
#Angular wavefunction using the Numerical Recipies routine
#

#There is no associated Legendre polynomial or shpherical harmonic in orthopolynom, so we have to do it ourselves

#Generates the symbolic expression for a Legendre polynomial
legendrePgen<-function(ll){
  
  poly.nom <- legendre.polynomials(ll, normalized=F)
  poly.nom <- list(poly.nom[[length(poly.nom)]]) #to gen polynomial values orthopolynom func expects a list
  
  return(poly.nom)
  
}

#Generates the symbolic expression for |m|^th derivative Legendre polynomial. Needed for Associated Legendre polynomials
mderivlegendrePgen<-function(ll,mm){
  
  polyd <- legendrePgen(l)
  if(abs(mm)>0){
    for(i in 1:abs(mm)){
      polyd <- polynomial.derivatives(polyd)
    }
  }
  
  return(polyd)
}

#A pure function for the associated Legendre polynomial
LegendrePfunc<-function(ll,mm){
  
  deriv.func <- polynomial.functions(mderivlegendrePgen(ll,mm))[[1]]
  
  Plmx <- function(xx){
    (-1)^abs(mm) * (1-xx^2)^(abs(mm)/2) * deriv.func(xx)
  }
  
  return(Plmx)
  
}

#Compute a numerical representation for the associated Legendre polynomial
LegendreP<-function(ll,mm,thetaa){
  
  xxx<-cos(thetaa)
  vals <- LegendrePfunc(ll,mm)(xxx)
  
  if(mm<0){
    vals <- (-1)^mm * factorial(ll-mm)/factorial(ll+mm) * vals
  }
  
  return(vals)
  
}

#Compute the (numerical) Phi wave function for the equitorial angle
PhiF<-function(mm,phii){
  
  return(exp((mm*complex(real = 0, imaginary = 1))*phii))
  
}



#---------------------------------------------------
# MCMC subroutines
#---------------------------------------------------
# Generate spline versions of the squared orbital component functions:
# This will save us from having to re-evaluate using the slower polynomial functions during MCMC
splined.sqorb.functs <- function(n, l, m, r.max = 20, num.knots = 1000){
  
  r     <- seq(from=1e-6, to=r.max, length.out = num.knots)
  theta <- seq(from=0,    to=pi,    length.out = num.knots)
  phi   <- seq(from=0,    to=2*pi,  length.out = num.knots)
  
  R.sq.vals <- sapply(r, function(r){Radial(n, l, r)^2})
  R.sq.func <- splinefun(r, R.sq.vals)
  #plot(r, r^2*R.den.func(r),typ="l")      # Check: is tail of Radial part long enough?
  #plot(r, log(r^2*R.den.func(r)),typ="l") # Check: is number of Radial nodes correct?
  
  Th.sq.vals <- sapply(theta, function(theta){sin(theta) * LegendreP(l, m, theta)^2})
  Th.sq.func <- splinefun(theta, Th.sq.vals)
  #plot(theta, sin(theta)*Th.den.func(theta), typ="l")
  #plot(theta, log(sin(theta)*Th.den.func(theta)), typ="l")
  
  Ph.sq.vals <- sapply(phi, function(phi){Re(PhiF(m, phi))^2})
  Ph.sq.func <- splinefun(phi, Ph.sq.vals)
  #plot(phi, Ph.den.func(phi), typ="l")
  #plot(phi, log(Ph.den.func(phi)), typ="l")
  
  spline.orb.func.info <- list(
    r,
    theta,
    phi,
    R.sq.func,
    Th.sq.func,
    Ph.sq.func
  )
  
  names(spline.orb.func.info) <- c(
    "r.axis",
    "theta.axis",
    "phi.axis",
    "R.sq.func",
    "Theta.sq.func",
    "Phi.sq.func"
  )

  return(spline.orb.func.info)
  
}

# Next step proposal for Metropolis algorithm, generalized to any number of parameters.
# CAUTION: uses Gaussian proposal
proposal <- function(theta.curr, proposal.wid){sapply(1:length(theta.curr), function(xx){rnorm(n = 1, mean = theta.curr[xx], sd = proposal.wid[xx])})}

# Orbital Ansatz using spline functions
# Log Density ansatz. Pass in splined function info to hide from user. ****Looks yucky though!
ansatz <- function(a.theta, R.den.func, Th.den.func, Ph.den.func, rax.ext){
  
  if((a.theta[1] > 0) & (a.theta[1] <= rax.ext)){
    val1 <- log( a.theta[1]^2 * R.den.func(a.theta[1])  )  # Replaces: Radial(n, l, a.theta[1])^2 
  } else {
    val1 <- log(0)
  }
  
  if( (a.theta[2] >= 0) & (a.theta[2] <= pi) ) {
    val2 <- log( sin(a.theta[2]) * Th.den.func(a.theta[2]) ) # Replaces: (LegendreP(l, m, a.theta[2]))^2
  } else {
    val2 <- log(0)
  }
  
  if( (a.theta[3] >= 0) & (a.theta[3] <= 2*pi)){
    val3 <- log( Ph.den.func(a.theta[3]) )                   # Replaces: Re(PhiF(m, a.theta[3]))^2
  } else {
    val3 <- log(0)
  }
  
  val <- val1 + val2 + val3
  
  return(val)
}

# Log Metropolis ratio. Pass in splined function info to hide code...... Looks yucky
r.metrop <- function(theta.curr, theta.prop, rsf, thsf, phsf, rext){
  val <- ansatz(theta.prop, rsf, thsf, phsf, rext) - ansatz(theta.curr, rsf, thsf, phsf, rext)
  if(is.nan(val)){
    val <- (-Inf)
  }
  return(val)
} 

# MCMC Metropolis routine:
# Should execute for any number of dimensions for a parameter vector.
# Requires that an ansatz (a not necessarily normalized pdf) be defined
sampler.gen <- function(num.iter=100, theta.init=c(0.5,0.5), proposal.width=c(0.5,0.5), spline.func.info){
  
  theta.current     <- theta.init
  ansatz.sample     <- array(NA,c(num.iter, length(theta.init))) # to theta sample from ansatz
  ansatz.sample[1,] <- theta.current
  accept.ps         <- array(NA,num.iter-1) # capture all acceptance probs just for checks
  accept.indcs      <- array(NA,num.iter-1) # capture all acceptance indicators just for checks
  
  # To pass into r.metrop and then to ansatz routine: YUK!!!!!!
  rsqf  <- spline.func.info$R.sq.func
  thsqf <- spline.func.info$Theta.sq.func
  phsqf <- spline.func.info$Phi.sq.func
  rmx   <- max(spline.func.info$r.axis)
  
  for(i in 2:num.iter){
    
    #print(i)
    
    # suggest new position
    theta.proposal <- proposal(theta.curr = theta.current, proposal.wid = proposal.width)
    
    # Accept proposal. Note log(ratio) and also pass in splined density parts. Can we do this better?? 
    accept.ratio <- r.metrop(
      theta.curr = theta.current, 
      theta.prop = theta.proposal, rsqf, thsqf, phsqf, rmx)
    
    acceptQ           <- log(runif(1)) < accept.ratio
    accept.ps[i-1]    <- min(1,exp(accept.ratio))     # Just for checking
    accept.indcs[i-1] <- acceptQ                      # Just for checking
    
    if(acceptQ) {
      # Update position
      theta.current <- theta.proposal
    }
    
    ansatz.sample[i,] <- theta.current
    
  } 
  
  return(list(ansatz.sample, accept.ps, accept.indcs))
  
}


#Sample the wave function's density and plot corresponding r,theta, phi (x,y,z) values
#Now uses MCMC Metropolis algorithm.

# Wrapper to run MCMC for the orbital densities:
sample.density3a<-function(nqn, lqn, mqn, orb.func.info, num.samples=5000, num.thin=1, num.burnin=0) {
  
  mpi  <- sampler.gen(
    theta.init       = runif(3, min = 0, max = 1), 
    proposal.width   = c(1,1,1), 
    num.iter         = num.samples, 
    spline.func.info = orb.func.info)
  
  mp   <- mpi[[1]] # Chain
  ap   <- mpi[[2]] # Acceptance probs
  ai   <- mpi[[3]] # Acceptance indicators
  
  # Remove burn-in and Thin:
  keep.idxs <- seq(num.burnin+1,nrow(mp),num.thin)
  mp2  <- mp[keep.idxs, ]
  
  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(mp2),function(xx){spherical.to.cartesian(mp2[xx,1],mp2[xx,2],mp2[xx,3])}))
  
  colnames(samp.coords) <- c("x","y","z")

  print(paste("The final sample size is:", nrow(samp.coords)))
  
  return(samp.coords)  
  
}

# Wrapper to run MCMC for the orbital densities:
sample.density3b<-function(nqn, lqn, mqn, orb.func.info, num.samples=5000, num.thin=1, num.burnin=0, jump.widths=c(1,1,1), printQ=FALSE) {
  
  mpi  <- sampler.gen(
    theta.init       = runif(3, min = 0, max = 1), 
    proposal.width   = jump.widths, 
    num.iter         = num.samples, 
    spline.func.info = orb.func.info)
  
  mp   <- mpi[[1]] # Chain
  ap   <- mpi[[2]] # Acceptance probs
  ai   <- mpi[[3]] # Acceptance indicators
  
  if(printQ==TRUE) {
    perc.1s    <- length(which(ap==1))/length(ap) * 100
    accpt.rate <- sum(ai)/length(ai) * 100 # Acceptance rate
    print(paste0("Higher move percentage: ", round(perc.1s,2), "%"))
    print(paste0("Acceptance rate:        ", round(accpt.rate,2), "%"))
  }
  
  # Remove burn-in and Thin:
  keep.idxs <- seq(num.burnin+1,nrow(mp),num.thin)
  mp2  <- mp[keep.idxs, ]
  
  # if(printQ==TRUE) {
  #   par(mfrow=c(3,1))
  #   acf(mp2[,1])
  #   acf(mp2[,2])
  #   acf(mp2[,3])
  # }
  
  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(mp2),function(xx){spherical.to.cartesian(mp2[xx,1],mp2[xx,2],mp2[xx,3])}))
  samp.coords <- cbind(samp.coords, mp2)
  
  colnames(samp.coords) <- c("x","y","z","r","theta","phi")
  
  print(paste("The final sample size is:", nrow(samp.coords)))
  
  return(samp.coords)  
  
}


# VERY stripped down wrapper (to a wrapper: sample.density3b) to run MCMC twice. First run is a very correlated cheat to get the inner radial nodes. 
# Second run is a regular run to get the rest of the orbital density
sample.orbital.density<-function(nqn, lqn, mqn, orb.func.info, printQ=TRUE) {
  
  print("====== First (cheat) pass to try and get inner nodes (if any) ======")
  density.samples.loc <- sample.density3b(nqn, lqn, mqn, orb.func.info, 
                                          num.samples=25000,               # Preset by experience
                                          num.thin = 40,                   # Preset by experience
                                          jump.widths = c(0.25,2.38,2.38), # Preset by experience
                                          num.burnin = 2,                  # Preset by experience
                                          printQ)
  x.loc <- density.samples.loc[,1]
  y.loc <- density.samples.loc[,2]
  z.loc <- density.samples.loc[,3]
  print("Chain 1 done.")
  for(i in 2:4){
    density.samples.loc <- sample.density3b(nqn, lqn, mqn, orb.func.info, 
                                            num.samples=25000,               # Preset by experience
                                            num.thin = 40,                   # Preset by experience
                                            jump.widths = c(0.25,2.38,2.38), # Preset by experience
                                            num.burnin = 2,                  # Preset by experience
                                            printQ)
    x.loc <- c(x.loc, density.samples.loc[,1])
    y.loc <- c(y.loc, density.samples.loc[,2])
    z.loc <- c(z.loc, density.samples.loc[,3])
    print(paste("Chain", i, "done."))
    
  }
  
  print("====== Second pass to get outer parts ======")
  density.samples2.loc <- sample.density3b(nqn, lqn, mqn, orb.func.info, 
                                           num.samples=125000,              # Preset by experience
                                           num.thin = 150,                  # Preset by experience
                                           jump.widths = c(2.38,2.38,2.38), # Preset by experience
                                           num.burnin = 100,                # Preset by experience
                                           printQ)
  
  x.loc <- c(x.loc, density.samples2.loc[,1])
  y.loc <- c(y.loc, density.samples2.loc[,2])
  z.loc <- c(z.loc, density.samples2.loc[,3])
  print("Chain 1 done.")
  
  for(i in 2:4) {

    density.samples2.loc <- sample.density3b(nqn, lqn, mqn, orb.func.info, 
                                             num.samples=125000,              # Preset by experience
                                             num.thin = 150,                  # Preset by experience
                                             jump.widths = c(2.38,2.38,2.38), # Preset by experience
                                             num.burnin = 100,                # Preset by experience
                                             printQ)
    
    x.loc <- c(x.loc, density.samples2.loc[,1])
    y.loc <- c(y.loc, density.samples2.loc[,2])
    z.loc <- c(z.loc, density.samples2.loc[,3])
    print(paste("Chain", i, "done."))
        
  }
  
  samp.coords <- cbind(x.loc, y.loc, z.loc)
  
  colnames(samp.coords) <- c("x","y","z")
  
  print(paste("The total final sample size is:", nrow(samp.coords)))
  
  return(samp.coords)  
  
}


# Simple parallel version of sample.orbital.density(). Run the chains in separate processes
sample.orbital.density.parallel<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE, num.processes=1) {
  
  
  registerDoMC(num.processes)
  print(paste("Using",getDoParWorkers(), "processes."))
  
  print("====== First (cheat) pass to try and get inner nodes (if any) ======")
  density.samples.loc <- foreach(i=1:4, .combine="rbind") %dopar% {  
    density.samples.loc <- sample.density3b(nqn, lqn, mqn, orb.func.info, 
                                            num.samples=25000,               # Preset by experience
                                            num.thin = 40,                   # Preset by experience
                                            jump.widths = c(0.25,2.38,2.38), # Preset by experience
                                            num.burnin = 2,                  # Preset by experience
                                            printQ)
  }
  x.loc <- density.samples.loc[,1]
  y.loc <- density.samples.loc[,2]
  z.loc <- density.samples.loc[,3]

  print("====== Second pass to get outer parts ======")
  density.samples2.loc <- foreach(i=1:4, .combine="rbind") %dopar% {
    sample.density3b(nqn, lqn, mqn, orb.func.info, 
                     num.samples=125000,              # Preset by experience
                     num.thin = 150,                  # Preset by experience
                     jump.widths = c(2.38,2.38,2.38), # Preset by experience
                     num.burnin = 100,                # Preset by experience
                     printQ)
  }
  
  x.loc <- c(x.loc, density.samples2.loc[,1])
  y.loc <- c(y.loc, density.samples2.loc[,2])
  z.loc <- c(z.loc, density.samples2.loc[,3])

  samp.coords <- cbind(x.loc, y.loc, z.loc)
  
  colnames(samp.coords) <- c("x","y","z")
  
  print(paste("The total final sample size is:", nrow(samp.coords)))
  
  return(samp.coords)  
  
}