#Numerov's recursive formula for the wave function (f is the wave function. np1/nm1 means n plus/minus one)
f.point.plus.one<-function(f.point,f.point.minus.one,increment,G.point.plus.one,G.point,G.point.minus.one) {
  return((2*f.point - f.point.minus.one + (5*G.point*f.point*increment^2)/6 + (G.point.minus.one*f.point.minus.one*increment^2)/12)/(1 - (G.point.plus.one*increment^2)/12))
}

#G function:
G<-function(disc.V, an.Energy) {
  return(2*(disc.V - an.Energy))
}


#Potential Energy function
V<-function(xax,potential.name, l=NULL) {
  
  if( !(potential.name %in% c("box", "harmonic", "anharmonic", "radial")) ){
    #print("Error!")
    #print("Your choices are: box, harmonic, anharmonic or laguerre.")
    stop("Bad choice in PE function. Your choices are: box, harmonic, anharmonic or radial.")
  }
  
  if(potential.name=="box"){
    return(rep(0,length(xax)))         #Particle in a box V
  }
  
  if(potential.name=="harmonic"){
    return((1/2) * xax^2)            #Harmonic oscillator V
  }
  
  if(potential.name=="anharmonic"){
    return(151.29 * (1-exp(-xax))^2) #Morse V, approx anharmonic oscillator 
    #return(7.61 * (1-exp(-(xax-74.1) * 0.0193))^2)
  }
  
  if(potential.name=="radial"){
    return(   0.5*((l*(l+1))/xax^2) -(1/xax))  #... used to send in l
  }
  
}


#Numerov's procedure:
numerov.procedure<-function(xaxis, dxaxis, PE.function.name, nodes.in.state, E.guess, num.iterations, delay.time) {
  
  Elow<-NULL  #Initialization
  Ehigh<-NULL
  Er<-E.guess #First iteration Enenrgy
  
  for(iter in 1:num.iterations) {
    
    if(PE.function.name=="radial"){
      Gr<-G(V(xaxis,PE.function.name,l), Er)
    } else{
      Gr<-G(V(xaxis,PE.function.name), Er)
    }
    
    
    psir<-numeric(length(xaxis))  #Initialize a psir for the iteration
    psir[1]<-0
    psir[2]<-(1*10^(-4))
    
    #Compute the psir for this iteration
    num.nodes<-0
    for(i in 3:length(xaxis)) {
      psir[i] <- f.point.plus.one(psir[i-1], psir[i-2], dxaxis, Gr[i], Gr[i-1], Gr[i-2])
      if(psir[i-1]*psir[i]<0) {
        num.nodes<-num.nodes+1
      }
      #print(psir[i])
    }
    
    #Check the energy and update it:
    if(num.nodes<=nodes.in.state) {
      Elow<-c(Elow,Er)
      if(is.null(Ehigh)) {
        print("E too low, but no Ehigh. Add 1 to E.")
        Er<-(Er+1)
      } else {
        print("E too low. Find avg. of max Elow and min Ehigh")
        print(paste("E bracket: [",max(Elow),",",min(Ehigh),"]"))
        Er<-(max(Elow)+min(Ehigh))/2
      }
      
    } else {
      print("E too high. Average with max Elow")
      Ehigh<-c(Ehigh,Er)
      print(paste("E bracket: [",max(Elow),",",min(Ehigh),"]"))
      Er<-(max(Elow) + min(Ehigh))/2
    }
    
    print(paste("Iteration: ",iter,"   E: ",Er,"   # of nodes: ",num.nodes,sep=""))
    plot(xaxis,psir)
    Sys.sleep(delay.time)
  }
  return(cbind(xaxis,psir))
}


#Appriximately normalize a 1D wave function
approx.normalize<-function(wf.info, include.jacobian=FALSE, plotQ=FALSE) {
  xx<-wf.info[,1]
  wfx<-wf.info[,2]
  
  #Fit a spline function to the wave function^2 values:
  if(include.jacobian == TRUE) {
    den.func<-splinefun(xx, xx^2 * wfx^2) # Include Jacobian chunk when normalizing radial wave functions
  } else {
    den.func<-splinefun(xx,wfx^2)
  }
  
  #Compute the integral needed in the normalization constant = Int Psi* Psi dx:
  intgpsp<-integrate(den.func,lower=min(xx), upper=max(xx))$value
  
  #Compute the normalization constant:
  N.const<-1/sqrt(intgpsp)
  if(plotQ==TRUE) {
    plot(xx,(den.func(xx))^2,ylab="|psi(x)|^2")
    print(paste("Approx. normalization const.:",N.const))
  }
  
  #Scale the wave function values with the normalization constant
  wf.info.normalized<-wf.info
  wf.info.normalized[,2]<-(N.const*wf.info.normalized[,2])
  
  return(wf.info.normalized)
}