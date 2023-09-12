#-------------------------------------------------------------------------------------
#This is a function to plot the Plank distribution as a function of wavenumber.
#
#ARGUEMENTS:
#nut.min: The smallest wavenumber to use. Use something at least little larger than 0.  
#nut.max: The largest wavenumber to use.
#Temp: Temperature in Kelvin
#--------------------------------------------------------------------------------------
planck.distribution.nutilde<-function(nut.min, nut.max, Temp) {
  
  #Make a nu-tilde (nut) axis. This is the x-axis.
  nut <- seq(from=nut.min, to=nut.max, length.out=2500)
  
  #Planck's dist as a function of nu-tilde (nut). This is the y-axis.
  rho <- 2*h*(cl^2)*(nut^3) * (1/(exp((h*cl*nut)/(kB*Temp))-1))
  
  #Make the plot.
  plot(nut, rho, typ="l", xlab="nu-tilde (m^-1)", ylab="Intensity", main="Planck's dist")
  
}


#Now USE the function planck.distribution.nutilde below:
