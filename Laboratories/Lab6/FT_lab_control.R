source("https://npetraco.github.io/CHE302/Laboratories/Lab6/ft_spectroscopy.func.R")

# Load the interferogram:
dat <- read.csv("https://npetraco.github.io/CHE302/Laboratories/Lab6/interferogram.csv", header=F)
dat <- as.matrix(dat)

# Plot the interferogram:
plot(1:length(dat[,1]), dat[,1], typ="l")

# Fourier transform the interfrogram:
spect<-fft.interferogram2(dat, plot.typ="%T")
spect<-fft.interferogram2(dat, plot.typ="Absorbance")

# What molecule could this be?