#source("ft_spectroscopy.func.R")
library(che302r)

# Load the interferogram:
dat <- read.csv("https://raw.githubusercontent.com/npetraco/CHE302/master/Laboratories/Lab6/interferogram.csv", header=T)
dat <- as.matrix(dat)

# Plot the interferogram:
plot(1:length(dat[,1]), dat[,1], typ="l")

# Fourier transform the interfrogram:
spect<-fft.interferogram2(dat, plot.typ="%T")
spect<-fft.interferogram2(dat, plot.typ="Absorbance")

# What molecule could this be?