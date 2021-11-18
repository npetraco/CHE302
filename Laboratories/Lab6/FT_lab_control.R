library(che302r)

# Load the interferogram:
dat <- read.csv("https://raw.githubusercontent.com/npetraco/CHE302/master/Laboratories/Lab6/interferogram.csv", header=T)
interferogram <- dat[,2]

# Plot the interferogram:
plot(1:length(interferogram), interferogram, typ="l")
plot(interferogram[2400:2600], typ="l") # Zoom

# Fourier transform the interfrogram:
spect<-fft.interferogram2(interferogram, plot.typ="%T")
spect<-fft.interferogram2(interferogram, plot.typ="Absorbance")

# What molecule could this be?
head(spect[[1]])
plot(spect[[2]])
plot(spect[[2]][,2])
peak.idxs <- which(spect[[2]][,2] > 1)
spect[[2]][peak.idxs,1]
