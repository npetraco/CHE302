library(che302r)

# Load the interferogram:
interferogram <- read.csv("https://raw.githubusercontent.com/npetraco/CHE302/master/Laboratories/Lab6/interferogram4.csv", header=T)
interferogram <- as.matrix(interferogram)

# Plot the interferogram:
plot(1:length(interferogram), interferogram, typ="l")

# Fourier transform the interfrogram:
spect.info <- fft.interferogram2(interferogram, plot.typ="Absorbance")

# Examine the spectrum:
spect <- spect.info[,c(1,2)]
plot.spectrum(spect)
find.peaks(spect, spectrum.typ = "Absorbance", deriv.tol = 0.0001, plotQ = T, peak.tol = 0.0025)

spect2 <- spect.info[,c(1,3)]
plot.spectrum(spect2)
find.peaks(spect2, spectrum.typ = "%T", deriv.tol = 0.0001, plotQ = T, peak.tol = 97.5)
