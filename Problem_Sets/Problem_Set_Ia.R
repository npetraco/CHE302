#library(rmarkdown)
library(che302r)

# 1.
lam <- 3e-3
pN  <- h/lam
pN

vN <- pN/mN
vN


# 2.
lam <- 350e-9
pph <- h/lam
pph

mH2 <- 2*mP + 2*me
mH2

pH2 <- pph
vH2 <- pH2/mH2
vH2


# 3.
dv <- 1000
dp <- mP*dv
dp

dx <-  hb/(2*dp)
dx


# 4.
lam <- 525e-9
Eph <- h * cl/lam
Eph # Energy per photon

# 500 W laser on for 10 s
ph <- 500 * (1/Eph) * 10 # number of photons spit out in 10s
ph
ph/N.A                   # mols of photons


# 5.
# a.
laser.pow  <- 5e6
pulse.time <- 2e-8
laser.pow * pulse.time

# b.
lam <- 1064e-9
Eph <- h * cl/lam
Eph # Energy per photon

# alternative for a.
# 5e6 W laser on for 2e-8 s
ph <- 5e6 * (1/Eph) * 2e-8 # number of photons spit out in 2e-8 s pulse
ph
Eph * ph                   # Energy in the 2e-8 s pulse

# c. Energy per pulse * 10 pulses
laser.pow * pulse.time * 10


# 6.
Phi    <- 2.14              # Work function in eV
lambda <- 300e-9            # Photon wavelength
#lambda <- 600e-9           # The other Photon's wavelength


PhiJ <- Phi * 1.602177e-19 # Convert Phi to J
Eph  <- h * cl/lambda      
Eph                        # Energy of the photon

KEe <- Eph - PhiJ
KEe                        # KE of the kicked out e-


# 7.
Phi    <- 3.84              # Work function in eV
lambda <- 170e-12           # Photon wavelength
#lambda <- 170e-6           # Photon wavelength

PhiJ <- Phi * 1.602177e-19 # Convert Phi to J
Eph  <- h * cl/lambda

# KE of the kicked out e-
KEe <- Eph - PhiJ
KEe

# Speed of the kicked out e-
ve <- sqrt(2*(Eph - PhiJ)/me)
ve


# 8.
# 1e- accelerated by 134V
# 1V = 1J/C
Eelec <- 134 * ec # ec = charge of an electron in Coulombs
Eelec

# From E=hc/lambda and p=h/lambda:
pelec <- Eelec/cl
pelec


# 9.
dx <- 150e-12
dp <- hb/(2*dx)
dp

dv <- dp/me
dv


# 10.
dt <- 1e-9
dE <- hb/(2*dt)
dE


# 11.
k100   <- planck.distribution(x.min = 0.1e-9, x.max = 200000e-9, Temp = 100, typ = "wavelength", plotQ = T)
k1000  <- planck.distribution(x.min = 0.1e-9, x.max = 200000e-9, Temp = 1000, typ = "wavelength", plotQ = T)
k10000 <- planck.distribution(x.min = 0.1e-9, x.max = 200000e-9, Temp = 10000, typ = "wavelength", plotQ = T)


# 13.
# a. True
# b. False
# c. False
# d. False
