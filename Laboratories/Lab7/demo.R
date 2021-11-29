library(colorscience)
library(che302r)

# SPD -> XYZ (Y scale??)-> xyz -> xy -> color diagram filled data

# Pigment reflectance spectrum [S(lambda)]
raw <- read.csv("/Users/karen2/latex/class/chem302/R/color/reflectance_data/Database Reflectance Spectra Pigments + acrylic binder on cardboard/vermilion nat.txt", sep = "\t")

lam.srd <- raw[,1]
int1    <- raw[,2]/100 # ?? Scale??
plot(lam.srd, int1, typ="l", xlab="lambda (nm)", ylab="Reflectance", main="Spectral Reflectance Distribution")

# Color matching function data:
data("ciexyz31")
lam.cie <- ciexyz31[,1] # 1nm resolution
xfd     <- ciexyz31[,2]
yfd     <- ciexyz31[,3]
zfd     <- ciexyz31[,4]
plot(lam.cie, xfd, ylab="x(lambda)", xlab="lambda (nm)", main="x-CMF", typ="l")
plot(lam.cie, yfd, ylab="y(lambda)", xlab="lambda (nm)", main="y-CMF", typ="l")
plot(lam.cie, zfd, ylab="z(lambda)", xlab="lambda (nm)", main="z-CMF", typ="l")

# Reference Illuminant data:
data("illuminantD65")
all.illd           <- illuminantD65
ill.lambda.min.idx <- which(all.illd[,1] == min(lam.cie))
ill.lambda.max.idx <- which(all.illd[,1] == max(lam.cie))
illd <- illuminantD65[ill.lambda.min.idx:ill.lambda.max.idx,2]/100
plot(lam.cie, illd, ylab="I(lambda)", xlab="lambda (nm)", main="D65 Standard Illuminant", typ="l")

# Re-sample all data to be at the worst (1nm) resolution
Ill <- splinefun(lam.cie, illd) # Reflectance spectrum in full resolution

Sfr1 <- splinefun(lam.srd, int1) # Reflectance spectrum in full resolution

x.cmf <- splinefun(lam.cie, xfd) # x CMF, 1nm resolution
y.cmf <- splinefun(lam.cie, yfd) # x CMF, 1nm resolution
z.cmf <- splinefun(lam.cie, zfd) # x CMF, 1nm resolution

plot(lam.cie, Sfr1(lam.cie),  ylab="S(lambda)", xlab="lambda (nm)", main="Spectral Reflectance Distribution", typ="l")
plot(lam.cie, Ill(lam.cie),   ylab="I(lambda)", xlab="lambda (nm)", main="D65 Standard Illuminant", typ="l")
plot(lam.cie, x.cmf(lam.cie), ylab="x(lambda)", xlab="lambda (nm)", main="x-CMF", typ="l")
plot(lam.cie, y.cmf(lam.cie), ylab="y(lambda)", xlab="lambda (nm)", main="y-CMF", typ="l")
plot(lam.cie, z.cmf(lam.cie), ylab="z(lambda)", xlab="lambda (nm)", main="z-CMF", typ="l")


plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*x.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,0.4))
plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*y.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,0.4))
plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*z.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,0.4))

plot(lam.cie, Ill(lam.cie)*y.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,1.1))
plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*x.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,1.1))
plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*y.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,1.1))
plot(lam.cie, Sfr1(lam.cie)*Ill(lam.cie)*z.cmf(lam.cie), ylab="", xlab="lambda (nm)", main="", typ="l", ylim=c(0,1.1))

d.lambda <- 1 # resolution is 1nm
N        <- sum(Ill(lam.cie)*y.cmf(lam.cie)*d.lambda)
X        <- (1/N)*sum(Sfr1(lam.cie)*Ill(lam.cie)*x.cmf(lam.cie)*d.lambda)
Y        <- (1/N)*sum(Sfr1(lam.cie)*Ill(lam.cie)*y.cmf(lam.cie)*d.lambda)
Z        <- (1/N)*sum(Sfr1(lam.cie)*Ill(lam.cie)*z.cmf(lam.cie)*d.lambda)

XYZ <- c(X,Y,Z)
XYZ

x    <- X/sum(c(X,Y,Z))
y    <- Y/sum(c(X,Y,Z))
xyY <- c(x,y,Y*100)
xyY

# Convert to closest sRGB and plot swatch:
XYZ2sRGB(XYZ)
plot(1, pch=16, cex=12, col=XYZ2sRGB(XYZ)$hex)
