dat <- read.csv("/Users/karen2/latex/class/chem302/Labratories/Lab 8/color_lab_data.csv",header = T)

attach(dat)
colnames(dat)

S <- unk.SRD2
I <- illuminant.SPD

d.lammbda <- 1 # resolution is 1nm
N    <- sum(y.cmf * I * d.lammbda)
X    <- (1/N) * sum(x.cmf * S * I * d.lammbda)
Y    <- (1/N) * sum(y.cmf * S * I * d.lammbda)
Z    <- (1/N) * sum(z.cmf * S * I * d.lammbda)

XYZ <- c(X,Y,Z)
XYZ

x    <- X/sum(c(X,Y,Z))
y    <- Y/sum(c(X,Y,Z))
xyY <- c(x,y,Y*100)
xyY

# Convert to closest sRGB and plot swatch:
XYZ2sRGB(XYZ)
plot(1, pch=16, cex=12, col=XYZ2sRGB(XYZ)$hex)

detach(dat)