require(bayesImageS)
require(tiff)
require(reshape)
source("GammaSAR.R")
source("imagematrix.R")

### Example of simulation of a MRF
set.seed(1234)
### This is the image size
mask <- matrix(1,512,512)
### This is the neighborhood
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
### This is the number of classes and the parameter set to the transition
k <- 5
beta <- log(1+sqrt(k))
res.sw <- swNoData(beta, k, neigh, blocks, niter=200)
### This is the output matrix
z <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))
image(z, xaxt = 'n', yaxt='n', col=rainbow(k), asp=1)

### Example of a simulation of SAR images with the same road network

# Setup
set.seet(123456789)
mask <- matrix(1,1500,1500)
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
k <- 5
beta <- log(1+sqrt(k))

# Two images of classes
res.sw <- swNoData(beta, k, neigh, blocks, niter=200)
Classes1 <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))
plot(imagematrix(equalize(Classes1)))
imagematrixPNG(imagematrix(equalize(Classes1)), name = "../../Images/PNG/Classes1.png")

res.sw <- swNoData(beta, k, neigh, blocks, niter=200)
Classes2 <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))
plot(imagematrix(equalize(Classes2)))
imagematrixPNG(imagematrix(equalize(Classes2)), name = "../../Images/PNG/Classes2.png")

# Classes plus roads
roads <- readTIFF("../../Images/TIF/10528735_15.tif")

Classes1pRoads <- pmax(Classes1, roads*(k+1))
plot(imagematrix(equalize(Classes1pRoads)))
imagematrixPNG(imagematrix(equalize(Classes1pRoads)), name="../../Images/PNG/Class1Roads.png")

Classes2pRoads <- pmax(Classes2, roads*(k+1))
plot(imagematrix(equalize(Classes2pRoads)))
imagematrixPNG(imagematrix(equalize(Classes2pRoads)), name="../../Images/PNG/Class2Roads.png")

## Distributions and data for Classes1pRoads
L <- 3

# Very dark roads

SARroads <- rgammaSAR(sum(Classes1pRoads==6), L, 10)

# Increasingly brighter classes

SARclass1 <- rgammaSAR(sum(Classes1pRoads==1), L, 50)
SARclass2 <- rgammaSAR(sum(Classes1pRoads==2), L, 100)
SARclass3 <- rgammaSAR(sum(Classes1pRoads==3), L, 150)
SARclass4 <- rgammaSAR(sum(Classes1pRoads==4), L, 200)
SARclass5 <- rgammaSAR(sum(Classes1pRoads==5), L, 500)

# The observed SAR image

SAR1 <- Classes1pRoads

SAR1[Classes1pRoads==1] <- SARclass1
SAR1[Classes1pRoads==2] <- SARclass2
SAR1[Classes1pRoads==3] <- SARclass3
SAR1[Classes1pRoads==4] <- SARclass4
SAR1[Classes1pRoads==5] <- SARclass5
SAR1[Classes1pRoads==6] <- SARroads

plot(imagematrix(equalize(SAR1)))

## Distributions and data for Classes2pRoads
L <- 1

# Dark roads

SARroads <- rgammaSAR(sum(Classes2pRoads==6), L, 30)

# Increasingly brighter classes

SARclass1 <- rgammaSAR(sum(Classes2pRoads==1), L, 40)
SARclass2 <- rgammaSAR(sum(Classes2pRoads==2), L, 60)
SARclass3 <- rgammaSAR(sum(Classes2pRoads==3), L, 80)
SARclass4 <- rgammaSAR(sum(Classes2pRoads==4), L, 100)
SARclass5 <- rgammaSAR(sum(Classes2pRoads==5), L, 200)

# The observed SAR image

SAR2 <- Classes2pRoads

SAR2[Classes2pRoads==1] <- SARclass1
SAR2[Classes2pRoads==2] <- SARclass2
SAR2[Classes2pRoads==3] <- SARclass3
SAR2[Classes2pRoads==4] <- SARclass4
SAR2[Classes2pRoads==5] <- SARclass5
SAR2[Classes2pRoads==6] <- SARroads

plot(imagematrix(equalize(SAR2)))
imagematrixPNG(imagematrix(equalize(SAR2)), name = "../../Images/PNG/1look.png")

plot(imagematrix(equalize(SAR1)))
imagematrixPNG(imagematrix(equalize(SAR1)), name = "../../Images/PNG/3looks.png")

### Densities for the SAR1 image

x <- seq(.001, 400, by=.01)
d30 <- dexp(x, rate=1/30)
d40 <- dexp(x, rate=1/40)
d60 <- dexp(x, rate=1/60)
d80 <- dexp(x, rate=1/80)
d100 <- dexp(x, rate=1/100)
d200 <- dexp(x, rate=1/200)

DensitiesSAR1.melt <-
  melt(data=data.frame(x, d30, d40, d60, d80, d100, d200), measure.vars=-1)
names(DensitiesSAR1.melt) <- c("x", "Distr", "Density")

ggplot(DensitiesSAR1.melt, aes(x=x, y=Density, col=Distr)) + 
  geom_line(size=2, alpha=.7) +
  theme_few(base_size = 30)
ggsave(file="../../Plots/PDF/DensitiesSAR1.pdf")  

d10 <- dgammaSAR(x, 3, 10)
d50 <- dgammaSAR(x, 3, 50)
d100 <- dgammaSAR(x, 3, 100)
d150 <- dgammaSAR(x, 3, 150)
d200 <- dgammaSAR(x, 3, 200)
d500 <- dgammaSAR(x, 3, 500)
  
DensitiesSAR2.melt <-
  melt(data=data.frame(x, d10, d50, d100, d150, d200, d500), measure.vars=-1)
names(DensitiesSAR2.melt) <- c("x", "Distr", "Density")

ggplot(DensitiesSAR2.melt, aes(x=x, y=Density, col=Distr)) + 
  geom_line(size=2, alpha=.7) +
  theme_few(base_size = 30)
ggsave(file="../../Plots/PDF/DensitiesSAR2.pdf") 

