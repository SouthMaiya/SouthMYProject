require(bayesImageS)
require(tiff)
source("GammaSAR.R")

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
image(Classes1, xaxt = 'n', yaxt='n', col=rainbow(k), asp=1)

res.sw <- swNoData(beta, k, neigh, blocks, niter=200)
Classes2 <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))
image(Classes2, xaxt = 'n', yaxt='n', col=rainbow(k), asp=1)

# Classes plus roads
roads <- readTIFF("../../Images/TIF/10528735_15.tif")

Classes1pRoads <- pmax(Classes1, roads*(k+1))
plot(imagematrix(equalize(Classes1pRoads)))

Classes2pRoads <- pmax(Classes2, roads*(k+1))
image(Classes2pRoads, xaxt = 'n', yaxt='n', col=rainbow(k+1), asp=1)

## Distributions
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
