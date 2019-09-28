require(bayesImageS)

set.seed(1234)
mask <- matrix(1,512,512)
neigh <- getNeighbors(mask, c(2,2,0,0))
blocks <- getBlocks(mask, 2)
k <- 5
beta <- log(1+sqrt(k))
res.sw <- swNoData(beta, k, neigh, blocks, niter=200)
z <- matrix(max.col(res.sw$z)[1:nrow(neigh)], nrow=nrow(mask))
image(z, xaxt = 'n', yaxt='n', col=rainbow(k), asp=1)

