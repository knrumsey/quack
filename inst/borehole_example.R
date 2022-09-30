library(GBASS)
library(BASS)
library(quack)
library(tictoc)
library(laGP)
library(BART)
library(lhs)
N <- 4000
M <- 5000

X <- smartLHS(N, 8)
y <- apply(X, 1, ff_borehole)
XX <- smartLHS(M, 8)
yy <- apply(XX, 1, ff_borehole)

# LEAP-GP
tic()
mod1 <- leapGP_build(X, y, verbose=TRUE)
t_leap <- toc()

yhat1 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat1[i] <- leapGP_predict(mod1, xx)
}
t_leap <- toc()

# LA-GP
t_la <- list(tic=0, toc=0)
yhat2 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat2[i] <- laGP(matrix(xx, nrow=1), start=floor(max(6, 0.1*sqrt(N))), end=ceiling(sqrt(nrow(X))),
                   X=X, Z=matrix(y, ncol=1))$mean
}
t_la2 <- toc()

# SLAP-GP
tic()
mod3 <- leapGP_synch(mod1)
t_slap <- toc()

rho_vec <- c(0, 0.25, 0.5, 0.8, 0.9, 0.95, 0.99, 1)
t_slap2 <- list()
yhat3 <- matrix(NA, nrow=length(rho_vec), ncol=M)

for(r in seq_along(rho_vec)){
  rho <- rho_vec[r]
  print(rho)

  #yhat3 <- rep(NA, M)
  hubs <- mod3$hubs
  tic()
  for(i in 1:M){
    xx <- XX[i,]
    emu <- slapGP(xx, X, y, hubs=hubs, rho=rho)
    yhat3[r,i] <- emu$pred
    hubs <- emu$hubs
  }
  t_slap2[[r]] <- toc()
}


#SLAP-GP FROM SCRATCH
rho_vec <- c(0, 0.25, 0.5, 0.8, 0.9, 0.95, 0.99, 1)
t_scratch2 <- list()
yhat4 <- matrix(NA, nrow=length(rho_vec), ncol=M)

for(r in seq_along(rho_vec)){
  rho <- rho_vec[r]
  print(rho)

  #yhat3 <- rep(NA, M)
  #hubs <- mod3$hubs
  hubs <- list()
  tic()
  for(i in 1:M){
    xx <- XX[i,]
    emu <- slapGP(xx, X, y, hubs=hubs, rho=rho)
    yhat4[r,i] <- emu$pred
    hubs <- emu$hubs
  }
  t_scratch2[[r]] <- toc()
}

rmse1 <- sqrt(mean((yy-yhat1)^2))
rmse2 <- sqrt(mean((yy-yhat2)^2))
rmse3 <- apply(yhat3, 1, function(yht, yy) sqrt(mean((yy-yht)^2)), yy=yy)
rmse4 <- apply(yhat4, 1, function(yht, yy) sqrt(mean((yy-yht)^2)), yy=yy)


#Make plot
library(RColorBrewer)
bob <- brewer.pal(5, "Set2")
plot(rho_vec, rmse4, type='o', pch=16, col=bob[1], lwd=2,
     xlab='rho', ylab="RMSE", cex.lab=1.4, ylim=c(0, 15))
lines(rho_vec, rmse3, lwd=2, col=bob[2])
points(rho_vec, rmse3, pch=16, col=bob[2])
abline(h=rmse1, col=bob[3], lwd=2)
abline(h=rmse2, col=bob[4], lwd=2)


#save(X, XX, mod1, mod2, mod3, yhat1, yhat2, yhat3, yhat4, t_la, t_la2, t_leap, t_leap2, t_slap, t_slap2,  t_scratch2,
#     file="data/gp_sims_40000")
