library(GBASS)
library(BASS)
library(quack)
library(tictoc)
library(laGP)
library(BART)

N <- 5001
X <- lhs::randomLHS(N, 3)
y <- apply(X, 1, GBASS::ff5)

# Make predictions
M <- 1000
XX <- lhs::randomLHS(M, 3)
yy <- apply(XX, 1, ff5)

# Fit a BASS model
tic()
mod1 <- bass(X, y)
t1b <- toc()

# Fit a leapGP model
tic()
mod2 <- leapGP_build(X, y, verbose=TRUE)
t2b <- toc()

# Make a slapGP model
mod3 <- leapGP_synch(mod2)

# Fit a BART model
tic()
mod5 <- gbart(X, y)
t5b <- toc()



# Bass predictions
yhat1 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat1[i] <- mean(predict(mod1, xx))
}
t1 <- toc()

# leapGP predictions
yhat2 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat2[i] <- leapGP_predict(mod2, xx)
}
t2 <- toc()

# slapGP predictions
yhat3 <- rep(NA, M)
hubs <- mod3$hub
tic()
for(i in 1:M){
  xx <- XX[i,]
  tmp <- slapGP(xx, X, y, hubs=hubs)
  yhat3[i] <- tmp$pred
  hubs <- tmp$hubs
}
t3 <- toc()

# laGP predictions
yhat4 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat4[i] <- laGP(matrix(xx, nrow=1), start=6, end=ceiling(sqrt(nrow(X))),
                   X=X, Z=matrix(y, ncol=1))$mean
}
t4 <- toc()


# BART predictions
yhat5 <- rep(NA, M)
tic()
for(i in 1:M){
  xx <- XX[i,]
  yhat5[i] <- mean(predict(mod5, matrix(xx, ncol=3)))
}
t5 <- toc()

# slagpGP(0.9999) predictions
mod6 <- leapGP_synch(mod2, rho=0.9999)
yhat6 <- rep(NA, M)
hubs <- mod6$hub
tic()
for(i in 1:M){
  xx <- XX[i,]
  tmp <- slapGP(xx, X, y, hubs=hubs, rho=0.9999)
  yhat6[i] <- tmp$pred
  hubs <- tmp$hubs
}
t6 <- toc()


#Make a table
TAB <- matrix(NA, nrow=5, ncol=3)
rownames(TAB) <- c("laGP", "leapGP", "slapGP", "BMARS", "BART")
colnames(TAB) <- c("RMSE", "time to train", "time to predict")
TAB[1,] <- round(c(sqrt(sum((yhat4-yy)^2)), 0, t4$toc - t4$tic), 3)
TAB[2,] <- round(c(sqrt(sum((yhat2-yy)^2)), t2b$toc - t2b$tic, t2$toc - t2$tic), 3)
TAB[3,] <- round(c(sqrt(sum((yhat3-yy)^2)), t2b$toc - t2b$tic, t3$toc - t3$tic), 3)
TAB[4,] <- round(c(sqrt(sum((yhat1-yy)^2)), t1b$toc - t1b$tic, t1$toc - t1$tic), 3)
TAB[5,] <- round(c(sqrt(sum((yhat5-yy)^2)), t5b$toc - t5b$tic, t5$toc - t5$tic), 3)

save(TAB, file=paste0("data/GP_comp_tab_", N,".rda"))

