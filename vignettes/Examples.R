## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----generating data----------------------------------------------------------
function.g1 <- function(x1) 24*(x1-1/2)^2-2
function.g2 <- function(x2) 2*pi*sin(pi*x2)-4
set.seed(140)
n <- 500
x1 <- runif(n)
x2 <- runif(n)
X <- cbind(x1, x2)
eps <- rnorm(n,0,sd=0.15)
regression <- function.g1(x1)+function.g2(x2)
y <- regression + eps

## ----bandw1-------------------------------------------------------------------
bandw <- c(0.05, 0.075)

## ----robust fit---------------------------------------------------------------
library(RBF)
point <- c(0.7, 0.6)
robust.fit <- backf.rob(y ~ X, point=point, windows=bandw, type = 'Tukey', degree=1)

## ----prediction---------------------------------------------------------------
robust.fit$prediction
c(function.g1(point[1]), function.g2(point[2]))

## ----plots-simu, fig.show="hold", out.width="33%"-----------------------------
lim.rob <- matrix(0, 2, 2)
functions.g <- cbind(function.g1(X[,1]), function.g2(X[,2]))
for(j in 1:2) {
  res <- y - robust.fit$alpha - robust.fit$g.matrix[,-j]
  lim.rob[,j] <- range(res)
  plot(X[,j], res, type='p', pch=19, col='gray45', xlab=colnames(X)[j], ylab='', 
       cex=1, ylim=lim.rob[,j])
  ord <- order(X[,j])
  lines(X[ord,j], robust.fit$g.matrix[ord,j], lwd=3, col='blue')
  lines(X[ord,j], functions.g[ord,j], lwd=3)
}

## ----contaminated responses---------------------------------------------------
eps <- rnorm(n,0,sd=0.15)
prop.cont <- 0.10
ou <- rbinom(n, size=1, prob=prop.cont) 
eps.ou <- eps
eps.ou[ ou == 1 ] <- rnorm(sum(ou),mean=15, sd=0.1)
yout <- regression + eps.ou

## ----bandw-out----------------------------------------------------------------
bandw <- c(0.05, 0.075)

## ----robust.fit.out-----------------------------------------------------------
robust.fit.out <- backf.rob(yout ~ X, point=point, windows=bandw, type = 'Tukey', degree=1)

## ----prediction2--------------------------------------------------------------
robust.fit.out$prediction
c(function.g1(point[1]), function.g2(point[2]))

## ----plots-simu2, fig.show="hold", out.width="33%"----------------------------
lim.rob <- matrix(0, 2, 2)
functions.g <- cbind(function.g1(X[,1]), function.g2(X[,2]))
for(j in 1:2) {
  res <- yout - robust.fit.out$alpha - robust.fit.out$g.matrix[,-j]
  lim.rob[,j] <- range(res)
  plot(X[,j], res, type='p', pch=19, col='gray45', xlab=colnames(X)[j], ylab='', 
       cex=1, ylim=lim.rob[,j])
  ord <- order(X[,j])
  lines(X[ord,j], robust.fit.out$g.matrix[ord,j], lwd=3, col='blue')
  lines(X[ord,j], functions.g[ord,j], lwd=3)
}

## ----rkfoldfunction-----------------------------------------------------------
backf.rob.cv <- function(k=5, Xp, yp, windows, degree, type) {
  n <- length(yp)
  k1 <- floor(n/k)
  ids <- rep(1:k, each=k1)
  if( length(ids) < n ) ids <- c(ids, 1:(n%%k))
  ids <- sample(ids)
  preds <- rep(NA, n)
  for(j in 1:k) {
    XX <- Xp[ids!=j,]
    yy <- yp[ids!=j]
    tmp <- backf.rob(yy ~ XX, point=Xp[ids==j,], windows=windows, degree=degree, type=type)
    preds[ids==j] <- rowSums(tmp$prediction) + tmp$alpha
  }
  preds.res <- preds - yp
  tmp.re <- RobStatTM::locScaleM(preds.res, na.rm=TRUE)
  sal <- tmp.re$mu^2 + tmp.re$disper^2
  return(sal)
}

## ----robust.5fold1, eval=FALSE------------------------------------------------
#  # Bandwidth selection with leave-one-out cross-validation
#  # This takes some time to compute (approx 2 minutes running
#  # R 3.6.3 on an Intel(R) Core(TM) i7-10900 CPU @ 2.90GHz)
#  h1 <- c(0.05, 0.075, 0.1, 0.125)
#  hh <- expand.grid(h1, h1)
#  nh <- nrow(hh)
#  rmspe <- rep(NA, nh)
#  system.time({
#  for(i in 1:nh){
#    print(i)
#    rmspe[i] <- backf.rob.cv(k=5, Xp=X, yp=y, windows=hh[i,], degree=1, type='Tukey')
#  }
#  })
#  i0 <- which.min(rmspe)
#  bandw <- hh[i0,]

## ----optbandw1----------------------------------------------------------------
bandw

## ----robust.5fold2, eval=FALSE------------------------------------------------
#  # Bandwidth selection with leave-one-out cross-validation
#  # This takes some time to compute (approx 2 minutes running
#  # R 3.6.3 on an Intel(R) Core(TM) i7-10900 CPU @ 2.90GHz)
#  h1 <- c(0.05, 0.075, 0.1, 0.125)
#  hh <- expand.grid(h1, h1)
#  nh <- nrow(hh)
#  rmspe <- rep(NA, nh)
#  system.time({
#  for(i in 1:nh){
#    print(i)
#    rmspe[i] <- backf.rob.cv(k=5, Xp=X, yp=yout, windows=hh[i,], degree=1, type='Tukey')
#  }
#  })
#  i0 <- which.min(rmspe)
#  bandw <- hh[i0,]

## ----optbandw2----------------------------------------------------------------
bandw

## ----intro--------------------------------------------------------------------
library(RBF)
data(airquality)
pairs(airquality[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col='gray30', cex=1.5)

## ----robust.leaveoneout, eval=FALSE-------------------------------------------
#  # Bandwidth selection with leave-one-out cross-validation
#  ## Without outliers
#  # This takes a long time to compute (approx 380 minutes running
#  # R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
#  ccs <- complete.cases(airquality)
#  x <- as.matrix( airquality[ccs, c('Solar.R', 'Wind', 'Temp')] )
#  y <- as.vector( airquality[ccs, 'Ozone'] )
#  a <- c(1/2, 1, 1.5, 2, 2.5, 3)
#  h1 <- a * sd(x[,1])
#  h2 <- a * sd(x[,2])
#  h3 <- a * sd(x[,3])
#  hh <- expand.grid(h1, h2, h3)
#  nh <- nrow(hh)
#  rmspe <- rep(NA, nh)
#  jbest <- 0
#  cvbest <- +Inf
#  # leave-one-out
#  n <- nrow(x)
#  for(i in 1:nh) {
#    # leave-one-out CV loop
#    preds <- rep(NA, n)
#    for(j in 1:n) {
#      tmp <- try( backf.rob(y ~ x, point = x[j, ],
#                            windows = hh[i, ], epsilon = 1e-6,
#                            degree = 1, type = 'Tukey', subset = c(-j) ))
#      if (class(tmp)[1] != "try-error") {
#        preds[j] <- rowSums(tmp$prediction) + tmp$alpha
#      }
#    }
#    pred.res <- preds - y
#    tmp.re <- RobStatTM::locScaleM(pred.res, na.rm=TRUE)
#    rmspe[i] <- tmp.re$mu^2 + tmp.re$disper^2
#    if( rmspe[i] < cvbest ) {
#      jbest <- i
#      cvbest <- rmspe[i]
#      print('Record')
#    }
#    print(c(i, rmspe[i]))
#  }
#  bandw <- hh[jbest,]

## ----bandw--------------------------------------------------------------------
bandw <- c(136.7285, 10.67314, 4.764985)

## ----rbfone-------------------------------------------------------------------
ccs <- complete.cases(airquality)
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, data=airquality,
                subset=ccs, windows=bandw, degree=1, type='Tukey')

## ----showfits, fig.show="hold", out.width="33%"-------------------------------
lim.cl <- lim.rob <- matrix(0, 2, 3)
x0 <- fit.full$Xp
for(j in 1:3) {
  re <- fit.full$y - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  lim.rob[,j] <- c(min(re), max(re))
  plot(re ~ x0[,j], type='p', pch=19, col='gray30', 
       xlab=colnames(x0)[j], ylab='', cex=1.5)
  oo <- order(x0[,j])
  lines(x0[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}

## ----plot.method, eval=FALSE--------------------------------------------------
#  # Plot each component
#  plot(fit.full, which=1:3)

## ----classicfits, fig.show="hold", out.width="33%"----------------------------
aircomplete <- airquality[ complete.cases(airquality), ]
library(gam)
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data=aircomplete)
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
# alpha.gam <- attr(fits, 'constant')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}

## ----outliers, fig.width=4, fig.height=4--------------------------------------
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=1.5, col='red')

## ----showouts-----------------------------------------------------------------
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)

## ----showouts2, fig.show="hold", out.width="33%"------------------------------
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x0[in.ro,j], pch=19, col='red', cex=1.5)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta')
}

## ----bothonclean, fig.show="hold", out.width="33%"----------------------------
# Run the classical backfitting algorithm without outliers
airclean <- aircomplete[-in.ro, ]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean)
fits2 <- predict(fit.gam2, type='terms')
# alpha.gam2 <- attr(fits2, 'constant')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=1.5)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue')
}

