## ----setup, include=FALSE----------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----read the dataset--------------------------------------------
data(Boston, package='MASS')
dd <- Boston[, c(1, 3, 5:8, 10:14)]
dd[, names(dd) != 'medv'] <- log( dd[, names(dd) != 'medv'] )

## ----loadpckg----------------------------------------------------
library(RBF)

## ----bandwidths--------------------------------------------------
bandw <- apply(dd[, names(dd) != 'medv'], 2, sd) / 2

## ----robustfit, cache=TRUE---------------------------------------
robust.fit <- backf.rob(medv ~ ., data = dd, degree = 0, type = 'Huber', 
                        windows = bandw)

## ----summary-----------------------------------------------------
summary(robust.fit)

## ----plot, out.width  = "45%"------------------------------------
plot(robust.fit)

## ----preds, cache=TRUE-------------------------------------------
po <- colMeans(dd[, names(dd) != 'medv'])
robust.fit1 <- backf.rob(medv ~ ., data = dd, degree = 0, type = 'Huber', 
                         windows = bandw, point = po)

## ----showpred----------------------------------------------------
robust.fit1$prediction

## ----outliers----------------------------------------------------
dd2 <- dd
dd2$medv[1:5]<- rep(400, 5)

## ----robustplotswithoutliers, cache=TRUE-------------------------
robust.fit.new <- backf.rob(medv ~ ., data = dd2, degree = 0, type = 'Huber', 
                            windows = bandw, point = po)
summary(robust.fit.new)
robust.fit.new$prediction

## ----robustplots2, warning=FALSE, out.width  = "45%"-------------
for(j in 1:10) {
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  oo <- order(dd2[,j])
  plot(dd2[oo,j], robust.fit.new$g.matrix[oo,j], type="l", lwd=2, col='blue', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], robust.fit$g.matrix[oo,j], lwd=2, col='green', lty=2)
}

## ----gam, warning=FALSE------------------------------------------
library(gam)
fit.gam <- gam(medv ~ lo(crim, span=1.62) + 
                 lo(indus, span=0.58) + 
                 lo(nox, span=0.15) + 
                 lo(rm, span=0.08) +
                 lo(age, span=0.46) + 
                 lo(dis, span=0.40) + 
                 lo(tax, span=0.30) + 
                 lo(ptratio, span=0.09) +
                 lo(black, span=0.58) + 
                 lo(lstat, span=0.45), data=dd)
fits <- predict(fit.gam, type='terms')
fit.gam.new <- gam(medv ~ lo(crim, span=1.62) + 
                 lo(indus, span=0.58) + 
                 lo(nox, span=0.15) + 
                 lo(rm, span=0.08) +
                 lo(age, span=0.46) + 
                 lo(dis, span=0.40) + 
                 lo(tax, span=0.30) + 
                 lo(ptratio, span=0.09) +
                 lo(black, span=0.58) + 
                 lo(lstat, span=0.45), data=dd2)
fits.new <- predict(fit.gam.new, type='terms')

## ----gamplots, out.width  = "45%"--------------------------------
for(j in 1:10) {
  oo <- order(dd2[,j])
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  plot(dd2[oo,j], fits.new[oo,j], type="l", lwd=2, col='purple', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], fits[oo,j], lwd=2, col='darkorange2', lty=2)
}

## ----scatterplot, out.width  = "75%"-----------------------------
data(airquality)
ccs <- complete.cases(airquality)
aircomplete <- airquality[ccs, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
pairs(aircomplete[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], pch=19, col='gray30')

## ----robustcv, warning=FALSE, cache=TRUE, eval=FALSE-------------
#  library(RBF)
#  
#  # Bandwidth selection with leave-one-out cross-validation
#  ## Without outliers
#  # This takes a long time to compute (approx 380 minutes running
#  # R 3.6.1 on an Intel(R) Core(TM) i7-4790 CPU @ 3.60GHz)
#  a <- c(1/2, 1, 1.5, 2, 2.5, 3)
#  h1 <- a * sd(aircomplete[,2])
#  h2 <- a * sd(aircomplete[,3])
#  h3 <- a * sd(aircomplete[,4])
#  hh <- expand.grid(h1, h2, h3)
#  nh <- nrow(hh)
#  rmspe <- rep(NA, nh)
#  jbest <- 0
#  cvbest <- +Inf
#  n <- nrow(aircomplete)
#  for(i in 1:nh) {
#    # leave-one-out CV loop
#    preds <- rep(NA, n)
#    for(j in 1:n) {
#      tmp <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp, point = aircomplete[j, -1],
#                            windows = hh[i, ], epsilon = 1e-6, data = aircomplete,
#                            degree = 1, type = 'Tukey', subset = c(-j) ))
#      if (class(tmp)[1] != "try-error") {
#        preds[j] <- rowSums(tmp$prediction) + tmp$alpha
#      }
#    }
#    tmp.re <- RobStatTM::locScaleM(preds - aircomplete$Ozone, na.rm=TRUE)
#    rmspe[i] <- tmp.re$mu^2 + tmp.re$disper^2
#    if( rmspe[i] < cvbest ) {
#      jbest <- i
#      cvbest <- rmspe[i]
#    }
#  }
#  (bandw <- hh[jbest,])

## ----bandw-------------------------------------------------------
bandw <- c(136.7285, 10.67314, 4.764985)

## ----fitfull-----------------------------------------------------
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, windows = bandw, 
                      epsilon = 1e-6, degree = 1, type = 'Tukey', 
                      subset = ccs, data = airquality)

## ----plotfitfull, out.width  = "45%"-----------------------------
plot(fit.full)

## ----gamcv, cache=TRUE, warning=FALSE, eval=FALSE----------------
#  library(gam)
#  a <- c(.3, .4, .5, .6, .7, .8, .9)
#  hh <- expand.grid(a, a, a)
#  nh <- nrow(hh)
#  jbest <- 0
#  cvbest <- +Inf
#  n <- nrow(aircomplete)
#  for(i in 1:nh) {
#    fi <- rep(0, n)
#    for(j in 1:n) {
#      tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
#                 + lo(Temp, span=hh[i,3]), data = aircomplete, subset=c(-j))
#      fi[j] <- as.numeric(predict(tmp, newdata=aircomplete[j, -1], type='response'))
#    }
#    ss <- mean((aircomplete$Ozone - fi)^2)
#    if(ss < cvbest) {
#      jbest <- i
#      cvbest <- ss
#    }
#  }
#  (hh[jbest,])
#  # Var1 Var2 Var3
#  # 131  0.7  0.7  0.5

## ----fitgam------------------------------------------------------
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data = aircomplete)

## ----plotrobgam, out.width  = "45%"------------------------------
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=2, col='magenta', lty=2)
}

## ----boxplot, out.width  = "50%"---------------------------------
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, col='red')
(in.ro)

## ----scatterplotpoints, out.width  = "75%"-----------------------
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2])

## ----plotoutred, out.width  = "45%"------------------------------
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red')
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=2, col='magenta', lty=2)
}

## ----cvgamclean, eval=FALSE--------------------------------------
#  airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
#  a <- c(.3, .4, .5, .6, .7, .8, .9)
#  hh <- expand.grid(a, a, a)
#  nh <- nrow(hh)
#  jbest <- 0
#  cvbest <- +Inf
#  n <- nrow(airclean)
#  for(i in 1:nh) {
#    fi <- rep(0, n)
#    for(j in 1:n) {
#      tmp <- gam(Ozone ~ lo(Solar.R, span=hh[i,1]) + lo(Wind, span=hh[i,2])
#                 + lo(Temp, span=hh[i,3]), data=airclean, subset=c(-j))
#      fi[j] <- as.numeric(predict(tmp, newdata=airclean[j,], type='response'))
#    }
#    ss <- mean((airclean$Ozone - fi)^2)
#    if(ss < cvbest) {
#      jbest <- i
#      cvbest <- ss
#    }
#  }
#  (hh[jbest,])
#  # # Var1 Var2 Var3
#  # # 40  0.7  0.8  0.3

## ----fitgam2-----------------------------------------------------
airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean) 

## ----finalplot, out.width  = "45%"-------------------------------
fits2 <- predict(fit.gam2, type='terms')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red')
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=2, col='magenta', lty=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=2, col='blue', lty=1)
}

## ----pred1-------------------------------------------------------
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}

## ----pred2, cache=TRUE, warning=FALSE, eval=FALSE----------------
#  dd <- airquality
#  dd <- dd[complete.cases(dd), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
#  # 100 runs of K-fold CV
#  M <- 100
#  # 5-fold
#  K <- 5
#  n <- nrow(dd)
#  # store (trimmed) TMSPE for robust and gam, and also
#  tmspe.ro <- tmspe.gam <- vector('numeric', M)
#  set.seed(123)
#  ii <- (1:n)%%K + 1
#  for(runs in 1:M) {
#    tmpro <- tmpgam <- vector('numeric', n)
#    ii <- sample(ii)
#    for(j in 1:K) {
#      fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp,
#                             point=dd[ii==j, -1], windows = bandw,
#                             epsilon = 1e-6, degree = 1, type = 'Tukey',
#                             subset = (ii!=j), data = dd)
#      tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
#      fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
#                       lo(Temp, span=.5), data = dd[ii!=j, ])
#      tmpgam[ ii == j ] <- predict(fit.gam, newdata=dd[ii==j, ], type='response')
#    }
#    tmspe.ro[runs] <- tms( dd$Ozone - tmpro, alpha=0.05)
#    tmspe.gam[runs] <- tms( dd$Ozone - tmpgam, alpha=0.05)
#  }

## ----pred.clean, cache=TRUE, warning=FALSE, eval=FALSE-----------
#  aq <- airquality
#  aq2 <- aq[complete.cases(aq), c('Ozone', 'Solar.R', 'Wind', 'Temp')]
#  airclean <- aq2[ -in.ro, ]
#  bandw <- c(138.2699, 10.46753, 4.828436)
#  M <- 100
#  K <- 5
#  n <- nrow(airclean)
#  mspe.ro <- mspe.gam <- tmspe.ro <- tmspe.gam <- vector('numeric', M)
#  set.seed(17)
#  ii <- (1:n)%%K + 1
#  for(runs in 1:M) {
#    tmpro <- tmpgam <- vector('numeric', n)
#    ii <- sample(ii)
#    for(j in 1:K) {
#      fit.full <- try( backf.rob(Ozone ~ Solar.R + Wind + Temp,
#                            point=airclean[ii==j, -1], windows = bandw,
#                            epsilon = 1e-6, degree = 1, type = 'Tukey',
#                            subset = (ii!=j), data = airclean) )
#      if (class(fit.full)[1] != "try-error") {
#        tmpro[ ii == j ] <- rowSums(fit.full$prediction) + fit.full$alpha
#      }
#      fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
#                       lo(Temp, span=.3), data = airclean, subset = (ii!=j) )
#      tmpgam[ ii == j ] <- predict(fit.gam, newdata=airclean[ii==j, ],
#                                   type='response')
#    }
#    tmspe.ro[runs] <- tms( airclean$Ozone - tmpro, alpha=0.05)
#    mspe.ro[runs] <- mean( ( airclean$Ozone - tmpro)^2, na.rm=TRUE)
#    tmspe.gam[runs] <- tms( airclean$Ozone - tmpgam, alpha=0.05)
#    mspe.gam[runs] <- mean( ( airclean$Ozone - tmpgam)^2, na.rm=TRUE)
#  }

