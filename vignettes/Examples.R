## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----read the dataset----------------------------------------------------------------------------------------------------------------------
data(Boston, package='MASS')
dd <- Boston[, c(1, 3, 5:8, 10:14)]
dd[, -11] <- log( dd[, names(dd) != 'medv'] )

## ----loadpckg------------------------------------------------------------------------------------------------------------------------------
library(RBF)

## ----bandwidths----------------------------------------------------------------------------------------------------------------------------
bandw <- apply(dd[, names(dd) != 'medv'], 2, sd) / 2

## ----summary-------------------------------------------------------------------------------------------------------------------------------
summary(robust.fit)

## ----plot----------------------------------------------------------------------------------------------------------------------------------
plot(robust.fit)

## ----showpred------------------------------------------------------------------------------------------------------------------------------
robust.fit1$prediction

## ----outliers------------------------------------------------------------------------------------------------------------------------------
dd2 <- dd
dd2$medv[1:5]<- rep(400, 5)

## ----robustplots2, warning=FALSE-----------------------------------------------------------------------------------------------------------
for(j in 1:10) {
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  oo <- order(dd2[,j])
  plot(dd2[oo,j], robust.fit.new$g.matrix[oo,j], type="l", lwd=5, col='blue', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], robust.fit$g.matrix[oo,j], lwd=5, col='green', lty=2)
}

## ----gam, warning=FALSE--------------------------------------------------------------------------------------------------------------------
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

## ----gamplots------------------------------------------------------------------------------------------------------------------------------
for(j in 1:10) {
  oo <- order(dd2[,j])
  name.x <- names(dd)[j] 
  name.y <- bquote(paste(hat('g')[.(j)]))
  plot(dd2[oo,j], fits.new[oo,j], type="l", lwd=5, col='purple', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(dd2[oo,j], fits[oo,j], lwd=5, col='darkorange2', lty=2)
}

## ----scatterplot---------------------------------------------------------------------------------------------------------------------------
data(airquality)
ccs <- complete.cases(airquality)
aircomplete <- airquality[ccs, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
pairs(aircomplete[, c('Ozone', 'Solar.R', 'Wind', 'Temp')], pch=19, col='gray30', cex=1.5)

## ----bandw---------------------------------------------------------------------------------------------------------------------------------
bandw <- c(136.7285, 10.67314, 4.764985)

## ----fitfull-------------------------------------------------------------------------------------------------------------------------------
fit.full <- backf.rob(Ozone ~ Solar.R + Wind + Temp, windows = bandw, 
                      epsilon = 1e-6, degree = 1, type = 'Tukey', 
                      subset = ccs, data = airquality)

## ----plotfitfull---------------------------------------------------------------------------------------------------------------------------
plot(fit.full)

## ----fitgam--------------------------------------------------------------------------------------------------------------------------------
fit.gam <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.7)+
                 lo(Temp, span=.5), data = aircomplete)

## ----plotrobgam----------------------------------------------------------------------------------------------------------------------------
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  #points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta', lty=2)
}

## ----boxplot-------------------------------------------------------------------------------------------------------------------------------
re.ro <- residuals(fit.full)
ou.ro <- boxplot(re.ro, col='gray80')$out
in.ro <- (1:length(re.ro))[ re.ro %in% ou.ro ]
points(rep(1, length(in.ro)), re.ro[in.ro], pch=20, cex=3, col='red')
(in.ro)

## ----scatterplotpoints---------------------------------------------------------------------------------------------------------------------
cs <- rep('gray30', nrow(aircomplete))
cs[in.ro] <- 'red'
os <- 1:nrow(aircomplete)
os2 <- c(os[-in.ro], os[in.ro])
pairs(aircomplete[os2, c('Ozone', 'Solar.R', 'Wind', 'Temp')], 
      pch=19, col=cs[os2], cex=1.5)

## ----plotoutred----------------------------------------------------------------------------------------------------------------------------
# Plot both fits (robust and classical) 
x <- as.matrix( aircomplete[ , c('Solar.R', 'Wind', 'Temp')] )
y <- as.vector( aircomplete[ , 'Ozone'] )
fits <- predict(fit.gam, type='terms')
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
  lines(x[oo,j], fits[oo,j], lwd=5, col='magenta', lty=2)
}

## ----cvgamclean, eval=FALSE----------------------------------------------------------------------------------------------------------------
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

## ----fitgam2-------------------------------------------------------------------------------------------------------------------------------
airclean <- aircomplete[-in.ro, c('Ozone', 'Solar.R', 'Wind', 'Temp')]
fit.gam2 <- gam(Ozone ~ lo(Solar.R, span=.7) + lo(Wind, span=.8)+
                  lo(Temp, span=.3), data=airclean) 

## ----finalplot-----------------------------------------------------------------------------------------------------------------------------
fits2 <- predict(fit.gam2, type='terms')
dd2 <- aircomplete[-in.ro, c('Solar.R', 'Wind', 'Temp')]
for(j in 1:3) {
  re <- fit.full$yp - fit.full$alpha - rowSums(fit.full$g.matrix[,-j])
  plot(re ~ x[,j], type='p', pch=20, col='gray45', xlab=colnames(x)[j], ylab='')
  points(re[in.ro] ~ x[in.ro,j], pch=20, col='red', cex=2)
  oo <- order(dd2[,j])
  lines(dd2[oo,j], fits2[oo,j], lwd=5, col='magenta', lty=2)
  oo <- order(x[,j])
  lines(x[oo,j], fit.full$g.matrix[oo,j], lwd=5, col='blue', lty=1)
}

## ----pred1---------------------------------------------------------------------------------------------------------------------------------
tms <- function(a, alpha=.1) {
  # alpha is the proportion to trim
  a2 <- sort(a^2, na.last=NA)
  n0 <- floor( length(a) * (1 - alpha) )
  return( mean(a2[1:n0], na.rm=TRUE) )
}

