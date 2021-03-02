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

