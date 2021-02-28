## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----library------------------------------------------------------------------------------------------------------------------------------
library("RBF")

## ----read the dataset---------------------------------------------------------------------------------------------------------------------
datos <- read.csv("housing.csv", header = FALSE, sep="")
x1 <- log(datos$V1)
x2 <- log(datos$V3)
x3 <- log(datos$V5)
x4 <- log(datos$V6)
x5 <- log(datos$V7)
x6 <- log(datos$V8)
x7 <- log(datos$V10)
x8 <- log(datos$V11)
x9 <- log(datos$V12)
x10 <- log(datos$V13)
y <- datos$V14

## ----bandwidths---------------------------------------------------------------------------------------------------------------------------
X <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
colnames(X) <- NULL #HAY QUE ARREGLAR LOS NOMBRES CUANDO VIENEN CON NOMBRES EN EL .R DEL PAQUETE
bandw <- (1/2)*apply(X,2,sd) 
bandw

## ----robust fit---------------------------------------------------------------------------------------------------------------------------
robust.fit <- backf.rob(y~X, windows=bandw)

## ----summary------------------------------------------------------------------------------------------------------------------------------
summary(robust.fit)

## ----plot---------------------------------------------------------------------------------------------------------------------------------
plot(robust.fit)

## ----prediction---------------------------------------------------------------------------------------------------------------------------
po <- colMeans(X)
robust.fit1 <- backf.rob(y~X, windows = bandw, point = po)
robust.fit1$prediction

## ----outliers-----------------------------------------------------------------------------------------------------------------------------
ynew <- y
#ynew[c(4,89,197,198,199,200,299,291,301,350)]<- rep(400,10)
ynew[1:5]<- rep(400,5)

## ----robustplotswithoutliers--------------------------------------------------------------------------------------------------------------
robust.fit.new <- backf.rob(ynew~X, windows = bandw, point = po)
summary(robust.fit.new)
robust.fit.new$prediction
plot(robust.fit.new)

## ----robustplots2, warning=FALSE----------------------------------------------------------------------------------------------------------
for(j in 1:10) {
  name.x <- bquote(paste('x')[.(j)]) 
  name.y <- bquote(paste(hat('g')[.(j)]))
  oo <- order(X[,j])
  plot(X[oo,j], robust.fit.new$g.matrix[oo,j], type="l", lwd=5, col='blue', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(X[oo,j], robust.fit$g.matrix[oo,j], lwd=5, col='green', lty=2)
}

## ----gam, warning=FALSE-------------------------------------------------------------------------------------------------------------------
library(gam)
dataset <- as.data.frame(cbind(y,X))
colnames(dataset) <- c("y","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
fit.gam <- gam(y ~ lo(x1, span=1.62)+lo(x2, span=0.58)+lo(x3, span=0.15)+lo(x4, span=0.08)+
                 lo(x5, span=0.46)+lo(x6, span=0.40)+lo(x7, span=0.30)+lo(x8, span=0.09)+
                 lo(x9, span=0.58)+lo(x10, span=0.45), data=dataset) #gam(y ~ lo(x1, span=2.16)+lo(x2, span=0.78)+lo(x3, span=0.20)+lo(x4, span=0.12)+ lo(x5, span=0.62)+lo(x6, span=0.54)+lo(x7, span=0.40)+lo(x8, span=0.12)+ lo(x9, span=0.77)+lo(x10, span=0.60) , data=dataset)
fits <- predict(fit.gam, type='terms')
dataset.new <- as.data.frame(cbind(ynew,X))
colnames(dataset.new) <- c("ynew","x1","x2","x3","x4","x5","x6","x7","x8","x9","x10")
fit.gam.new <- gam(ynew ~  lo(x1, span=1.62)+lo(x2, span=0.58)+lo(x3, span=0.15)+
                     lo(x4, span=0.08)+lo(x5, span=0.46)+lo(x6, span=0.40)+lo(x7, span=0.30)+
                     lo(x8, span=0.09)+ lo(x9, span=0.58)+lo(x10, span=0.45), data=dataset.new)
fits.new <- predict(fit.gam.new, type='terms')

## ----gamplots-----------------------------------------------------------------------------------------------------------------------------
for(j in 1:10) {
  oo <- order(X[,j])
  name.x <- bquote(paste('x')[.(j)])
  name.y <- bquote(paste(hat('g')[.(j)]))
  plot(X[oo,j], fits.new[oo,j], type="l", lwd=5, col='purple', lty=1, 
       xlab=name.x, ylab=name.y)
  lines(X[oo,j], fits[oo,j], lwd=5, col='darkorange2', lty=2)
}

