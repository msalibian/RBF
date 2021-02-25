## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----library------------------------------------------------------------------
library("RBF")

## ----read the dataset---------------------------------------------------------
datos <- read.csv("housing.csv", header = FALSE, sep="")
x1 <- datos$V1
x2 <- datos$V3
x3 <- datos$V5
x4 <- datos$V6
x5 <- datos$V7
x6 <- datos$V8
x7 <- datos$V10
x8 <- datos$V11
x9 <- datos$V12
x10 <- datos$V13
X <- log(cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10))
colnames(X) <- NULL
y <- datos$V14

## ----bandwidths---------------------------------------------------------------
bandw <- 1/2*apply(X,2,sd)
bandw

## ----robust fit---------------------------------------------------------------
robust.fit <- backf.rob( y~X, windows=bandw)

## ----summary------------------------------------------------------------------
summary(robust.fit)

## ----plot---------------------------------------------------------------------
plot(robust.fit)

