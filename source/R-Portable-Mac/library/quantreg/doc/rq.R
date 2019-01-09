### R code from vignette source 'rq.Rnw'

###################################################
### code chunk number 1: R options
###################################################
options(width = 60)
options(SweaveHooks = list(fig = function() par(mar=c(3,3,1,0.5),mgp = c(2,1,0))))


###################################################
### code chunk number 2: rq.Rnw:126-127
###################################################
library(quantreg)


###################################################
### code chunk number 3: rq.Rnw:137-139 (eval = FALSE)
###################################################
## help(package="quantreg")
## help(rq)


###################################################
### code chunk number 4: rq.Rnw:183-185
###################################################
data(engel)
fit1 <- rq(foodexp ~ income, tau = .5, data = engel)


###################################################
### code chunk number 5: rq.Rnw:197-198
###################################################
fit1


###################################################
### code chunk number 6: rq.Rnw:204-205
###################################################
summary(fit1)


###################################################
### code chunk number 7: rq.Rnw:213-215
###################################################
r1 <- resid(fit1)
c1 <- coef(fit1)


###################################################
### code chunk number 8: rq.Rnw:269-270
###################################################
summary(fit1,se = "nid")


###################################################
### code chunk number 9: engelplot
###################################################
getOption("SweaveHooks")[["fig"]]()
library(quantreg)
data(engel)
attach(engel)
plot(income,foodexp,cex=.25,type="n",xlab="Household Income", ylab="Food Expenditure")
points(income,foodexp,cex=.5,col="blue")
abline(rq(foodexp~income,tau=.5),col="blue")
abline(lm(foodexp~income),lty=2,col="red") #the dreaded ols line
taus <- c(.05,.1,.25,.75,.90,.95)
for( i in 1:length(taus)){
        abline(rq(foodexp~income,tau=taus[i]),col="gray")
        }


###################################################
### code chunk number 10: engelcoef
###################################################
xx <- income - mean(income)
fit1 <- summary(rq(foodexp~xx,tau=2:98/100))
fit2 <- summary(rq(foodexp~xx,tau=c(.05, .25, .5, .75, .95)))


###################################################
### code chunk number 11: engelcoefplot
###################################################
pdf("engelcoef.pdf",width=6.5,height=3.5)
plot(fit1,mfrow = c(1,2))
dev.off()


###################################################
### code chunk number 12: engeltable
###################################################
latex(fit2, caption="Engel's Law", transpose=TRUE)


###################################################
### code chunk number 13: rqProcess
###################################################
z <- rq(foodexp~income,tau=-1)


###################################################
### code chunk number 14: eqfs
###################################################
getOption("SweaveHooks")[["fig"]]()
x.poor <- quantile(income,.1) #Poor is defined as at the .1 quantile of the sample distn
x.rich <- quantile(income,.9) #Rich is defined as at the .9 quantile of the sample distn
ps <- z$sol[1,]
qs.poor <- c(c(1,x.poor)%*%z$sol[4:5,])
qs.rich <- c(c(1,x.rich)%*%z$sol[4:5,])
#now plot the two quantile functions to compare
par(mfrow = c(1,2))
plot(c(ps,ps),c(qs.poor,qs.rich), type="n",
     xlab = expression(tau), ylab = "quantile")
plot(stepfun(ps,c(qs.poor[1],qs.poor)), do.points=FALSE, add=TRUE)
plot(stepfun(ps,c(qs.poor[1],qs.rich)), do.points=FALSE, add=TRUE,
     col.hor = "gray", col.vert = "gray")
## now plot associated conditional density estimates
## weights from ps (process)
ps.wts <- (c(0,diff(ps)) + c(diff(ps),0)) / 2
ap <- akj(qs.poor, z=qs.poor, p = ps.wts)
ar <- akj(qs.rich, z=qs.rich, p = ps.wts)
plot(c(qs.poor,qs.rich),c(ap$dens,ar$dens),type="n",
     xlab= "Food Expenditure", ylab= "Density")
lines(qs.rich, ar$dens, col="gray")
lines(qs.poor, ap$dens, col="black")
legend("topright", c("poor","rich"), lty = c(1,1), col=c("black","gray"))


###################################################
### code chunk number 15: engellogplot
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(income,foodexp,log="xy",xlab="Household Income", ylab="Food Expenditure")
taus <- c(.05,.1,.25,.75,.90,.95)
abline(rq(log10(foodexp)~log10(income),tau=.5),col="blue")
abline(lm(log10(foodexp)~log10(income)),lty = 3,col="red")
for( i in 1:length(taus)){
       abline(rq(log10(foodexp)~log10(income),tau=taus[i]),col="gray")
       }


###################################################
### code chunk number 16: rq.Rnw:550-553
###################################################
fit1 <- rq(foodexp~income,tau=.25)
fit2 <- rq(foodexp~income,tau=.50)
fit3 <- rq(foodexp~income,tau=.75)


###################################################
### code chunk number 17: rq.Rnw:563-564
###################################################
anova(fit1, fit2, fit3)


###################################################
### code chunk number 18: gastest
###################################################
source("gasprice.R")
x <- gasprice
n <- length(x)
p <- 5 # lag length
X <- cbind(x[(p-1):(n-1)], x[(p-2):(n-2)], x[(p-3):(n-3)], x[(p-4):(n-4)])
y <- x[p:n]
T1 <- KhmaladzeTest(y ~ X,taus = -1, nullH="location")
T2 <- KhmaladzeTest(y ~ X,taus = 10:290/300, nullH="location",se="ker")


###################################################
### code chunk number 19: Frank
###################################################
	n <- 200
	df <- 8
	delta <- 8
	set.seed(4003)
	x <- sort(rt(n,df))
	u <- runif(n)
	v <- -log(1-(1-exp(-delta))/(1+exp(-delta*pt(x,df))*((1/u)-1)))/delta
	y <- qt(v,df)


###################################################
### code chunk number 20: Frankplot
###################################################
getOption("SweaveHooks")[["fig"]]()
	plot(x,y,col="blue",cex = .25)
	us <- c(.25,.5,.75)
	for(i in 1:length(us)){
		u <- us[i]
		v <- -log(1-(1-exp(-delta))/
			(1+exp(-delta*pt(x,df))*((1/u)-1)))/delta
		lines(x,qt(v,df))
		}
	Dat <- NULL
	Dat$x <- x
	Dat$y <- y
	deltas <- matrix(0,3,length(us))
	FrankModel <- function(x,delta,mu,sigma,df,tau){
        	z <- qt(-log(1-(1-exp(-delta))/
			(1+exp(-delta*pt(x,df))*((1/tau)-1)))/delta,df)
		mu + sigma*z
		}
	for(i in 1:length(us)){
		tau = us[i]
		fit <- nlrq(y~FrankModel(x,delta,mu,sigma,df=8,tau=tau),
			data=Dat,tau= tau, start=list(delta=5,
			mu = 0, sigma = 1),trace=TRUE)
		lines(x, predict(fit, newdata=x), lty=2, col="green")
		deltas[i,] <- coef(fit)
		}


###################################################
### code chunk number 21: lprq
###################################################
"lprq" <-
function(x, y, h, m=50 , tau=.5)
{
        xx <- seq(min(x),max(x),length=m)
        fv <- xx
        dv <- xx
        for(i in 1:length(xx)) {
                z <- x - xx[i]
                wx <- dnorm(z/h)
                r <- rq(y~z, weights=wx, tau=tau, ci=FALSE)
                fv[i] <- r$coef[1.]
                dv[i] <- r$coef[2.]
        }
        list(xx = xx, fv = fv, dv = dv)
}


###################################################
### code chunk number 22: mcycle1
###################################################
getOption("SweaveHooks")[["fig"]]()
library(MASS)
data(mcycle)
attach(mcycle)
plot(times,accel,xlab = "milliseconds", ylab = "acceleration")
hs <- c(1,2,3,4)
for(i in hs){
        h = hs[i]
        fit <- lprq(times,accel,h=h,tau=.5)
        lines(fit$xx,fit$fv,lty=i)
        }
legend(45,-70,c("h=1","h=2","h=3","h=4"),lty=1:length(hs))


###################################################
### code chunk number 23: regspline
###################################################
getOption("SweaveHooks")[["fig"]]()
library(splines)
plot(times,accel,xlab = "milliseconds", ylab = "acceleration",type="n")
points(times,accel,cex = .75)
X <- model.matrix(accel ~ bs(times, df=15))
for(tau in 1:3/4){
	fit <- rq(accel ~ bs(times, df=15), tau=tau, data=mcycle)
	accel.fit <- X %*% fit$coef
	lines(times,accel.fit)
	}


###################################################
### code chunk number 24: nprq
###################################################
data(Mammals)
attach(Mammals)


###################################################
### code chunk number 25: mammals
###################################################
getOption("SweaveHooks")[["fig"]]()
x <- log(weight)
y <- log(speed)
plot(x,y, xlab="Weight in log(Kg)", ylab="Speed in log(Km/hour)",type="n")
points(x[hoppers],y[hoppers],pch = "h", col="red")
points(x[specials],y[specials],pch = "s", col="blue")
others <- (!hoppers & !specials)
points(x[others],y[others], col="black",cex = .75)
fit <- rqss(y ~ qss(x, lambda = 1),tau = .9)
plot(fit, add = TRUE)


###################################################
### code chunk number 26: cobar
###################################################
getOption("SweaveHooks")[["fig"]]()
data(CobarOre)
fit <- rqss(z ~ qss(cbind(x,y), lambda = .01, ndum=100),data = CobarOre)
plot(fit, axes = FALSE, xlab = "", ylab = "")
rm(list=ls())


