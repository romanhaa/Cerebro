### R code from vignette source 'regression-examples.Rnw'

###################################################
### code chunk number 1: regression-examples.Rnw:11-14
###################################################
library("stats")
library("graphics")
library("flexmix")


###################################################
### code chunk number 2: start
###################################################
options(width=70, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
set.seed(1802)
library("lattice")
ltheme <- canonical.theme("postscript", FALSE)
lattice.options(default.theme=ltheme)


###################################################
### code chunk number 3: NregFix
###################################################
set.seed(2807)
library("flexmix")
data("NregFix", package = "flexmix")
Model <- FLXMRglm(~ x2 + x1)
fittedModel <- stepFlexmix(y ~ 1, model = Model, nrep = 3, k = 3,
  data = NregFix, concomitant = FLXPmultinom(~ w))             
fittedModel <- relabel(fittedModel, "model", "x1")
summary(refit(fittedModel))


###################################################
### code chunk number 4: diffModel
###################################################
Model2 <- FLXMRglmfix(fixed = ~ x2, nested = list(k = c(1, 2), 
  formula = c(~ 0, ~ x1)), varFix = TRUE)
fittedModel2 <- flexmix(y ~ 1, model = Model2, 
  cluster = posterior(fittedModel), data = NregFix, 
  concomitant = FLXPmultinom(~ w))             
BIC(fittedModel)
BIC(fittedModel2)


###################################################
### code chunk number 5: artificial-example
###################################################
plotNregFix <- NregFix
plotNregFix$w <- factor(NregFix$w, levels = 0:1, labels = paste("w =", 0:1))
plotNregFix$x2 <- factor(NregFix$x2, levels = 0:1,
  labels = paste("x2 =", 0:1))
plotNregFix$class <- factor(NregFix$class, levels = 1:3, labels = paste("Class", 1:3))
print(xyplot(y ~ x1 | x2*w, groups = class, data = plotNregFix, cex = 0.6,
  auto.key = list(space="right"), layout = c(2,2)))


###################################################
### code chunk number 6: refit
###################################################
summary(refit(fittedModel2))


###################################################
### code chunk number 7: beta-glm
###################################################
data("betablocker", package = "flexmix")
betaGlm <- glm(cbind(Deaths, Total - Deaths) ~ Treatment, 
  family = "binomial", data = betablocker)
betaGlm


###################################################
### code chunk number 8: beta-fix
###################################################
betaMixFix <- stepFlexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center,
  model = FLXMRglmfix(family = "binomial", fixed = ~ Treatment), 
  k = 2:4, nrep = 3, data = betablocker)
betaMixFix


###################################################
### code chunk number 9: beta-fig
###################################################
library("grid")
betaMixFix_3 <- getModel(betaMixFix, "3")
betaMixFix_3 <- relabel(betaMixFix_3, "model", "Intercept")
betablocker$Center <- with(betablocker, factor(Center, levels = Center[order((Deaths/Total)[1:22])]))
clusters <- factor(clusters(betaMixFix_3), labels = paste("Cluster", 1:3))
print(dotplot(Deaths/Total ~ Center | clusters, groups  = Treatment, as.table = TRUE,
  data  = betablocker, xlab = "Center", layout = c(3, 1), scales = list(x = list(draw = FALSE)),
  key = simpleKey(levels(betablocker$Treatment), lines = TRUE, corner = c(1,0))))
betaMixFix.fitted <- fitted(betaMixFix_3)
for (i in 1:3) {
  seekViewport(trellis.vpname("panel", i, 1))
  grid.lines(unit(1:22, "native"), unit(betaMixFix.fitted[1:22, i], "native"), gp = gpar(lty = 1))
  grid.lines(unit(1:22, "native"), unit(betaMixFix.fitted[23:44, i], "native"), gp = gpar(lty = 2))
}


###################################################
### code chunk number 10: beta-full
###################################################
betaMix <- stepFlexmix(cbind(Deaths, Total - Deaths) ~ Treatment | Center,
  model = FLXMRglm(family = "binomial"), k = 3, nrep = 3, 
  data = betablocker)
summary(betaMix)


###################################################
### code chunk number 11: default-plot
###################################################
print(plot(betaMixFix_3, mark = 1, col = "grey", markcol = 1))


###################################################
### code chunk number 12: parameters
###################################################
parameters(betaMix)


###################################################
### code chunk number 13: cluster
###################################################
table(clusters(betaMix))


###################################################
### code chunk number 14: predict
###################################################
predict(betaMix, 
  newdata = data.frame(Treatment = c("Control", "Treated")))


###################################################
### code chunk number 15: fitted
###################################################
fitted(betaMix)[c(1, 23), ]


###################################################
### code chunk number 16: refit
###################################################
summary(refit(getModel(betaMixFix, "3")))


###################################################
### code chunk number 17: mehta-fix
###################################################
data("Mehta", package = "flexmix")
mehtaMix <- stepFlexmix(cbind(Response, Total - Response)~ 1 | Site, 
  model = FLXMRglmfix(family = "binomial", fixed = ~ Drug), 
  control = list(minprior = 0.04), nrep = 3, k = 3, data = Mehta)
summary(mehtaMix)


###################################################
### code chunk number 18: mehta-fig
###################################################
Mehta$Site <- with(Mehta, factor(Site, levels = Site[order((Response/Total)[1:22])]))
clusters <- factor(clusters(mehtaMix), labels = paste("Cluster", 1:3)) 
print(dotplot(Response/Total ~ Site | clusters, groups  = Drug,  layout = c(3,1),
              data  = Mehta, xlab = "Site", scales = list(x = list(draw = FALSE)),
              key = simpleKey(levels(Mehta$Drug), lines = TRUE, corner = c(1,0))))
mehtaMix.fitted <- fitted(mehtaMix)
for (i in 1:3) {
  seekViewport(trellis.vpname("panel", i, 1))
  sapply(1:nlevels(Mehta$Drug), function(j) 
  grid.lines(unit(1:22, "native"), unit(mehtaMix.fitted[Mehta$Drug == levels(Mehta$Drug)[j], i], "native"), gp = gpar(lty = j)))
}


###################################################
### code chunk number 19: mehta-full
###################################################
mehtaMix <- stepFlexmix(cbind(Response, Total - Response) ~ Drug | Site, 
  model = FLXMRglm(family = "binomial"), k = 3, data = Mehta, nrep = 3,
  control = list(minprior = 0.04))
summary(mehtaMix)


###################################################
### code chunk number 20: mehta-sub
###################################################
Mehta.sub <- subset(Mehta, Site != 15)
mehtaMix <- stepFlexmix(cbind(Response, Total - Response) ~ 1 | Site, 
  model = FLXMRglmfix(family = "binomial", fixed = ~ Drug),
  data = Mehta.sub, k = 2, nrep = 3)
summary(mehtaMix)


###################################################
### code chunk number 21: tribolium
###################################################
data("tribolium", package = "flexmix")
TribMix <- stepFlexmix(cbind(Remaining, Total - Remaining) ~ 1, 
  k = 2:3, model = FLXMRglmfix(fixed = ~ Species, family = "binomial"),
  concomitant = FLXPmultinom(~ Replicate), data = tribolium)


###################################################
### code chunk number 22: wang-model
###################################################
modelWang <- FLXMRglmfix(fixed = ~ I(Species == "Confusum"), 
  family = "binomial") 
concomitantWang <- FLXPmultinom(~ I(Replicate == 3))
TribMixWang <- stepFlexmix(cbind(Remaining, Total - Remaining) ~ 1,
  data = tribolium, model = modelWang, concomitant = concomitantWang, 
  k = 2) 
summary(refit(TribMixWang))


###################################################
### code chunk number 23: tribolium
###################################################
clusters <- factor(clusters(TribMixWang), labels = paste("Cluster", 1:TribMixWang@k))
print(dotplot(Remaining/Total ~ factor(Replicate) | clusters, groups  = Species, 
              data  = tribolium[rep(1:9, each = 3) + c(0:2)*9,], xlab = "Replicate",
              auto.key = list(corner = c(1,0))))


###################################################
### code chunk number 24: subset
###################################################
TribMixWangSub <- stepFlexmix(cbind(Remaining, Total - Remaining) ~ 1, 
  k = 2, data = tribolium[-7,], model = modelWang,
  concomitant = concomitantWang)


###################################################
### code chunk number 25: trypanosome
###################################################
data("trypanosome", package = "flexmix")
TrypMix <- stepFlexmix(cbind(Dead, 1-Dead) ~ 1, k = 2, nrep = 3, 
  data = trypanosome,  model = FLXMRglmfix(family = "binomial", 
  fixed = ~ log(Dose)))
summary(refit(TrypMix))


###################################################
### code chunk number 26: trypanosome
###################################################
tab <- with(trypanosome, table(Dead, Dose))
Tryp.dat <- data.frame(Dead = tab["1",], Alive = tab["0",], 
                       Dose = as.numeric(colnames(tab)))
plot(Dead/(Dead+Alive) ~ Dose, data = Tryp.dat)
Tryp.pred <- predict(glm(cbind(Dead, 1-Dead) ~ log(Dose), family = "binomial", data = trypanosome), newdata=Tryp.dat, type = "response")
TrypMix.pred <- predict(TrypMix, newdata = Tryp.dat, aggregate = TRUE)[[1]]
lines(Tryp.dat$Dose, Tryp.pred, lty = 2)
lines(Tryp.dat$Dose, TrypMix.pred, lty = 3)
legend(4.7, 1, c("GLM", "Mixture model"), lty=c(2, 3), xjust=0, yjust=1)


###################################################
### code chunk number 27: fabric-fix
###################################################
data("fabricfault", package = "flexmix")
fabricMix <- stepFlexmix(Faults ~ 1, model = FLXMRglmfix(family="poisson", 
  fixed = ~ log(Length)), data = fabricfault, k = 2, nrep = 3)
summary(fabricMix)
summary(refit(fabricMix))
Lnew <- seq(0, 1000, by = 50)
fabricMix.pred <- predict(fabricMix, newdata = data.frame(Length = Lnew))


###################################################
### code chunk number 28: fabric-fix-nested
###################################################
fabricMix2 <- flexmix(Faults ~ 0, data = fabricfault, 
  cluster = posterior(fabricMix), 
  model = FLXMRglmfix(family = "poisson", fixed = ~ log(Length),
  nested = list(k=c(1,1), formula=list(~0,~1))))
summary(refit(fabricMix2))
fabricMix2.pred <- predict(fabricMix2, 
  newdata = data.frame(Length = Lnew))


###################################################
### code chunk number 29: fabric-fig
###################################################
plot(Faults ~ Length, data = fabricfault)
sapply(fabricMix.pred, function(y) lines(Lnew, y, lty = 1))
sapply(fabricMix2.pred, function(y) lines(Lnew, y, lty = 2))
legend(190, 25, paste("Model", 1:2), lty=c(1, 2), xjust=0, yjust=1)


###################################################
### code chunk number 30: patent
###################################################
data("patent", package = "flexmix")
ModelPat <- FLXMRglm(family = "poisson")
FittedPat <- stepFlexmix(Patents ~ lgRD, k = 3, nrep = 3, 
  model = ModelPat, data = patent, concomitant = FLXPmultinom(~ RDS))             
summary(FittedPat)


###################################################
### code chunk number 31: patent-fixed
###################################################
ModelFixed <- FLXMRglmfix(family = "poisson", fixed = ~ lgRD)
FittedPatFixed <- flexmix(Patents ~ 1, model = ModelFixed, 
  cluster = posterior(FittedPat), concomitant = FLXPmultinom(~ RDS), 
  data = patent)             
summary(FittedPatFixed)


###################################################
### code chunk number 32: Poisson
###################################################
lgRDv <- seq(-3, 5, by = 0.05)
newdata <- data.frame(lgRD = lgRDv)
plotData <- function(fitted) {
  with(patent, data.frame(Patents = c(Patents, 
                            unlist(predict(fitted, newdata = newdata))),
                          lgRD = c(lgRD, rep(lgRDv, 3)),
                          class = c(clusters(fitted), rep(1:3, each = nrow(newdata))),
                          type = rep(c("data", "fit"), c(nrow(patent), nrow(newdata)*3))))
}
plotPatents <- cbind(plotData(FittedPat), which = "Wang et al.")
plotPatentsFixed <- cbind(plotData(FittedPatFixed), which = "Fixed effects")
plotP <- rbind(plotPatents, plotPatentsFixed)
rds <- seq(0, 3, by = 0.02)
x <- model.matrix(FittedPat@concomitant@formula, data = data.frame(RDS = rds))
plotConc <- function(fitted) {
  E <- exp(x%*%fitted@concomitant@coef)
  data.frame(Probability = as.vector(E/rowSums(E)),
                         class = rep(1:3, each = nrow(x)),
                         RDS = rep(rds, 3))
}
plotConc1 <- cbind(plotConc(FittedPat), which = "Wang et al.")
plotConc2 <- cbind(plotConc(FittedPatFixed), which = "Fixed effects")
plotC <- rbind(plotConc1, plotConc2)
print(xyplot(Patents ~ lgRD | which, data = plotP, groups=class, xlab = "log(R&D)",
             panel = "panel.superpose", type = plotP$type,
             panel.groups = function(x, y, type = "p", subscripts, ...)
             {
               ind <- plotP$type[subscripts] == "data"
               panel.xyplot(x[ind], y[ind], ...)
               panel.xyplot(x[!ind], y[!ind], type = "l", ...)
             },
             scales = list(alternating=FALSE), layout=c(1,2), as.table=TRUE),
      more=TRUE, position=c(0,0,0.6, 1))
print(xyplot(Probability ~ RDS | which, groups = class, data = plotC, type = "l",
             scales = list(alternating=FALSE), layout=c(1,2), as.table=TRUE),
      position=c(0.6, 0.01, 1, 0.99))


###################################################
### code chunk number 33: seizure
###################################################
data("seizure", package = "flexmix")
seizMix <- stepFlexmix(Seizures ~ Treatment * log(Day), data = seizure, 
  k = 2, nrep = 3, model = FLXMRglm(family = "poisson", 
  offset = log(seizure$Hours)))
summary(seizMix)
summary(refit(seizMix))


###################################################
### code chunk number 34: seizure
###################################################
seizMix2 <- flexmix(Seizures ~ Treatment * log(Day/27), 
  data = seizure, cluster = posterior(seizMix),
  model = FLXMRglm(family = "poisson", offset = log(seizure$Hours)))
summary(seizMix2)
summary(refit(seizMix2))


###################################################
### code chunk number 35: seizure
###################################################
seizMix3 <- flexmix(Seizures ~ log(Day/27)/Treatment, data = seizure, 
  cluster = posterior(seizMix), model = FLXMRglm(family = "poisson", 
  offset = log(seizure$Hours)))
summary(seizMix3)
summary(refit(seizMix3))


###################################################
### code chunk number 36: seizure
###################################################
plot(Seizures/Hours~Day, pch = c(1,3)[as.integer(Treatment)], data=seizure)
abline(v=27.5, lty=2, col="grey")
legend(140, 9, c("Baseline", "Treatment"), pch=c(1, 3), xjust=1, yjust=1)
matplot(seizure$Day, fitted(seizMix)/seizure$Hours, type="l", add=TRUE, lty = 1, col = 1)
matplot(seizure$Day, fitted(seizMix3)/seizure$Hours, type="l", add=TRUE, lty = 3, col = 1)
legend(140, 7, paste("Model", c(1,3)), lty=c(1, 3), xjust=1, yjust=1)


###################################################
### code chunk number 37: salmonella
###################################################
data("salmonellaTA98", package = "flexmix")
salmonMix <- stepFlexmix(y ~ 1, data = salmonellaTA98, k = 2, nrep = 3,
  model = FLXMRglmfix(family = "poisson", fixed = ~ x + log(x + 10)))


###################################################
### code chunk number 38: salmonella
###################################################
salmonMix.pr <- predict(salmonMix, newdata=salmonellaTA98)
plot(y~x, data=salmonellaTA98, 
     pch=as.character(clusters(salmonMix)), 
     xlab="Dose of quinoline", ylab="Number of revertant colonies of salmonella",
     ylim=range(c(salmonellaTA98$y, unlist(salmonMix.pr))))
for (i in 1:2) lines(salmonellaTA98$x, salmonMix.pr[[i]], lty=i)


###################################################
### code chunk number 39: regression-examples.Rnw:923-927
###################################################
SI <- sessionInfo()
pkgs <- paste(sapply(c(SI$otherPkgs, SI$loadedOnly), function(x) 
                     paste("\\\\pkg{", x$Package, "} ", 
                           x$Version, sep = "")), collapse = ", ")


