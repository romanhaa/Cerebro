### R code from vignette source 'mixture-regressions.Rnw'

###################################################
### code chunk number 1: mixture-regressions.Rnw:63-73
###################################################
options(width=60, prompt = "R> ", continue = "+  ", useFancyQuotes = FALSE)
library("graphics")
library("stats")
library("flexmix")
library("lattice")
ltheme <- canonical.theme("postscript", FALSE)
lattice.options(default.theme=ltheme)
data("NPreg", package = "flexmix")
data("dmft", package = "flexmix")
source("myConcomitant.R")


###################################################
### code chunk number 2: mixture-regressions.Rnw:500-503
###################################################
par(mfrow=c(1,2))
plot(yn~x, col=class, pch=class, data=NPreg)
plot(yp~x, col=class, pch=class, data=NPreg)


###################################################
### code chunk number 3: mixture-regressions.Rnw:510-517
###################################################
set.seed(1802)
library("flexmix")
data("NPreg", package = "flexmix")
Model_n <- FLXMRglm(yn ~ . + I(x^2))
Model_p <- FLXMRglm(yp ~ ., family = "poisson")
m1 <- flexmix(. ~ x, data = NPreg, k = 2, model = list(Model_n, Model_p),
  control = list(verbose = 10))


###################################################
### code chunk number 4: mixture-regressions.Rnw:558-559
###################################################
print(plot(m1))


###################################################
### code chunk number 5: mixture-regressions.Rnw:598-600
###################################################
m1.refit <- refit(m1)
summary(m1.refit, which = "model", model = 1)


###################################################
### code chunk number 6: mixture-regressions.Rnw:605-612
###################################################
print(plot(m1.refit, layout = c(1,3), bycluster = FALSE,
      main = expression(paste(yn *tilde(" ")* x + x^2))),
      split= c(1,1,2,1), more = TRUE)
print(plot(m1.refit, model = 2, 
           main = expression(paste(yp *tilde(" ")* x)), 
           layout = c(1,2), bycluster = FALSE), 
      split = c(2,1,2,1))


###################################################
### code chunk number 7: mixture-regressions.Rnw:643-648
###################################################
Model_n2 <- FLXMRglmfix(yn ~ . + 0, nested = list(k = c(1, 1), 
  formula = c(~ 1 + I(x^2), ~ 0)))
m2 <- flexmix(. ~ x, data = NPreg, cluster = posterior(m1), 
  model = list(Model_n2, Model_p))
m2


###################################################
### code chunk number 8: mixture-regressions.Rnw:653-654
###################################################
c(BIC(m1), BIC(m2))


###################################################
### code chunk number 9: mixture-regressions.Rnw:672-676
###################################################
data("betablocker", package = "flexmix")
betaGlm <- glm(cbind(Deaths, Total - Deaths) ~ Treatment, 
  family = "binomial", data = betablocker)
betaGlm


###################################################
### code chunk number 10: mixture-regressions.Rnw:693-696
###################################################
betaMixFix <- stepFlexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center,
  model = FLXMRglmfix(family = "binomial", fixed = ~ Treatment), 
  k = 2:4, nrep = 5, data = betablocker)


###################################################
### code chunk number 11: mixture-regressions.Rnw:705-706
###################################################
betaMixFix


###################################################
### code chunk number 12: mixture-regressions.Rnw:713-715
###################################################
betaMixFix_3 <- getModel(betaMixFix, which = "BIC")
betaMixFix_3 <- relabel(betaMixFix_3, "model", "Intercept")


###################################################
### code chunk number 13: mixture-regressions.Rnw:728-729
###################################################
parameters(betaMixFix_3)


###################################################
### code chunk number 14: mixture-regressions.Rnw:737-750
###################################################
library("grid")
betablocker$Center <- with(betablocker, factor(Center, levels = Center[order((Deaths/Total)[1:22])]))
clusters <- factor(clusters(betaMixFix_3), labels = paste("Cluster", 1:3))
print(dotplot(Deaths/Total ~ Center | clusters, groups  = Treatment, as.table = TRUE,
              data  = betablocker, xlab = "Center", layout = c(3, 1), 
              scales  = list(x = list(cex = 0.7, tck = c(1, 0))),
              key = simpleKey(levels(betablocker$Treatment), lines = TRUE, corner = c(1,0))))
betaMixFix.fitted <- fitted(betaMixFix_3)
for (i in 1:3) {
  seekViewport(trellis.vpname("panel", i, 1))
  grid.lines(unit(1:22, "native"), unit(betaMixFix.fitted[1:22, i], "native"), gp = gpar(lty = 1))
  grid.lines(unit(1:22, "native"), unit(betaMixFix.fitted[23:44, i], "native"), gp = gpar(lty = 2))
}


###################################################
### code chunk number 15: mixture-regressions.Rnw:769-775
###################################################
betaMix <- stepFlexmix(cbind(Deaths, Total - Deaths) ~ Treatment | Center,
  model = FLXMRglm(family = "binomial"), k = 3, nrep = 5, 
  data = betablocker)
betaMix <- relabel(betaMix, "model", "Treatment")
parameters(betaMix)
c(BIC(betaMixFix_3), BIC(betaMix))


###################################################
### code chunk number 16: mixture-regressions.Rnw:795-796
###################################################
print(plot(betaMixFix_3, nint = 10, mark = 1, col = "grey", layout = c(3, 1)))


###################################################
### code chunk number 17: mixture-regressions.Rnw:805-806
###################################################
print(plot(betaMixFix_3, nint = 10, mark = 2, col = "grey", layout = c(3, 1)))


###################################################
### code chunk number 18: mixture-regressions.Rnw:820-821
###################################################
table(clusters(betaMix))


###################################################
### code chunk number 19: mixture-regressions.Rnw:826-828
###################################################
predict(betaMix, 
  newdata = data.frame(Treatment = c("Control", "Treated")))


###################################################
### code chunk number 20: mixture-regressions.Rnw:834-836
###################################################
betablocker[c(1, 23), ]
fitted(betaMix)[c(1, 23), ]


###################################################
### code chunk number 21: mixture-regressions.Rnw:846-847
###################################################
summary(refit(betaMix))


###################################################
### code chunk number 22: mixture-regressions.Rnw:858-865
###################################################
ModelNested <- FLXMRglmfix(family = "binomial", nested = list(k = c(2, 1),
  formula = c(~ Treatment, ~ 0)))
betaMixNested <- flexmix(cbind(Deaths, Total - Deaths) ~ 1 | Center,
  model = ModelNested, k = 3, data = betablocker, 
  cluster = posterior(betaMix))
parameters(betaMixNested)
c(BIC(betaMix), BIC(betaMixNested), BIC(betaMixFix_3))


###################################################
### code chunk number 23: mixture-regressions.Rnw:876-877
###################################################
data("bioChemists", package = "flexmix")


###################################################
### code chunk number 24: mixture-regressions.Rnw:908-912
###################################################
data("bioChemists", package = "flexmix")
Model1 <- FLXMRglm(family = "poisson")
ff_1 <- stepFlexmix(art ~ ., data = bioChemists, k = 1:3, model = Model1)
ff_1 <- getModel(ff_1, "BIC")


###################################################
### code chunk number 25: mixture-regressions.Rnw:929-931
###################################################
print(plot(refit(ff_1), bycluster = FALSE, 
           scales = list(x = list(relation = "free"))))


###################################################
### code chunk number 26: mixture-regressions.Rnw:938-942
###################################################
Model2 <- FLXMRglmfix(family = "poisson", fixed = ~ kid5 + mar + ment)
ff_2 <- flexmix(art ~ fem + phd, data = bioChemists, 
  cluster = posterior(ff_1), model = Model2)
c(BIC(ff_1), BIC(ff_2))


###################################################
### code chunk number 27: mixture-regressions.Rnw:950-951
###################################################
summary(refit(ff_2))


###################################################
### code chunk number 28: mixture-regressions.Rnw:958-962
###################################################
Model3 <- FLXMRglmfix(family = "poisson", fixed = ~ kid5 + mar + ment)
ff_3 <- flexmix(art ~ fem, data = bioChemists, cluster = posterior(ff_2),
  model = Model3)
c(BIC(ff_2), BIC(ff_3))


###################################################
### code chunk number 29: mixture-regressions.Rnw:970-971
###################################################
print(plot(refit(ff_3), bycluster = FALSE, scales = list(x = list(relation = "free"))))


###################################################
### code chunk number 30: mixture-regressions.Rnw:981-987
###################################################
Model4 <- FLXMRglmfix(family = "poisson", fixed = ~ kid5 + mar + ment)
ff_4 <- flexmix(art ~ 1, data = bioChemists, cluster = posterior(ff_2),
  concomitant = FLXPmultinom(~ fem), model = Model4)
parameters(ff_4)
summary(refit(ff_4), which = "concomitant")
BIC(ff_4)


###################################################
### code chunk number 31: mixture-regressions.Rnw:996-1000
###################################################
Model5 <- FLXMRglmfix(family = "poisson", fixed = ~ kid5 + ment + fem)
ff_5 <- flexmix(art ~ 1, data = bioChemists, cluster = posterior(ff_2),
  model = Model5)
BIC(ff_5)


###################################################
### code chunk number 32: mixture-regressions.Rnw:1006-1013
###################################################
pp <- predict(ff_5, newdata = data.frame(kid5 = 0, 
  mar = factor("Married", levels = c("Single", "Married")),
  fem = c("Men", "Women"),  ment = mean(bioChemists$ment)))
matplot(0:12, sapply(unlist(pp), function(x) dpois(0:12, x)), 
        type = "b", lty = 1, xlab = "Number of articles", ylab = "Probability")
legend("topright", paste("Comp.", rep(1:2, each = 2), ":",
  c("Men", "Women")), lty = 1, col = 1:4, pch = paste(1:4), bty = "n")


###################################################
### code chunk number 33: mixture-regressions.Rnw:1362-1367
###################################################
data("dmft", package = "flexmix")
Model <- FLXMRziglm(family = "poisson")
Fitted <- flexmix(End ~ log(Begin + 0.5) + Gender + Ethnic + Treatment, 
  model = Model, k = 2 , data = dmft, control = list(minprior = 0.01))
summary(refit(Fitted))


###################################################
### code chunk number 34: refit (eval = FALSE)
###################################################
## print(plot(refit(Fitted), components = 2, box.ratio = 3))


###################################################
### code chunk number 35: mixture-regressions.Rnw:1396-1397
###################################################
print(plot(refit(Fitted), components = 2, box.ratio = 3))


###################################################
### code chunk number 36: mixture-regressions.Rnw:1442-1449
###################################################
Concomitant <- FLXPmultinom(~ yb)
MyConcomitant <- myConcomitant(~ yb)

m2 <- flexmix(. ~ x, data = NPreg, k = 2, model = list(Model_n, Model_p), 
  concomitant = Concomitant)
m3 <- flexmix(. ~ x, data = NPreg, k = 2, model = list(Model_n, Model_p), 
  cluster = posterior(m2), concomitant = MyConcomitant)


###################################################
### code chunk number 37: mixture-regressions.Rnw:1451-1453
###################################################
summary(m2)
summary(m3)


###################################################
### code chunk number 38: mixture-regressions.Rnw:1458-1462
###################################################
determinePrior <- function(object) {
  object@concomitant@fit(object@concomitant@x, 
    posterior(object))[!duplicated(object@concomitant@x), ]
}


###################################################
### code chunk number 39: mixture-regressions.Rnw:1465-1467
###################################################
determinePrior(m2)
determinePrior(m3)


###################################################
### code chunk number 40: mixture-regressions.Rnw:1509-1513
###################################################
SI <- sessionInfo()
pkgs <- paste(sapply(c(SI$otherPkgs, SI$loadedOnly), function(x) 
                     paste("\\\\pkg{", x$Package, "} ", 
                           x$Version, sep = "")), collapse = ", ")


