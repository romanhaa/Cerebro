### R code from vignette source 'bootstrapping.Rnw'

###################################################
### code chunk number 1: bootstrapping.Rnw:11-34
###################################################
options(useFancyQuotes = FALSE)
digits <- 3
Colors <- c("#E495A5", "#39BEB1")
critical_values <- function(n, p = "0.95") {
  data("qDiptab", package = "diptest")
  if (n %in% rownames(qDiptab)) {
    return(qDiptab[as.character(n), p])
  }else {
    n.approx <- as.numeric(rownames(qDiptab)[which.min(abs(n - as.numeric(rownames(qDiptab))))])
    return(sqrt(n.approx)/sqrt(n) * qDiptab[as.character(n.approx), p])
  }
}
library("graphics")
library("flexmix")
combine <- function(x, sep, width) {
    cs <- cumsum(nchar(x))
      remaining <- if (any(cs[-1] > width)) combine(x[c(FALSE, cs[-1] > width)], sep, width)
      c(paste(x[c(TRUE, cs[-1] <= width)], collapse= sep), remaining)
}
prettyPrint <- function(x, sep = " ", linebreak = "\n\t", width = getOption("width")) {
  x <- strsplit(x, sep)[[1]]
  paste(combine(x, sep, width), collapse = paste(sep, linebreak, collapse = ""))
}


###################################################
### code chunk number 2: bootstrapping.Rnw:95-100
###################################################
cat(prettyPrint(gsub("boot_flexmix", "boot", 
                     prompt(flexmix:::boot_flexmix, 
                            filename = NA)$usage[[2]]), sep = ", ", 
                linebreak = paste("\n", paste(rep(" ", 2), collapse = ""), sep= ""),
                width = 70))


###################################################
### code chunk number 3: bootstrapping.Rnw:195-200
###################################################
library("flexmix")
Component_1 <- list(Model_1 = list(coef = c(1, -2), sigma = sqrt(0.1)))
Component_2 <- list(Model_1 = list(coef = c(2,  2), sigma = sqrt(0.1)))
ArtEx.mix <- FLXdist(y ~ x, k = rep(0.5, 2),
  components = list(Component_1, Component_2))


###################################################
### code chunk number 4: bootstrapping.Rnw:211-216
###################################################
ArtEx.data <- data.frame(x = rep(0:1, each = 100/2))
set.seed(123)
ArtEx.sim <- rflexmix(ArtEx.mix, newdata = ArtEx.data)
ArtEx.data$y <- ArtEx.sim$y[[1]]
ArtEx.data$class <- ArtEx.sim$class


###################################################
### code chunk number 5: bootstrapping.Rnw:225-230
###################################################
par(mar = c(5, 4, 2, 0) + 0.1)
plot(y ~ x, data = ArtEx.data, pch = with(ArtEx.data, 2*class + x))
pars <- list(matrix(c(1, -2, 2,  2), ncol = 2),
  matrix(c(1,  3, 2, -3), ncol = 2))
for (i in 1:2) apply(pars[[i]], 2, abline, col = Colors[i])


###################################################
### code chunk number 6: bootstrapping.Rnw:238-241
###################################################
set.seed(123)
ArtEx.fit <- stepFlexmix(y ~ x, data = ArtEx.data, k = 2, nrep = 5, 
  control = list(iter = 1000, tol = 1e-8, verbose = 0))


###################################################
### code chunk number 7: bootstrapping.Rnw:246-248
###################################################
summary(ArtEx.fit)
parameters(ArtEx.fit)


###################################################
### code chunk number 8: bootstrapping.Rnw:256-258
###################################################
ArtEx.refit <- refit(ArtEx.fit)
summary(ArtEx.refit)


###################################################
### code chunk number 9: bootstrapping.Rnw:280-284
###################################################
set.seed(123)
ArtEx.bs <- boot(ArtEx.fit, R = 15, sim = "parametric")
if ("boot-output.rda" %in% list.files()) load("boot-output.rda")
ArtEx.bs


###################################################
### code chunk number 10: bootstrapping.Rnw:300-301
###################################################
print(plot(ArtEx.bs, ordering = "coef.x", col = Colors))


###################################################
### code chunk number 11: bootstrapping.Rnw:318-331
###################################################
require("diptest")
parameters <- parameters(ArtEx.bs)
Ordering <- factor(as.vector(apply(matrix(parameters[,"coef.x"], 
  nrow = 2), 2, order)))
Comp1 <- parameters[Ordering == 1,]
Comp2 <- parameters[Ordering == 2,]
dip.values.art <- matrix(nrow = ncol(parameters), ncol = 3, 
  dimnames=list(colnames(parameters),
  c("Aggregated", "Comp 1", "Comp 2")))
dip.values.art[,"Aggregated"] <- apply(parameters, 2, dip)
dip.values.art[,"Comp 1"] <- apply(Comp1, 2, dip)
dip.values.art[,"Comp 2"] <- apply(Comp2, 2, dip)
dip.values.art


###################################################
### code chunk number 12: bootstrapping.Rnw:373-379
###################################################
data("seizure", package = "flexmix")
model <- FLXMRglm(family = "poisson", offset = log(seizure$Hours))
control <- list(iter = 1000, tol = 1e-10, verbose = 0)
set.seed(123)
seizMix <- stepFlexmix(Seizures ~ Treatment * log(Day), 
  data = seizure, k = 2, nrep = 5, model = model, control = control)


###################################################
### code chunk number 13: bootstrapping.Rnw:387-392
###################################################
par(mar = c(5, 4, 2, 0) + 0.1)
plot(Seizures/Hours~Day, data=seizure, pch = as.integer(seizure$Treatment))
abline(v = 27.5, lty = 2, col = "grey")
matplot(seizure$Day, fitted(seizMix)/seizure$Hours, type="l",
  add = TRUE, col = 1, lty = 1, lwd = 2)


###################################################
### code chunk number 14: bootstrapping.Rnw:414-418
###################################################
set.seed(123)
seizMix.bs <- boot(seizMix, R = 15, sim = "parametric")
if ("boot-output.rda" %in% list.files()) load("boot-output.rda")
print(plot(seizMix.bs, ordering = "coef.(Intercept)", col = Colors))


###################################################
### code chunk number 15: bootstrapping.Rnw:425-430
###################################################
parameters <- parameters(seizMix.bs)
Ordering <- factor(as.vector(apply(matrix(parameters[,"coef.(Intercept)"], 
  nrow = 2), 2, order)))
Comp1 <- parameters[Ordering == 1,]
Comp2 <- parameters[Ordering == 2,]


###################################################
### code chunk number 16: bootstrapping.Rnw:439-446
###################################################
dip.values.art <- matrix(nrow = ncol(parameters), ncol = 3, 
  dimnames = list(colnames(parameters), 
  c("Aggregated", "Comp 1", "Comp 2")))
dip.values.art[,"Aggregated"] <- apply(parameters, 2, dip)
dip.values.art[,"Comp 1"] <- apply(Comp1, 2, dip)
dip.values.art[,"Comp 2"] <- apply(Comp2, 2, dip)
dip.values.art


###################################################
### code chunk number 17: bootstrapping.Rnw:461-466 (eval = FALSE)
###################################################
## set.seed(123)
## ArtEx.bs <- boot(ArtEx.fit, R = 200, sim = "parametric")
## set.seed(123)
## seizMix.bs <- boot(seizMix, R = 200, sim = "parametric")
## save(ArtEx.bs, seizMix.bs, file = "boot-output.rda")


