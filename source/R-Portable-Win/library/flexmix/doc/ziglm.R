setClass("FLXMRziglm", contains = "FLXMRglm")

FLXMRziglm <- function(formula = . ~ .,
                       family = c("binomial", "poisson"), ...)
{
  family <- match.arg(family)
  new("FLXMRziglm", FLXMRglm(formula, family, ...),
      name = paste("FLXMRziglm", family, sep=":"))
}

setMethod("FLXgetModelmatrix", signature(model="FLXMRziglm"),
          function(model, data, formula, lhs=TRUE, ...)
{
  model <- callNextMethod(model, data, formula, lhs)
  if (attr(terms(model@fullformula), "intercept") == 0)
    stop("please include an intercept")
  model
})

setMethod("FLXremoveComponent", signature(model = "FLXMRziglm"),
          function(model, nok, ...)
{
  if (1 %in% nok) as(model, "FLXMRglm") else model
})

setMethod("FLXmstep", signature(model = "FLXMRziglm"),
          function(model, weights, components, ...)
{
  coef <- c(-Inf, rep(0, ncol(model@x)-1))
  names(coef) <- colnames(model@x)
  comp.1 <- with(list(coef = coef, df = 0, offset = NULL,
                 family = model@family), eval(model@defineComponent))
  c(list(comp.1),
    FLXmstep(as(model, "FLXMRglm"), weights[, -1, drop=FALSE],
             components[-1]))
})

