mymclust <- function (formula = .~., diagonal = TRUE) 
{    
    retval <- new("FLXMC", weighted = TRUE,
                  formula = formula, dist = "mvnorm",
                  name = "my model-based clustering")
    retval@defineComponent <- function(para) {
        logLik <- function(x, y) {
            mvtnorm::dmvnorm(y, mean = para$center,
                             sigma = para$cov, log = TRUE)
        }
        predict <- function(x) {
            matrix(para$center, nrow = nrow(x),
                   ncol = length(para$center), byrow = TRUE)
        }
        new("FLXcomponent", parameters =
            list(center = para$center, cov = para$cov),
            df = para$df, logLik = logLik, predict = predict)
    }
    retval@fit <- function(x, y, w, ...) {
        
        para <- cov.wt(y, wt = w)[c("center", "cov")]
        df <- (3 * ncol(y) + ncol(y)^2)/2
        
        if (diagonal) {
            para$cov <- diag(diag(para$cov))
            df <- 2 * ncol(y)
        }
        retval@defineComponent(c(para, df = df))
    }
    retval
}
