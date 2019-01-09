myConcomitant <-
  function(formula = ~ 1) {
    z <- new("FLXP",
             name = "myConcomitant",
             formula = formula)
    z@fit <- function(x, y, w, ...) {
      if (missing(w) || is.null(w)) w <- rep(1, length(x))
      f <- as.integer(factor(apply(x, 1, paste,
                                   collapse = "")))
      AVG <- apply(w*y, 2, tapply, f, mean)
      (AVG/rowSums(AVG))[f,,drop=FALSE]
    }
    z@refit <- function(x, y, w, ...) {
      if (missing(w) || is.null(w)) w <- rep(1, length(x))
      f <- as.integer(factor(apply(x, 1, paste,
                                   collapse = "")))
      AVG <- apply(w*y, 2, tapply, f, mean)
      (AVG/rowSums(AVG))
    }
    z
  }

