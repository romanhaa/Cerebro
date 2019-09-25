test_FilterRules_construct <- function() {

  ## as a simple character vector
  filts <- c("peaks", "promoters")
  parsedFilts <- list(peaks = expression(peaks),
                      promoters = expression(promoters))

  filters <- FilterRules()
  checkTrue(validObject(filters))
  checkIdentical(as.list(filters), list())
  
  filters <- FilterRules(filts)
  checkTrue(validObject(filters))
  checkIdentical(as.list(filters), parsedFilts)
  checkIdentical(active(filters), structure(rep(TRUE, 2), names=filts))

  ## with functions and expressions
  filts <- c(parsedFilts, list(find_eboxes = function(rd) rep(FALSE, nrow(rd))))
  filters <- FilterRules(filts, active = FALSE)
  checkTrue(validObject(filters))
  filts$find_eboxes <- new("FilterClosure", filts$find_eboxes)
  checkIdentical(as.list(filters), filts)
  checkIdentical(active(filters), structure(rep(FALSE, 3), names=names(filts)))

  ## direct, quoted args (character literal parsed)
  filters <- FilterRules(under_peaks = peaks, in_promoters = "promoters")
  filts <- list(under_peaks = expression(peaks),
                in_promoters = expression(promoters))
  checkTrue(validObject(filters))
  checkIdentical(as.list(filters), filts)
  ## mix them up
  filters <- FilterRules(filts, diffexp = de)
  checkTrue(validObject(filters))
  checkIdentical(as.list(filters), c(filts, list(diffexp = expression(de))))
  filts <- as.list(filters)
  
  checkException(FilterRules(c(filts, 1)), silent = TRUE)
  checkException(FilterRules(filts, active = filts), silent = TRUE)
  checkException(FilterRules(list(find_eboxes = function() NULL)),
                 silent = TRUE)
}

test_FilterRules_append <- function() {
  filts <- c("peaks", "promoters")
  filts2 <- c("introns", "exons")

  filters <- FilterRules(filts)
  filters2 <- FilterRules(filts2, active=FALSE)

  both <- append(filters, filters2)
  checkTrue(validObject(both))
  bothFilts <- structure(list(quote(peaks), quote(promoters),
                              quote(introns), quote(exons)),
                         names = c(filts, filts2))
  checkIdentical(unlist(as.list(both)), bothFilts)
  bothActive <- structure(c(TRUE, TRUE, FALSE, FALSE), names = names(bothFilts))
  checkIdentical(active(both), bothActive)
  both <- c(filters, filters2)
  checkTrue(validObject(both))
  checkIdentical(unlist(as.list(both)), bothFilts)
  checkIdentical(active(both), bothActive)

  filters[["cons"]] <- "cons"
  filts <- list(peaks = quote(peaks), promoters = quote(promoters))
  filts <- c(filts, cons = quote(cons))
  checkIdentical(unlist(as.list(filters)), filts)
  filters[["cons"]] <- quote(cons)
  checkIdentical(unlist(as.list(filters)), filts)
  filters[["cons"]] <- expression(cons)
  checkIdentical(unlist(as.list(filters)), filts)
  fun <- function(rd) rep(FALSE, nrow(rd))
  filters[[4]] <- fun
  filts <- c(filts, X = new("FilterClosure", fun))
  checkIdentical(unlist(as.list(filters)), filts)
  
  checkException(filters[[]] <- "threeprime", silent = TRUE)
  checkException(filters[[1]] <- 2, silent = TRUE)
  checkException(filters[[1]] <- list(quote(foo), quote(bar)), silent = TRUE)
}

test_FilterRules_subset <- function() {
  filts <- c("peaks", "promoters", "introns")
  filters <- FilterRules(filts)

  checkIdentical(sapply(unlist(filters[1:2]), deparse),
                 structure(filts[1:2], names = filts[1:2]))
  checkIdentical(sapply(unlist(filters[]),deparse),
                 structure(filts, names = filts))
}

test_FilterRules_active <- function() {
  filts <- c("peaks", "promoters", "introns")
  filters <- FilterRules(filts)

  ## set the active state directly
  
  active(filters) <- FALSE
  checkIdentical(active(filters), structure(rep(FALSE, 3), names = filts))
  active(filters) <- TRUE
  checkIdentical(active(filters), structure(rep(TRUE, 3), names = filts))
  active(filters) <- c(FALSE, FALSE, TRUE)
  checkIdentical(active(filters),
                 structure(c(FALSE, FALSE, TRUE), names = filts))
  active(filters)["promoters"] <- TRUE
  checkIdentical(active(filters),
                 structure(c(FALSE, TRUE, TRUE), names = filts))
  checkException(active(filters) <- rep(FALSE, 2), silent = TRUE)
  checkException(active(filters) <- rep(FALSE, 5), silent = TRUE)
  checkException(active(filters)["introns"] <- NA, silent = TRUE)
  
  ## toggle the active state by name or index
  
  active(filters) <- c(NA, 2) # NA's are dropped
  checkIdentical(active(filters),
                 structure(c(FALSE, TRUE, FALSE), names = filts))
  active(filters) <- c("peaks", NA) 
  checkIdentical(active(filters),
                 structure(c(TRUE, FALSE, FALSE), names = filts))
  checkException(active(filters) <- "foo", silent = TRUE)
  checkException(active(filters) <- 15, silent = TRUE)
}

test_FilterRules_annotation <- function() {
  filts <- c("peaks", "promoters")
  filters <- FilterRules(filts)
  mcols(filters) <- DataFrame(a = 1:2)
  checkIdentical(mcols(filters)[,1], 1:2)
  checkIdentical(mcols(filters[2:1])[,1], 2:1)
  checkIdentical(mcols(c(filters,filters))[,1], rep(1:2,2))
  checkIdentical(mcols(append(filters,filters))[,1], rep(1:2,2))
}
