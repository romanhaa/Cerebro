
## splitting
test_DataFrame_splitting <- function() {
  data(swiss)
  rn <- rownames(swiss)
  sw <- DataFrame(swiss, row.names=rn)
  swisssplit <- split(swiss, swiss$Education)

  ## split
  swsplit <- split(sw, sw[["Education"]])
  checkTrue(validObject(swsplit))
  checkIdentical(as.list(lapply(swsplit, as.data.frame)), swisssplit)
  checkTrue(validObject(split(DataFrame(IRanges(1:26, 1:26), LETTERS),
                              letters)))
}

