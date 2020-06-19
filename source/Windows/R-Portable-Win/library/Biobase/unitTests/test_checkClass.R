testCheckData_S3ExtendedDataFrame<- function() {
    df <- data.frame(foo=1:2)
    class(df) = c("bar", class(df))
    Biobase:::checkClass(df, "data.frame")
    checkTrue(TRUE)
}
