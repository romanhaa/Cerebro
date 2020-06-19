#some simple test classes to see if things are working at all

library(methods)

setClass("foo", representation(a="numeric", b="numeric"))

setClass("bar", representation(c="numeric", b="character"))

setClass("baz", "foo", representation(d="list", e="logical"))


#load up Biobase so we can test

library(Biobase)

f1 <- new("foo", a=10,b=15)
l1 <- list(a=f1)

nc1 <- new("container", x=l1, content="foo", locked=FALSE)

b1 <- new("bar", b="AAA", c=10)

#this should fail...

 nc1[[1]] <- b1

 nc1[[1]]

 v<-nc1[1]
