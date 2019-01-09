nb <- list()
for (i in 1:34)
  nb <- c(nb,list(scan(file="libsrc/prabclus/data/nb.txt",
                   skip=i-1,nlines=1, quiet=TRUE)))
remove(i)
save(nb,file="libsrc/prabclus/data/nb.rda")
