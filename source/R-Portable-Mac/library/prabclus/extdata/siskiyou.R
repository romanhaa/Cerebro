siskiyou <- as.matrix(read.table("libsrc/prabclus/data/LeiMik1.txt"))
siskiyou.nb <- list()
for (i in 1:6)
  siskiyou.nb <- c(siskiyou.nb,list(scan(file="libsrc/prabclus/data/LeiMik1NB.txt",
                   skip=i-1,nlines=1, quiet=TRUE)))
remove(i)
siskiyou.groups <- scan("libsrc/prabclus/data/LeiMik1G.txt")
save(siskiyou,siskiyou.nb,siskiyou.groups,file="libsrc/prabclus/data/siskiyou.rda")
