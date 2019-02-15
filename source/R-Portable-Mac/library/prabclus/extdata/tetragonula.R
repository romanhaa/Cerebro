tetragonula <- read.table("libsrc/prabclus/data/Heterotrigona_indoFO.txt",colClasses="character")
tetragonula.coord <- read.table("libsrc/prabclus/data/Franck04koord.txt")
save(tetragonula,tetragonula.coord,file="libsrc/prabclus/data/tetragonula.rda")
