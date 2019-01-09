veronica <- as.matrix(read.table("libsrc/prabclus/data/MartinezOrtega04AFLP.txt"))
veronica.coord <- as.matrix(read.table("libsrc/prabclus/data/MartinezKoord.txt"))
save(veronica,veronica.coord,file="libsrc/prabclus/data/veronica.rda")

