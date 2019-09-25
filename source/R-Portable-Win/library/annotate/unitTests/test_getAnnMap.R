## Set up 
require(org.Hs.eg.db)
require(annotate)
require(RUnit)

## For testing 
test_getAnnMap <- function(){
  ## test for a map that exist
  map <- getAnnMap("CHRLOC","org.Hs.eg.db")
  checkTrue( class(map) == "AnnDbMap" )
  ## and test for a map that does not (but which is available via select)
  map2 <- getAnnMap("ONTOLOGY","org.Hs.eg.db")
  checkTrue( class(map2) == "FlatBimap" )
}
