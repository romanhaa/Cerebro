## unit test to just check that mapIds is behaving itself
require(org.Hs.eg.db)

## constants
k <- head(keys(org.Hs.eg.db, 'ENTREZID'))

# trace("mapIds", tracer=browser, signature ="AnnotationDb")

### The actual tests.

## Default is currently 'first'
test_mapIds_default <- function(){
    res <- mapIds(org.Hs.eg.db, keys=k, column='ALIAS', keytype='ENTREZID')
    checkTrue(length(res) == length(k))
    checkTrue(res[[1]] == "A1B")
    checkTrue(class(res)=='character')
}

## test other return types.

## "list"
test_mapIds_CharacterList <- function(){
    res <- mapIds(org.Hs.eg.db, keys=k, column='ALIAS', keytype='ENTREZID', 
                  multiVals="list")
    checkTrue(length(res) == length(k))
    checkTrue(res[[1]][1] == "A1B")
    checkTrue(class(res)=='list')
}

## "CharacterList"
test_mapIds_CharacterList <- function(){
    res <- mapIds(org.Hs.eg.db, keys=k, column='ALIAS', keytype='ENTREZID', 
       multiVals="CharacterList")
    checkTrue(length(res) == length(k))
    checkTrue(res[[1]][1] == "A1B")
    checkTrue(class(res)=='SimpleCharacterList')
}

## "NAMultiples"
test_mapIds_NAMultiples <- function(){
    res <- mapIds(org.Hs.eg.db, keys=k, column='PFAM', keytype='ENTREZID', 
       multiVals="asNA")
    checkTrue(length(res) == length(k))
    checkTrue(res[['10']] == "PF00797")
    checkTrue(class(res)=='character')
}

## "filterMultiples"
test_mapIds_filterMultiples <- function(){
    res <- mapIds(org.Hs.eg.db, keys=k, column='PFAM', keytype='ENTREZID', 
       multiVals="filter")
    checkTrue(length(res) < length(k)) ## multi matchers means should be shorter
    checkTrue(res[['10']] == "PF00797")
    checkTrue(res[['1']] == "PF13895")
    checkTrue(class(res)=='character')
}


## CUSTOM function
test_mapIds_FUN <- function(){
    last <- function(x){
    x[[length(x)]]
    }
    res <- mapIds(org.Hs.eg.db, keys=k, column='ALIAS', keytype='ENTREZID', 
       multiVals=last)
    checkTrue(length(res) == length(k))
    checkTrue(res[[1]] == "A1BG")
    checkTrue(class(res)=='character')
}
