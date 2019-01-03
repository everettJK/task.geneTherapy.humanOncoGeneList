library(stringr)
library(tidyr)
library(dplyr)
options(stringsAsFactors = FALSE)

Sys.setlocale('LC_ALL','C') 
atlas <- unique(toupper(unlist(lapply(LETTERS, function(x){
  p <- readLines(paste0('http://atlasgeneticsoncology.org/Indexbyalpha/idxa_', x, '.html'))
  i <- grep('^<TABLE', p)
  message(paste(x, ' - ', length(p), ' text lines retrieved. TABLE indexes: ', paste(i, collapse = ',')))
  o <- p[(i[2]):(i[3]-5)]
  m <- str_match_all(o, '>\\s*([^<^\\(^\\s)]+).+</TR>')
  m <- unlist(lapply(m, '[', 2))
  m[!is.na(m)]
}))))

lymphoma <- toupper(read.table(url('http://www.bushmanlab.org/assets/doc/humanLymph.tsv'), header = TRUE, sep='\t')$symbol)
sanger   <- toupper(read.table('SangerCosmid_20180502.tsv', header = TRUE, sep = '\t')$Gene.Symbol)

## PMID23539594 <- unique(toupper(unlist(lapply(list.files(pattern = 'PMID23539594'), function(x){ scan(x, what = 'character') }))))
## CanGenes <- unique(toupper(unlist(lapply(list.files(pattern = 'PMID16959974'), function(x){ scan(x, what = 'character') }))))
## misc  <- unique(toupper(read.table(url('http://www.bushmanlab.org/assets/doc/miscellaneous.tsv'), header = TRUE, sep = '\t')$symbol))


# Read in a custom HGNC archive file and split appart both the previous symbol and synonym columns.
hgnc <- read.table('HGNC.20180503.tsv', sep = '\t', header = TRUE, quote = '', fill = TRUE, strip.white = TRUE, comment.char = '')
hgnc$Previous.Symbols <- stringr::str_split(hgnc$Previous.Symbols, '\\s*,\\s*')
hgnc <- unnest(hgnc, Previous.Symbols)
hgnc$Synonyms <- str_split(hgnc$Synonyms, '\\s*,\\s*')
hgnc <- unnest(hgnc, Synonyms)

hgnc$Approved.Symbol  <- toupper(hgnc$Approved.Symbol)
hgnc$Previous.Symbols <- toupper(hgnc$Previous.Symbols)
hgnc$Synonyms         <- toupper(hgnc$Synonyms)

prevAllOnco <- read.table('allOnco_Feb2017.tsv', sep = '\t', header = TRUE, quote = '', fill = TRUE, strip.white = TRUE, comment.char = '')

# AMP19  From Atlas: AML1 (now RUNX1) fusion.  RUNX1 on onocoGeneList
# LMDRA  From Atlas: HGNC name: LRMDA
# FGA7   From Atlas: FGA7 (fused gene 7 to AML1(now RUNX1)) RUNX1 on onocoGeneList

allOnco <- unique(toupper(c(prevAllOnco$symbol, atlas, lymphoma, sanger)))
allOnco <- c(allOnco[! allOnco %in% c('AMP19', 'LMDRA', 'FGA7')], 'LRMDA')

allOnco <- data.frame(symbol = allOnco,
                      prevSymbols = NA,
                      synonyms = NA,
                      name = NA,
                      refSeqName = NA)

# Update previous symbols with current symbols.
invisible(lapply(allOnco[! allOnco$symbol %in% hgnc$Approved.Symbol,]$symbol, function(x){
  g <- hgnc[match(x, hgnc$Previous.Symbols),]$Approved.Symbol
  allOnco[match(x, allOnco$symbol),]$symbol <<- g
}))

# Symbol updates may lead to duplicates.
allOnco <- allOnco[! duplicated(allOnco),]



allOnco$prevSymbols  <- sapply(allOnco$symbol, function(x){
  g <- subset(hgnc, hgnc$Approved.Symbol == x)$Previous.Symbols
  g <- unique(g[nchar(g) > 0])
  paste0(g, collapse = ',')
})

allOnco$synonyms <- sapply(allOnco$symbol, function(x){
  g <- subset(hgnc, hgnc$Approved.Symbol == x)$Synonyms
  g <- unique(g[nchar(g) > 0])
  paste0(g, collapse = ',')
})

allOnco$name <- hgnc[match(allOnco$symbol, hgnc$Approved.Symbol),]$Approved.Name


ncbiRefSeq <- read.table('ncbiRefSeq.txt', '\t', header = FALSE)
names(ncbiRefSeq) <- c('bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts',
              'exonEnds', 'score', 'name2', 'cdsStartStat', 'cdsEndStat',
              'exonFrames')

allOnco$refSeqName <- unname(sapply(allOnco$symbol, function(x){
  if(x %in% toupper(ncbiRefSeq$name2)){
    return(x)
  } else {
    prev <- base::intersect(subset(hgnc, Approved.Symbol == x)$Previous.Symbols, toupper(ncbiRefSeq$name2))[1]
    syn  <- base::intersect(subset(hgnc, Approved.Symbol == x)$Synonyms, toupper(ncbiRefSeq$name2))[1]
    if(! is.na(prev)) return(prev)
    # if(! is.na(syn)) return(syn)
    return(NA)
  }
}))

write.table(allOnco, file = 'allOncoNew.tsv', sep = '\t')


