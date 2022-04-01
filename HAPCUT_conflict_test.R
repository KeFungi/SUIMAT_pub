library(tidyverse)
library(magrittr)

str_vec <- function(x){str_split(x,"")[[1]]}

#read fasta alignment from file
read_ali <- function(inpath, checkalign=TRUE){
  rawin <- read_lines(inpath)
  outtable <- tibble()
  nseq <- 0
  
  for ( l in rawin){
    if (startsWith(l,">") == TRUE){
      nseq <- nseq + 1
      seqname <- str_sub(l,2)
      seq <- ""
      
      outtable <- 
        rbind(outtable,
              tibble(seqname=seqname,seq=seq))
    }
    
    if (startsWith(l,">") != TRUE){
      outtable[nseq,2] = paste0(outtable[nseq,2], l)
    }
  }
  
  if (checkalign == TRUE) {
    seqlenlist <- outtable %>% mutate(seqlen=nchar(seq)) %>% pull(seqlen)
    if (all(seqlenlist==seqlenlist[1])){
      attr(outtable, "length") <- seqlenlist[1]
    }
    else {
      warning("sequence in the same length, not alignmened")
    }
  }
  
  attr(outtable, "n") <- nrow(outtable)
  return(outtable)
}

#transpose alignment to a table where the cols are seuqnce names, rows are
#sites in the sequence
ali_to_table <- function(ali, with_all_site=FALSE){
  newtable <- tibble(site = 1:attr(ali, "length"))
  
  for ( i in 1:nrow(ali)){
    newcol <- tibble(!!paste(ali$seqname[i]) := str_vec(ali$seq[i]))
    newtable <- cbind(newtable, newcol)
  }
  
  if (with_all_site){
    sitetable <- mutate_all(newtable, .funs = ~cumsum(.x!="-"))
    newtable <- left_join(newtable, sitetable, by = c("site"="site"), suffix = c("", ".site"))
  }
  
  return(newtable)
}

read_HAPCUT_output <- function(path){
  x <- readLines(path) %>% .[1] %>% str_split(pattern = " ") %>% .[[1]]
  outtable <- tibble(len=x[5], phased=x[7], span=x[9], fragments=x[11])
  return(outtable)
}

pad_len <- 1000
spcodel=c("Rhisa1", "Rhivi1", "Suiame1", "Suibr2", "Suigr1" )
report <- tibble()

for (spcode in spcodel) {
  fasta_inpath <- paste0("HAPCUT/", spcode, "_HapCut_compare.fasta")
  hapcut_inpath <- paste0("HAPCUT/", spcode, "_HapCut")
  ali <- read_ali(fasta_inpath,checkalign = TRUE)
  seqtable <- ali_to_table(ali)
  seqtable <-
    seqtable %>%
    set_colnames(c("site","hapcut1", "hapcut2", "denovo1", "denovo2"))
  
  seqtable_diff <-
    seqtable %>%
    .[-1:-pad_len,] %>%
    head(-pad_len) %>%
    filter(hapcut1!=hapcut2)
  
  comb1 <- seqtable_diff %>% filter((hapcut1==denovo2)&(hapcut2==denovo1))
  comb2 <- seqtable_diff %>% filter((hapcut1==denovo1)&(hapcut2==denovo2))
  bestcom <- list(comb1, comb2)[[which.min(c(nrow(comb1), nrow(comb2)))]]
  
  HAPCUT_output <-
    read_HAPCUT_output(hapcut_inpath)
  
  report <-
    tibble(spcode=spcode %>% str_extract("^[^_]+"),
           ali_len=nrow(seqtable),
           phased_site=nrow(seqtable_diff),
           conflict_site=nrow(bestcom)
           ) %>%
    bind_cols(HAPCUT_output) %>%
    bind_rows(report, .)
}

write_csv(report, "HAPCUT_summary.csv")
