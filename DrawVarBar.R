library(tidyverse)
library(R.utils)
library(cowplot)
library(grid)

init.seqence.bar <-
  function(len, seqwidth){
    # draw a sequence
    xincre <- 1/len
    grid.rect(x=0.5, y=0.5, width=1, height=seqwidth)

    vp <- viewport(width = 0.95, height = 0.95)
    pushViewport(vp)

    #functions to draw in sequence
    create.site.view <-
      function(site){
        site.vp <- viewport(x=xincre*(site-0.5), y=0.5, width=xincre, height=seqwidth)
        return(site.vp)
      }

    add.object.site <-
      function(object, site){
        site.vp <- create.site.view(site)
        pushViewport(site.vp)
        object()
        upViewport()
      }

    assign("add.object.site", add.object.site, envir = .GlobalEnv)
  }

redbar <-
  function() grid.rect(x=0.5, y=0.5, width=1, height=1, gp=gpar(fill="red", col="red"))

greenbar <-
  function() grid.rect(x=0.5, y=0.5, width=1, height=1, gp=gpar(fill="green", col="green"))

bluebar <-
  function() grid.rect(x=0.5, y=0.5, width=1, height=1, gp=gpar(fill="blue", col="blue"))

str_vec <-
  function(x){str_split(x,"")[[1]]}

read_ali <-
  function(inpath, checkalign=TRUE){
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

DrawSeqBar_denovo <-
  function(inpath, revertCol = FALSE){
    #import sequence alignment
    ali <- read_ali(inpath)

    #check number of sequence
    if (attr(ali,"n") != 3) stop("wrong number of sequence")
    n <- attr(ali,"n")
    len <- attr(ali,"length")

    #put sequences into table
    sequence_table <-
      tibble(site=1:len,
             ref=str_vec(ali$seq[1]),
             seq1=str_vec(ali$seq[2]),
             seq2=str_vec(ali$seq[3]))

    #generate
    sequence_table %<>%
      filter((ref==seq1)|(ref==seq2)) %>%
      filter(!((ref==seq1)&(seq1==seq2))) %>%
      rowwise() %>%
      mutate(color=ifelse(ref==seq1, 1, 2))

    if (revertCol) sequence_table %<>% mutate(color=ifelse(ref==seq1, 2, 1))

    #crate sequence bar
    grid.newpage()
    vp <- viewport(width = .9)
    pushViewport(vp)
    init.seqence.bar(len, 0.2)

    #add bar
    mapply(FUN=function(site, color){
      if (color==1) add.object.site(redbar, site)
      if (color==2) add.object.site(greenbar, site)
    },
    site=sequence_table$site,
    color=sequence_table$color,
    SIMPLIFY=FALSE)

    seq_graph <- grid.grab()
    return(seq_graph)
  }


DrawSeqBar_ref <-
  function(inpath, revertCol = FALSE){
    #import sequence alignment
    ali <- read_ali(inpath)

    #check number of sequence
    if (attr(ali,"n") != 2) stop("wrong number of sequence")
    n <- attr(ali,"n")
    len <- attr(ali,"length")

    #put sequences into table
    sequence_table <-
      tibble(site=1:len,
             ref=str_vec(ali$seq[1]),
             seq1=str_vec(ali$seq[2]))

    #generate
    sequence_table %<>%
      filter(ref!=seq1) %>%
      rowwise() %>%
      mutate(color=1)

    if (revertCol) sequence_table %<>% mutate(color=ifelse(ref==seq1, 2, 1))

    #crate sequence bar
    grid.newpage()
    vp <- viewport(width = .9)
    pushViewport(vp)
    init.seqence.bar(len, 0.2)

    #add bar
    mapply(FUN=function(site, color){
      if (color==1) add.object.site(bluebar, site)
    },
    site=sequence_table$site,
    color=sequence_table$color,
    SIMPLIFY=FALSE)

    seq_graph <- grid.grab()
    return(seq_graph)
  }

DrawSeqBar_both <-
  function(inpath, revertCol = FALSE){

    #import sequence alignment
    ali <- read_ali(inpath)

    #check number of sequence
    if (attr(ali,"n") != 4) stop("wrong number of sequence")
    n <- attr(ali,"n")
    len <- attr(ali,"length")

    # all counts
    all_table <-
      tibble(site=1:len,
             ref=str_vec(ali$seq[1]),
             refalt=str_vec(ali$seq[2]),
             seq1=str_vec(ali$seq[3]),
             seq2=str_vec(ali$seq[4])
      )

    counts <-
      all_table %>%
      mutate(refvar=ref!=refalt,
             denovovar=seq1!=seq2) %>%
      mutate(not_ref=!refvar & denovovar,
             not_denovovar=refvar & !denovovar
      ) %>%
      summarise(n_refvar=sum(refvar),
                n_denovovar=sum(denovovar),
                n_not_ref=sum(not_ref),
                n_not_denovovar=sum(not_denovovar)
      ) %>%
      mutate(in_ref=round(1-n_not_ref/n_denovovar, 3)*100,
             in_denovo=round(1-n_not_denovovar/n_refvar, 3)*100)


    #put sequences into table
    ref_table <-
      tibble(site=1:len,
             ref=str_vec(ali$seq[1]),
             seq1=str_vec(ali$seq[2]))

    #generate
    ref_table %<>%
      filter(ref!=seq1) %>%
      rowwise() %>%
      mutate(color=1)

    if (revertCol) ref_table %<>% mutate(color=ifelse(ref==seq1, 2, 1))

    #crate sequence bar
    grid.newpage()
    vp <- viewport(width = .6)

    grid.text(label =
                paste0(counts$n_refvar, "(", counts$in_ref, "%)"),
              x = 0.81, y=0.5, just = "left",
              gp = gpar(fontsize=10))

    grid.text(label = "mapping",
              x = 0.19, y=0.5, just = "right",
              gp = gpar(fontsize=10))

    pushViewport(vp)
    init.seqence.bar(len, 0.8)


    #add bar
    mapply(FUN=function(site, color){
      if (color==1) add.object.site(bluebar, site)
    },
    site=ref_table$site,
    color=ref_table$color,
    SIMPLIFY=FALSE)

    ref_graph <- grid.grab()


    #put sequences into table
    denovo_table <-
      tibble(site=1:len,
             ref=str_vec(ali$seq[1]),
             seq1=str_vec(ali$seq[3]),
             seq2=str_vec(ali$seq[4]))

    #generate
    denovo_table %<>%
      filter((ref==seq1)|(ref==seq2)) %>%
      filter(!((ref==seq1)&(seq1==seq2))) %>%
      rowwise() %>%
      mutate(color=ifelse(ref==seq1, 1, 2))

    if (revertCol) denovo_table %<>% mutate(color=ifelse(ref==seq1, 2, 1))

    #crate sequence bar
    grid.newpage()
    vp <- viewport(width = .6)

    grid.text(label =
                counts$n_denovovar,
              x = 0.81, y=0.5, just = "left",
              gp = gpar(fontsize=10))

    grid.text(label =
                paste0(counts$n_denovovar, "(", counts$in_denovo, "%)"),
              x = 0.81, y=0.5, just = "left",
              gp = gpar(fontsize=10))

    grid.text(label = "de novo",
              x = 0.19, y=0.5, just = "right",
              gp = gpar(fontface="italic", fontsize=10))

    pushViewport(vp)
    init.seqence.bar(len, .8)

    #add bar
    mapply(FUN=function(site, color){
      if (color==1) add.object.site(redbar, site)
      if (color==2) add.object.site(greenbar, site)
    },
    site=denovo_table$site,
    color=denovo_table$color,
    SIMPLIFY=FALSE)

    denovo_graph <- grid.grab()

    return(list(refplot=ref_graph, denovoplot=denovo_graph))
  }
