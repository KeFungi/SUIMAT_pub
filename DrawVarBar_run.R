source("DrawVarBar.R")

haplolist <- c("Suiame1", "Suibr2", "Suigr1", "Rhivi1", "Rhisa1")

haplolist_name <-
  c("S. americanus", "S. brevipes", "S. weaverae", "R. vinicolor", "R. salebrosus") %>%
  map(function(x){
    grid.newpage()

    grid.text(label = x,
              x = 0.2, y=0.4, just = "left",
              gp = gpar(fontsize=14, fontface="bold.italic"))
    name_graph <- grid.grab()
  }
  )

inpath_list <-
  sapply(haplolist,
         FUN=function(x) paste0("HDalignments/", x, "_MATA_allseq_ali.fasta"))

cowlist <- list()

for (i in seq_along(inpath_list)){
  newplot <- DrawSeqBar_both(inpath_list[[i]])
  cowlist <- c(cowlist, list(haplolist_name[[i]]), list(newplot$refplot), list(newplot$denovoplot))
}

plot_grid(plotlist = cowlist,
          greedy = FALSE,
          ncol = 1)

ggsave2("VarBars.pdf", height = 5, width = 5.5)




