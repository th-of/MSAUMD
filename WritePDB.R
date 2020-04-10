# Writing AMBER compliant PDB structure complexes
# This code works ONLY if your structure files (PDB) are in the working directory AND the filename is a a single capital letter [A-Z] with the extension .pdb
library("Rpdb")
library("gWidgets")
options(guiToolkit="RGtk2")


for (i in list.files(path = ".", pattern = "[A-Z]{1}.pdb")){
  assign(paste0(i), read.pdb(i, ATOM = TRUE, HETATM = TRUE, CRYST1 = TRUE, CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1))
}

redraw <- function(...){
  a <- Txyz(A.pdb, x = svalue(xval))
  aggr <- merge(a, B.pdb)
  return(visualize(aggr, type = "l"))
}

window <- gwindow("Translation along the x-axis", visible=FALSE)
group <- ggroup(cont = window, expand = TRUE)
label <- glabel("x-axis", cont = group)
xval <- gslider(from=-100, to=100, by=5, value=0, expand=TRUE, container = group, handler = redraw)
visible(window) <- TRUE


# gakb[["atoms"]][["chainid"]] <- "B"
# 
# gakab <- merge(A.pdb, B.pdb)
# 
# visualize(gakab, type = "s")
# 
# write.pdb(gakab, file = "/home/thomas/MD/GarKS complex/gakab_fix.pdb")
# 
# write("END", file = "/home/thomas/MD/GarKS complex/C_.pdb", append = TRUE)
# 
# complex <- read.pdb("/home/thomas/MD/GarKS complex/complex.pdb", ATOM = TRUE, HETATM = TRUE, CRYST1 = TRUE, CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1)
# 
# complex <- merge(gaka, complex)
# 
# visualize(complex)
# 
# gakc_trans <- Txyz(gakc, x=20)


