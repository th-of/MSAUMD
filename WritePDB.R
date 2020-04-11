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
  return(c(rgl.clear(type = "all"), visualize(aggr, add = TRUE)))
}

redraw_y <- function(...){
  a <- Txyz(A.pdb, y = svalue(yval))
  aggr <- merge(a, B.pdb)
  return(c(rgl.clear(type = "all"), visualize(aggr, add = TRUE)))
}

redraw_z <- function(...){
  a <- Txyz(A.pdb, z = svalue(zval))
  aggr <- merge(a, B.pdb)
  return(c(rgl.clear(type = "all"), visualize(aggr, add = TRUE)))
}

z_view <- function(...){
  return(view3d(theta = 0, phi = 0))
}
x_view <- function(...){
  return(view3d(theta = 90, phi = 0))
}
y_view <- function(...){
  return(view3d(theta = 0, phi = 90))
}

window <- gwindow("Translation", visible=FALSE)
group <- ggroup(horizontal = FALSE, cont = window, expand = TRUE)

label <- glabel("x-axis", cont = group)
xval <- gslider(from=-100, to=100, by=5, value=0, expand=TRUE, container = group, handler = redraw)

label <- glabel("y-axis", cont = group)
yval <- gslider(from=-100, to=100, by=5, value=0, expand=TRUE, container = group, handler = redraw_y)

label <- glabel("z-axis", cont = group)
zval <- gslider(from=-100, to=100, by=5, value=0, expand=TRUE, container = group, handler = redraw_z)

z_set <- gbutton(text = "Z view", border = TRUE, container = group, handler = z_view)
x_set <- gbutton(text = "X view", border = TRUE, container = group, handler = x_view)
y_set <- gbutton(text = "Y view", border = TRUE, container = group, handler = y_view)

act <- gaction("merge", handler=function(...) {print(svalue(xval))})
b3 <- gbutton(action=act, cont=group)

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


