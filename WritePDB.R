# Writing AMBER compliant PDB structure complexes
# This code works ONLY if your structure files (PDB) are in the working directory AND the filename is a a single capital letter [A-Z] with the extension .pdb
library("Rpdb")
library("gWidgets")
options(guiToolkit="RGtk2")

readpdb <- function(x){
  read.pdb(x, ATOM = TRUE, HETATM = TRUE, CRYST1 = TRUE, CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1)
}

fils <- list.files(path = ".", pattern = "[A-Z]{1}.pdb")
pdbs <- lapply(fils, readpdb)

savepdb <- function(...){
  write.pdb(pdbs, file = "testcomplex.pdb")
}

buttonpress <- function(h, ...){
  rgl.clear(type = "all")
  struc1 <- Txyz(pdbs[[1]], x = svalue(xval), y = svalue(yval), z = svalue(zval))
  print(svalue(zval))
  struc2 <- merge(struc1, pdbs[[2]])
  pdbs <- pdbs[-2]
  pdbs[[1]] <- struc2
  pdbs <<- pdbs
}
redraw <- function(...){
  a <- Txyz(pdbs[[1]], x = svalue(xval))
  aggr <- merge(a, pdbs[[2]])
  return(c(rgl.clear(type = "all"), visualize(aggr, add = TRUE)))
}

redraw_y <- function(...){
  a <- Txyz(pdbs[[1]], y = svalue(yval))
  aggr <- merge(a, pdbs[[2]])
  return(c(rgl.clear(type = "all"), visualize(aggr, add = TRUE)))
}

redraw_z <- function(...){
  a <- Txyz(pdbs[[1]], z = svalue(zval))
  aggr <- merge(a, pdbs[[2]])
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

merger <- gaction("Merge", handler=buttonpress)
mergebtn <- gbutton(action=merger, cont=group)

saves <- gaction("Save", handler=savepdb)
mergebtn <- gbutton(action=saves, cont=group)

visible(window) <- TRUE
