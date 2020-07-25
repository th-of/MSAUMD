## Mempass
# cpptraj
# parm prmtop
# trajin nc
# trajin filename.pdb pdb multi start x stop y
# https://math.stackexchange.com/questions/1472049/check-if-a-point-is-inside-a-rectangular-shaped-area-3d  
library("bio3d")

files <- list.files("./pdbs")
files <- files[order(nchar(files), files)]
l <- vector("list", length(files))

## Calculate the vertices of a rectangle formed by the membrane lipids (from the first pdb only)
a <- read.pdb(paste0("./pdbs/", files[1]), ATOM.only = TRUE)
frame <- as.data.frame(a$atom)
frame2 <- frame[frame$resid == 'OPE' | frame$resid == 'OPG', ]

x_min <- min(frame2$x)
x_max <- max(frame2$x)

y_min <- min(frame2$y)
y_max <- max(frame2$y)

z_min <- min(frame2$z)
z_max <- max(frame2$z)

P_1 <- c(x_max, y_min, z_min)
P_2 <- c(x_max, y_max, z_min)
P_3 <- c(x_min, y_max, z_min)
P_4 <- c(x_min, y_min, z_min)

P_5 <- c(x_max, y_min, z_max)
P_6 <- c(x_max, y_max, z_max)
P_7 <- c(x_min, y_max, z_max)
P_8 <- c(x_min, y_min, z_max)


for (i in 1:length(files)){
  a <- read.pdb(paste0("./pdbs/", files[i]), ATOM.only = TRUE)
  frame <- as.data.frame(a$atom)
  frame <- frame[which(frame$resid == 'WAT'),]
  l[[i]] <- frame
}

#Water molecule atom with the lowest z
# which.min(l[[1]]$z) = 34878
# The original index value of this atom corresponding to 'serial' in VMD:
# row.names(l[[1]][34878,])

