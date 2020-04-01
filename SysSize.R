# Calculate size of system
# Thomas Oftedal (thof@nmbu.no)

library("Rpdb")
struc <- read.pdb("bilayer_GakA.pdb", ATOM = TRUE)
newdata <- struc$atoms[which(struc$atoms$resname == 'WAT'), 9:11]
cat(sprintf("%.6f %.6f %.6f", max(newdata$x1)-min(newdata$x1), max(newdata$x2)-min(newdata$x2), max(newdata$x3)-min(newdata$x3)))

