# Calculate size of system
# Thomas Oftedal (thof@nmbu.no)
syssize <- function(){					
args <- commandArgs(trailingOnly = TRUE)

library("Rpdb")
struc <- read.pdb(args[1], ATOM = TRUE)
newdata <- struc$atoms[which(struc$atoms$resname == 'WAT'), 9:11]
cat(sprintf("%.6f %.6f %.6f", max(newdata$x1)-min(newdata$x1), max(newdata$x2)-min(newdata$x2), max(newdata$x3)-min(newdata$x3)))
}
syssize()