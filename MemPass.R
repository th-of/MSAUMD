## Mempass
# cpptraj
# parm prmtop
# trajin nc
# trajout filename.pdb pdb multi start x stop y
library("bio3d")

files <- list.files("./pdbs")
files <- files[order(nchar(files), files)]
l <- vector("list", length(files))

## Calculate the thickness of the membrane using a mean cluster analysis on phosphorus atoms
a <- read.pdb(paste0("./pdbs/", files[1]), ATOM.only = TRUE)
frame <- as.data.frame(a$atom)
# frame2 <- frame[frame$resid == 'OPE' | frame$resid == 'OPG', ]
frame3 <- frame[frame$elety == 'P', ]
fit <- kmeans(frame3$z, 2)
memthick <- aggregate(frame3$z, by=list(fit$cluster), FUN=mean)$x[1:2]

#
{
# x_min <- min(frame2$x)
# x_max <- max(frame2$x)
# 
# y_min <- min(frame2$y)
# y_max <- max(frame2$y)
# 
# z_min <- min(frame2$z)
# z_max <- max(frame2$z)
# 
# P_1 <- c(x_max, y_min, z_min)
# P_2 <- c(x_max, y_max, z_min)
# P_3 <- c(x_min, y_max, z_min)
# P_4 <- c(x_min, y_min, z_min)
# 
# P_5 <- c(x_max, y_min, z_max)
# P_6 <- c(x_max, y_max, z_max)
# P_7 <- c(x_min, y_max, z_max)
# P_8 <- c(x_min, y_min, z_max)

}
#

for (i in 1:length(files)){
  a <- read.pdb(paste0("./pdbs/", files[i]), ATOM.only = TRUE)
  frame <- as.data.frame(a$atom)
  frame <- frame[which((frame$elety == 'O') & (frame$resid == 'WAT')),]
  l[[i]] <- frame
}

counts <- c()
count <- l[[50]][13259,][1,11]
keep <- c()

for (n in 1:length(l[[1]]$z)){
  for (m in 1:length(l)){
    counts[m] <- l[[m]][n,11]
  }
  if (any(counts < 25) && any(counts > 45) && any(counts > 25 & counts < 45)){
    counts[101] <- as.integer(row.names(l[[m]][n,]))
    keep <- rbind(keep, counts, deparse.level = 0)
  }
}

for (k in 1:(length(keep[,1]))){
  keepNotJumped <- keep[k,1:99] - keep[k,2:100]
  if (any(abs(keepNotJumped) > 10) == FALSE){
    print(keep[k,101])
  }
}
keepNotJumped <- keep[1,1:99] - keep[1,2:100]


# waters <- l[[1:100]]$z
# criteria <- any(count < 20) && any(count > 40) && any(count > 20 & count < 40)
# criteria <- any(vector < 20) && any(vector > 40) && any(vector > 20 & vector < 40)

# criteria <- sum(vector[vector < 15 & vector > 45 & (vector > 15 & vector < 45)]) > 3
#Water molecule atom with the lowest z
# which.min(l[[1]]$z) = 34878
# The original index value of this atom corresponding to 'serial' in VMD:
# row.names(l[[1]][34878,])

