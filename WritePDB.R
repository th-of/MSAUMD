# Writing AMBER compliant PDB structure complexes
library("Rpdb")
struc <- read.pdb(, ATOM = TRUE, HETATM = TRUE, CRYST1 = TRUE, CONECT = TRUE, TITLE = TRUE, REMARK = TRUE, MODEL = 1)