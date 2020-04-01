#!/bin/bash

export AMBERHOME=/home/thomas/amber18
source $AMBERHOME/amber.sh

## name of MD run, model or project; CHANGE SEQUENCE IN TLEAP !!
runname=N_GakA

## GakA, GakB or GakC
gakpeptide=GakB

cat <<EOF > gakselect.R
gakselect <- function(){					
args <- commandArgs(trailingOnly = TRUE)
library("seqinr")
if (args[1] == "GakA"){
  gak <- toupper(aaa(strsplit("MGAIIKAGAKIVGKGVLGGGASWLGWNVGEKIWK", "")[[1]]))
} else if (args[1] == "GakB"){
  gak <- toupper(aaa(strsplit("MGAIIKAGAKIIGKGLLGGAAGGATYGGLKKIFG", "")[[1]]))
} else {
  gak <- toupper(aaa(strsplit("MGAIIKAGAKIVGKGALTGGGVWLAEKLFGGK", "")[[1]]))
}
cat(" ", paste0("N", gak[1]), gak[2:(length(gak)-1)], paste0("C", tail(gak, 1)), " ")
}
gakselect()
EOF

peptide=$(Rscript gakselect.R $gakpeptide)

## Build preproduction system:
cat <<EOF > ${runname}_preprod.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc 				
source leaprc.water.opc							
${runname} = sequence {$peptide}									      
solvateShell ${runname} OPCBOX 20.0 0.8				        
addions ${runname} Cl- 0
savepdb ${runname} ${runname}_atomcount.pdb
saveamberparm ${runname} ${runname}_preprod.prmtop ${runname}_preprod.inpcrd					
quit
EOF
tleap -f ${runname}_preprod.in

a=$(grep -n -m 1 "TER" ${runname}_atomcount.pdb |cut -f1 -d:)
atomcount=$(expr $a - 2)
echo 'Atoms in pdb' $atomcount

mv leap.log preprodsalt_leap.log

## min.in
cat <<EOF > min.in
minimise ${runname}						
 &cntrl							
  ig=-1							
  iwrap =1					
  imin=1,ntmin=0, maxcyc=2500,ncyc=500,			
  cut=10,ntb=1,				 		
  ntc=1,ntf=1,					
  ntpr=10,						
 /						
EOF

## heat.in
cat <<EOF > heat.in
heat ${runname}							
 &cntrl							
  ig=-1							
  iwrap =1					
  imin=0,irest=0,ntx=1,					
  nstlim=20000,dt=0.002,					
  ntc=2,ntf=2,					
  cut=12, ntb=1, ntp=0,					
  ntpr=500, ntwx=500, ntwr=1000,			
  ntt=3, gamma_ln=1.0,				
  tempi=0.0, temp0=300.0,				
  tol=1.0e-8,jfastw=0, nmropt=1,			
 /							
 &ewald						
  dsum_tol=1.0e-6,				 	
  order=4, skinnb=2.0, vdwmeth=1,			
 /							
 &wt TYPE='TEMP0', istep1=0, istep2=20000,		
  value1=0.1, value2=300.0, /			
 &wt TYPE='END' /
EOF

## density.in
cat <<EOF > density.in
density ${runname}						
 &cntrl							
  ig=-1							
  iwrap =1					
  imin=0,irest=1,ntx=5,					
  nstlim=250000,dt=0.002,				
  ntc=2,ntf=2,					
  cut=12, ntb=2, ntp=1, taup=1.0,			
  ntpr=500, ntwx=500, ntwr=1000, ntwprt=$atomcount,	
  ntt=3, gamma_ln=1.0,				
  temp0=300.0,						
  tol=1.0e-8,jfastw=0,				
 /						
 &ewald							
  dsum_tol=1.0e-6,					
  order=4, skinnb=2.0, vdwmeth=1,			
 /
EOF

## equil.in
cat <<EOF > equil.in
equil ${runname}						
 &cntrl							
  ig=-1							
  iwrap =1					
  imin=0,irest=1,ntx=5,					
  nstlim=250000,dt=0.002,				
  ntc=2,ntf=2,					
  cut=12, ntb=1, ntp=0, tautp=10.0,			
  ntpr=500, ntwx=500, ntwr=1000, ntwprt=$atomcount,	
  ntt=1,							
  temp0=300.0,					
  tol=1.0e-8,jfastw=0,				
 /						
 &ewald							
  dsum_tol=1.0e-6,				
  order=4, skinnb=2.0, vdwmeth=1,			
 /						
EOF

## preprod.in
cat <<EOF > preprod.in
equil ${runname}									
 &cntrl										
  ig=-1										
  iwrap =1									
  imin=0,irest=1,ntx=5,								
  nstlim=2500000,dt=0.002,							
  ntc=2,ntf=2,								
  cut=12, ntb=1, ntp=0, tautp=10.0,						
  ntpr=1000, ntwx=1000, ntwr=10000, ntwv=-1, ioutfm=1, ntwprt=$atomcount,	
  ntt=1,										
  temp0=300.0,									
  tol=1.0e-8,jfastw=0,								
 /									
 &ewald										
  dsum_tol=1.0e-6,								
  order=4, skinnb=2.0, vdwmeth=1,						
 /
EOF

## prod.in
cat <<EOF > prod.in	
equil ${runname}									
 &cntrl										
  ig=-1										
  iwrap =1									
  imin=0,irest=1,ntx=5,								
  nstlim=2500000,dt=0.002,							
  ntc=2,ntf=2,									
  cut=12, ntb=1, ntp=0, tautp=10.0,						
  ntpr=1000, ntwx=1000, ntwr=10000, ntwv=-1, ioutfm=1, ntwprt=$atomcount,
  ntt=1,									
  temp0=300.0,								
  tol=1.0e-8,jfastw=0,							
 /										
 &ewald										
  dsum_tol=1.0e-6,								
  order=4, skinnb=2.0, vdwmeth=1,						
 /
EOF
date +"%T"
echo 'Start preproduction minimization'
pmemd.cuda -O -i min.in -o preminmdout -p ${runname}_preprod.prmtop -c ${runname}_preprod.inpcrd -r ${runname}_premin.xyz
date +"%T"
echo 'Preproduction heating'
pmemd.cuda -O -i heat.in -o preheatmdout -p ${runname}_preprod.prmtop -c ${runname}_premin.xyz  -r ${runname}_preheat.xyz -x ${runname}_preheat.nc 
date +"%T"
echo 'Preproduction density equilibriation'
pmemd.cuda -O -i density.in -o predensitymdout -p ${runname}_preprod.prmtop -c ${runname}_preheat.xyz -r ${runname}_predens.xyz -x ${runname}_predens.nc
date +"%T"
echo 'Preproduction equilibriation'
pmemd.cuda -O -i equil.in -o preequilmdout -p ${runname}_preprod.prmtop -c ${runname}_predens.xyz -r ${runname}_preequil.xyz -x ${runname}_preequil.nc 
date +"%T"
echo 'Preproduction'
pmemd.cuda -O -i preprod.in -o preprodmdout -p ${runname}_preprod.prmtop -c ${runname}_preequil.xyz -x ${runname}_preprod.nc 

# extract the preproduction structure for production run:
cat <<EOF > extract.in
parm ${runname}_preprod.prmtop				
parmstrip :WAT,Cl-				
trajin ${runname}_preprod.nc					
trajout ${runname}.pdb pdb onlyframes 2500 		
run
EOF
cpptraj -i extract.in

## Resolvate
cat <<EOF > solvated.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc				
source leaprc.water.opc						
${runname} = loadpdb ${runname}.pdb						
solvateOct ${runname} OPCBOX 8.0 						
addions ${runname} Cl- 0			 	
saveamberparm ${runname} ${runname}.prmtop ${runname}.inpcrd
savepdb ${runname} ${runname}_production_structure.pdb						
quit
EOF
tleap -f solvated.in

#mv tleap.log tleap_log_prod.log

date +"%T"
echo 'Production minimization'
pmemd.cuda -O -i min.in -o minmdout -p ${runname}.prmtop -c ${runname}.inpcrd -r ${runname}_min.xyz -x ${runname}_min.nc
date +"%T"
echo 'Production heating'
pmemd.cuda -O -i heat.in -o heatmdout -p ${runname}.prmtop -c ${runname}_min.xyz -r ${runname}_heat.xyz -x ${runname}_heat.nc
date +"%T"
echo 'Production density equilibriation'
pmemd.cuda -O -i density.in -o densitymdout -p ${runname}.prmtop -c ${runname}_heat.xyz -r ${runname}_dens.xyz -x ${runname}_dens.nc
date +"%T"
echo 'Production equilibriation'
pmemd.cuda -O -i equil.in -o equilmdout -p ${runname}.prmtop -c ${runname}_dens.xyz -r ${runname}_equil.xyz -x ${runname}_equil.nc
date +"%T"
echo 'Production run 1'
pmemd.cuda -O -i prod.in -o prodmdout -p ${runname}.prmtop -c ${runname}_equil.xyz -r ${runname}_prod1.xyz -x ${runname}_prod1.nc 

# Repeat

for j in {2..210}
do
    date +"%T"
    echo 'Production run' $j
	pmemd.cuda -O -i prod.in -o prod"$j"mdout -p ${runname}.prmtop -c ${runname}_prod"$(($j-1))".xyz -r ${runname}_prod"$j".xyz -x ${runname}_prod"$j".nc 
done


