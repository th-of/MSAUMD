#!/bin/bash
set -e
#export AMBERHOME=/home/thomas/amber18
#source $AMBERHOME/amber.sh

#Unpacked structure file:
#pdb=unpacked.pdb
## packmol-memgen --pdb ${pdb} --lipids DOPE:DOPG --ratio 3:1

pdbfilename=bilayer_preprod_gaka.pdb
runname=GakA_mem_lastframe

## Find the number of atoms in peptide (for ntwprt):

## Generate AMBER input topology (prmtop) and input coordinates (inpcrd)
cat << EOF > system.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc 				
source leaprc.water.opc
source leaprc.lipid17
system = loadpdb ${pdbfilename}
set system box { 79.41899871826172 79.981998443603514 87.02600097656250 }
saveamberparm system ${runname}.prmtop ${runname}.inpcrd
quit
EOF
tleap -f system.in

a=$(grep -n -m 1 "TER" ${pdbfilename} |cut -f1 -d:)
atomcount=$(expr $a - 3)
echo $atomcount

## AMBER minimization input
cat <<EOF > min_1.in
minimisation ${runname}_1
&cntrl
  ig=-1
  iwrap=1,ntr=1
  restraint_wt=1.0, restraintmask='!(:WAT,K+)',
  imin=1,ntmin=1,maxcyc=20000,ncyc=2500,
  cut=10,ntb=1,
  ischeme=1,ithermostat=1,therm_par=1,
  ntc=1,ntf=1,
  ntpr=10,
 /
EOF

cat <<EOF > min_2.in
minimisation ${runname}_2
&cntrl
  ig=-1
  iwrap=1,ntr=1
  restraint_wt=1.0, restraintmask='(:1-34)',
  imin=1,ntmin=1,maxcyc=20000,ncyc=2500,
  cut=10,ntb=1,
  ntc=1,ntf=1,
  ntpr=10,
 /
EOF

cat <<EOF > min_2.in
minimisation ${runname}_2
&cntrl
  ig=-1
  iwrap=1,ntr=1
  restraint_wt=1.0, restraintmask='!(:WAT,K+)',
  imin=1,ntmin=1,maxcyc=20000,ncyc=2500,
  cut=10,ntb=1,
  ntc=1,ntf=1,
  ntpr=10,
 /
EOF

cat <<EOF > min_5.in
minimisation ${runname}_5
&cntrl
  ig=-1
  iwrap=1
  imin=1,ntmin=0,maxcyc=20000,ncyc=2500,
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
  nstlim=2500,dt=0.001,
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
echo 'Production minimization'
mpirun -np 4 sander.MPI -O -i min.in -o minmdout -p ${runname}.prmtop -c ${runname}.inpcrd -r ${runname}_min.xyz -x ${runname}_min.nc
date +"%T"
echo 'Production heating'
sander -O -i heat.in -o heatmdout -p ${runname}.prmtop -c ${runname}_min.xyz -r ${runname}_heat.xyz -x ${runname}_heat.nc
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

for j in {2..20}
do
    date +"%T"
    echo 'Production run' $j
	pmemd.cuda -O -i prod.in -o prod"$j"mdout -p ${runname}.prmtop -c ${runname}_prod"$(($j-1))".xyz -r ${runname}_prod"$j".xyz -x ${runname}_prod"$j".nc 
done
