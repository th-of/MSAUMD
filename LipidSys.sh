#!/bin/bash
set -e
#export AMBERHOME=/home/thomas/amber18
#source $AMBERHOME/amber.sh

#Unpacked structure file:
#pdb=unpacked.pdb
## packmol-memgen --pdb ${pdb} --lipids DOPE:DOPG --ratio 3:1

pdbfilename=GakA.pdb
runname=GakA_mem

# Generate the lipid bilayer system:

packmol-memgen --pdb ${pdbfilename} --lipids DOPE:DOPG --ratio 3:1 --overwrite

## Calculate box size:

dims=$(./vmd_box_dims.sh -i bilayer_${pdbfilename} -s water)
dim1=$(echo $dims | cut -d , -f 1)
dim2=$(echo $dims | cut -d , -f 2)
dim3=$(echo $dims | cut -d , -f 3)

## Generate AMBER input topology (prmtop) and input coordinates (inpcrd)
cat << EOF > system.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc 				
source leaprc.water.opc
source leaprc.lipid17
system = loadpdb ${pdbfilename}
set system box { ${dim1}${dim2}${dim3} }
saveamberparm system ${runname}.prmtop ${runname}.inpcrd
quit
EOF
tleap -f system.in

a=$(grep -n -m 1 "TER" ${pdbfilename} |cut -f1 -d:)
atomcount=$(expr $a - 3)
echo $atomcount

## minimization
cat <<EOF > min.in
Memgen AMBER minimization with restraints
 &cntrl
  nmropt = 0,
  imin = 1,
  maxcyc = 500,
  ncyc = 250,
 &end
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
