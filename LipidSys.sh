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

dims=$(Rscript SysSize.R bilayer_${pdbfilename})

## Generate AMBER input topology (prmtop) and input coordinates (inpcrd)
cat << EOF > system.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc 				
source leaprc.water.opc
source leaprc.lipid17
system = loadpdb bilayer_${pdbfilename}
set system box { ${dims} }
saveamberparm system ${runname}.prmtop ${runname}.inpcrd
quit
EOF
tleap -f system.in

#a=$(grep -n -m 1 "TER" bilayer_${pdbfilename} |cut -f1 -d:)
#atomcount=$(expr $a - 3)
#echo $atomcount

## minimization
cat <<EOF > min.in
Lipid minimization $(date)
 &cntrl
  imin=1,       ! Minimize the initial structure
  maxcyc=10000, ! Maximum number of cycles for minimization
  ncyc=5000,    ! Switch from steepest descent to conjugate gradient minimization after ncyc cycles
  ntb=1,        ! Constant volume
  ntp=0,        ! No pressure scaling
  ntf=1,        ! Complete force evaluation
  ntc=1,        ! No SHAKE
  ntpr=50,      ! Print to mdout every ntpr steps
  ntwr=2000,    ! Write a restart file every ntwr steps
  cut=10.0,     ! Nonbonded cutoff in Angstroms
 /
EOF

date +"%T"
echo 'Minimization'
mpirun -np 4 pmemd.MPI -O -i min.in -o minmdout -p ${runname}.prmtop -c ${runname}.inpcrd -r min.rst


## heat.in
cat <<EOF > heat1.in
Lipid heating $(date)
 &cntrl
  imin=0,                       ! Molecular dynamics
  ntx=1,                        ! Positions read formatted with no initial velocities
  irest=0,                      ! No restart
  ntc=2,                        ! SHAKE on for bonds with hydrogen
  ntf=2,                        ! No force evaluation for bonds with hydrogen
  tol=0.0000001,                ! SHAKE tolerance
  nstlim=2500,                  ! Number of MD steps
  ntt=3,                        ! Langevin thermostat
  gamma_ln=1.0,                 ! Collision frequency for Langevin thermostat
  ntr=1,                        ! Restrain atoms using a harmonic potential
  restraint_wt=10               ! degree of restraint
  restraintmask=':OL,PC,PGR',   ! Restrain all except oleoyl lipid tails / acyl chains
  ig=-1,                        ! Random seed for Langevin thermostat
  ntpr=100,
  ntwr=10000,
  ntwx=100,                     ! Write to trajectory file every ntwx steps
  dt=0.002,                     ! Timestep (ps)
  nmropt=1,                     ! NMR restraints will be read (See TEMP0 control below)
  ntb=1,
  ntp=0,
  cut=10.0,
  ioutfm=1,                     ! Write a binary (netcdf) trajectory
  ntxo=2,                       ! Write binary restart files
 /
 &wt 
  type='TEMP0',   ! Varies the target temperature TEMP0
  istep1=0,       ! Initial step
  istep2=2500,    ! Final step
  value1=0.0,     ! Initial temp0 (K)
  value2=100.0 /  ! final temp0 (K)
 &wt type='END' / ! End of varying conditions
Hold lipid fixed
RES :OL,PC,PGR
END
END               ! End GROUP input
EOF

date +"%T"
echo 'Heating 1'
pmemd.cuda -O -i heat1.in -o heatmdout -p ${runname}.prmtop -c min.rst -r heat.rst -ref min.rst -x heat.nc

## heat2.in
cat <<EOF > heat2.in	
heat2 ${date}
Lipid 128 heating 303K
 &cntrl
  imin=0,
  ntx=5,                        ! Positions and velocities read formatted
  irest=1,                      ! Restart calculation
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=50000,                 ! Number of MD steps
  ntt=3,
  gamma_ln=1.0, 
  ntr=1,
  restraint_wt=10               ! degree of restraint
  restraintmask=':OL,PC,PGR',   ! Restrain all except oleoyl lipid tails / acyl chains
  ig=-1,
  ntpr=100,
  ntwr=10000,
  ntwx=100,
  dt=0.002,
  nmropt=1,
  ntb=2,                        ! Constant pressure periodic boundary conditions
  ntp=2,                        ! Anisotropic pressure coupling
  taup=2.0,                     ! Pressure relaxation time (ps)
  cut=10.0,
  ioutfm=1,
  ntxo=2,
 /
 &wt
  type='TEMP0',
  istep1=0,
  istep2=50000,
  value1=100.0,
  value2=303.0 /
 &wt type='END' /
Hold lipids fixed
RES :OL,PC,PGR
END
END
EOF

date +"%T"
echo 'Heating 2'
pmemd.cuda -O -i heat2.in -o heat2mdout -p ${runname}.prmtop -c heat.rst -r heat2.rst -ref heat.rst -x heat2.nc

cat <<EOF > equil.in
Boundary dimensions equilibriation
 &cntrl
  imin=0,
  ntx=5,
  irest=1, 
  ntc=2,
  ntf=2,
  tol=0.0000001,
  nstlim=250000,
  ntt=3,
  gamma_ln=1.0,
  temp0=303.0,
  ntpr=5000,
  ntwr=5000,
  ntwx=5000,
  dt=0.002,
  ig=-1, 
  ntb=2,
  ntp=2,
  cut=10.0,
  ioutfm=1,
  ntxo=2,
 /
 /
 &ewald
  skinnb=5, ! Increase skinnb to avoid skinnb errors
 /
EOF

pmemd.cuda -O -i equil.in -o equil1mdout -p ${runname}.prmtop -c heat2.rst -r equil1.rst -x equil1.nc

# Repeat equilibriation
skinnb=5
j=2

while [ $j -lt 9 ] && [ $skinnb -lt 20 ] 
do
    date +"%T"
    echo 'Production run' $j
	pmemd.cuda -O -i equil.in -o equil"$j"mdout -p ${runname}.prmtop -c equil"$(($j-1))".rst -r equil"$j".rst -x equil"$j".nc
	if [ $? -eq 0 ]
	then
	 j=$(expr ${j} + 1)
	 exit 0
	else
	 grep -rl skinnb=${skinnb} equil.in | xargs sed -i 's/skinnb=${skinnb}/skinnb="$(($j+1))"/g'
	 skinnb=$(expr ${skinnb} + 1)
	fi
done

echo "Equilibriation finished with a skinnb value of $skinnb."
