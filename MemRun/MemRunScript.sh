#!/bin/bash
set -e

export AMBERHOME=/home/thomas/amber18
source $AMBERHOME/amber.sh

runname=GakA

packmol-memgen --pdb ${runname} --lipids DOPE:DOPG --ratio 3:1 --overwrite

cat << EOF > system.in
source leaprc.protein.ff19SB						
loadamberparams frcmod.ions1lm_126_iod_opc 				
source leaprc.water.opc
source leaprc.lipid17
system = loadpdb bilayer_${runname}
set system box { 100 100 100 }
saveamberparm system ${runname}.prmtop ${runname}.inpcrd
quit
EOF
tleap -f system.in

cat <<EOF > min.in
Minimization input file in explicit solvent
 &cntrl
    imin=1, maxcyc=5000,  ncyc=2500, 
    cut=9.0, ntpr=100,  ntxo=2, 
    ntr=1, nmropt=1,
    watnam='WAT',
    owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA.inpcrd
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
10.0
RES 1 34
END
Membrane posres
2.5
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
2.5
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil1.in
A NVT simulation for common production-level simulations
 &cntrl
    imin=0, irest=0, ntx=1,
    ntt=3, gamma_ln=1.0,
    tempi=303.15, temp0=303.15,
    cut=9.0, nstlim=125000, dt=0.001,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
!   ntwv=-1,     ! Uncomment to also print velocities to trajectory
!   ntwf=-1,      ! Uncomment to also print forces to trajectory
    ntxo=2, ioutfm=1, iwrap=0,
    ntr=1, nmropt=1,
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_min.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
10.0
RES 1 34
END
Membrane posres
2.5
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
2.5
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil2.in
A NVT simulation for common production-level simulations
 &cntrl
    imin=0, irest=1, ntx=5,
    ntt=3, gamma_ln=1.0,
    temp0=303.15,
    cut=9.0, nstlim=125000, dt=0.001,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
!   ntwv=-1,     ! Uncomment to also print velocities to trajectory
!   ntwf=-1,      ! Uncomment to also print forces to trajectory
    ntxo=2, ioutfm=1, iwrap=0,
    ntr=1, nmropt=1,
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_equil1.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
5.0
RES 1 34
END
Membrane posres
2.5
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
2.5
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil3.in
A NPT simulation for common production-level simulations
 &cntrl
    imin=0, irest=1, ntx=5,
    ntt=3, gamma_ln=1.0,
    temp0=303.15,
    cut=9.0, nstlim=125000, dt=0.001,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
!   ntwv=-1,     ! Uncomment to also print velocities to trajectory
!   ntwf=-1,      ! Uncomment to also print forces to trajectory
    ntxo=2, ioutfm=1, iwrap=0,
    barostat=1, ntp=3, pres0=1.0, taup=1.0
    csurften=3, gamma_ten=0.0, 
    ninterface=2, ntr=1, nmropt=1,
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_equil2.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
2.5
RES 1 34
END
Membrane posres
1
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
1
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil4.in
A NPT simulation for common production-level simulations
 &cntrl
    imin=0, irest=1, ntx=5,
    ntt=3, gamma_ln=1.0,
    temp0=303.15,
    cut=9.0, nstlim=125000, dt=0.002,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
!   ntwv=-1,     ! Uncomment to also print velocities to trajectory
!   ntwf=-1,      ! Uncomment to also print forces to trajectory
    ntxo=2, ioutfm=1, iwrap=0,
    barostat=1, ntp=3, pres0=1.0, taup=1.0
    csurften=3, gamma_ten=0.0, 
    ninterface=2, ntr=1, nmropt=1,
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_equil3.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
1
RES 1 34
END
Membrane posres
0.5
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
0.5
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil5.in
A NPT simulation for common production-level simulations
 &cntrl
    imin=0, irest=1, ntx=5,
    ntt=3, gamma_ln=1.0,
    temp0=303.15,
    cut=9.0, nstlim=250000, dt=0.002,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
!   ntwv=-1,     ! Uncomment to also print velocities to trajectory
!   ntwf=-1,      ! Uncomment to also print forces to trajectory
    ntxo=2, ioutfm=1, iwrap=0,
    barostat=1, ntp=3, pres0=1.0, taup=1.0
    csurften=3, gamma_ten=0.0, 
    ninterface=2, ntr=1, nmropt=1,
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_equil4.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
0.5
RES 1 34
END
Membrane posres
0.1
FIND
P31 * * PC
SEARCH
RES 1 9999
END
Membrane posres2
0.1
FIND
P31 * * OL
SEARCH
RES 1 9999
END
END
EOF

cat <<EOF > equil6.in
A NPT simulation for common production-level simulations
 &cntrl
    imin=0, irest=1, ntx=5,
    ntt=3, gamma_ln=1.0,
    temp0=303.15,
    cut=9.0, nstlim=250000, dt=0.002,
    ntc=2, ntf=2, ntpr=1000, ntwx=5000,
    ntwr=10000,
    ntxo=2, ioutfm=1, iwrap=0,
    barostat=1, ntp=3, pres0=1.0, taup=1.0
    csurften=3, gamma_ten=0.0, 
    ninterface=2, ntr=1, 
    watnam='WAT', owtnm='O',
 /
 &wt
    type='END'
 /
DISANG=GakA_equil5.xyz
LISTIN=POUT
LISTOUT=POUT
&end
Protein posres
0.1
RES 1 34
END
END
EOF

date +"%T"
echo 'Minimization'
mpirun -np 2 pmemd.MPI -O -i min.in -o minmdout -p ${runname}.prmtop -c ${runname}.inpcrd -r ${runname}_min.xyz -x ${runname}_min.nc

date +"%T"
echo 'Equilibriation'
pmemd.cuda -O -i equil1.in -o equil1mdout -p ${runname}.prmtop -c ${runname}_min.xyz -r ${runname}_equil1.xyz -x ${runname}_equil1.nc 

for j in {2..6}
do
    date +"%T"
    echo 'Production run' $j
	pmemd.cuda -O -i equil${j}.in -o equil${j}mdout -p ${runname}.prmtop -c ${runname}_equil"$(($j-1))".xyz -r ${runname}_equil"$j".xyz -x ${runname}_equil"$j".nc 
done

















