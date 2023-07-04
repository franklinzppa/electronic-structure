#!/bin/sh
NAME="a-graphene"

for AVAL in  4.55 4.57 4.59 4.61 4.63 4.65 4.67 4.69 4.71 4.73 4.75
do
cat > ${NAME}_${AVAL}_.in << EOF
 &control
    calculation = 'scf',
    restart_mode = 'from_scratch',
    prefix      = 'graphene',
    outdir      = './tmp/graphene',
    pseudo_dir  = './',
    verbosity   = 'high'
 /
 &system
    ibrav     = 4,
    celldm(1) = $AVAL,
    celldm(3) = 3.0,
    nat  = 2,
    ntyp = 1,
    ecutwfc = 40,
    ecutrho = 400,
    ! occupations = 'smearing',
    ! smearing = 'gaussian',
    ! degauss = 0.01
 /
 &electrons
    conv_thr = 1.0d-8
    mixing_beta = 0.7
 /

ATOMIC_SPECIES
   C  12.0107 C.pz-n-rrkjus_psl.0.1.UPF
   
ATOMIC_POSITIONS alat
   C    0.000000    0.0000000   0.000000
   C    0.000000    0.5773503   0.000000
   
K_POINTS automatic
   15 15 1  0 0 0
EOF

pw.x < ${NAME}_${AVAL}_.in > ${NAME}_${AVAL}_.out
echo ${NAME}_${AVAL}
grep ! ${NAME}_${AVAL}_.out

done
