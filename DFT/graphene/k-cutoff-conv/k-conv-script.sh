#!/bin/sh
NAME="k-cut-graphene"

for KVAL in  06 09 12 18 24
do
cat > ${NAME}_${KVAL}_.in << EOF
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
    celldm(1) = 4.654,
    celldm(3) = 3.0,
    nat  = 2,
    ntyp = 1,
    ecutwfc = 40,
    ecutrho = 400,
    occupations = 'smearing',
    smearing = 'gaussian',
    degauss = 0.01
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
   $KVAL $KVAL 1   0 0 0
EOF

pw.x < ${NAME}_${KVAL}_.in > ${NAME}_${KVAL}_.out
echo ${NAME}_${KVAL}
grep ! ${NAME}_${KVAL}_.out

done
