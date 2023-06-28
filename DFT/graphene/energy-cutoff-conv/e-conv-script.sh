#!/bin/sh
NAME="e-cut-graphene"

for CUTOFF in  10 15 20 25 30 35 40 45 50 60
do
cat > ${NAME}_${CUTOFF}.in << EOF
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
    ecutwfc = $CUTOFF,
    ecutrho = $((CUTOFF * 10)),
    nbnd = 8
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
   4 4 1   0 0 0
EOF

pw.x < ${NAME}_${CUTOFF}.in > ${NAME}_${CUTOFF}.out
echo ${NAME}_${CUTOFF}
grep ! ${NAME}_${CUTOFF}.out

done
