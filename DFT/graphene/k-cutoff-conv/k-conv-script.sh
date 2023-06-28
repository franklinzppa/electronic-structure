#!/bin/sh
NAME="k-cut-graphene"

for K in 04 06 08 10 12 16 
do
cat > ${NAME}_${K}.in << EOF
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
    ecutwfc = 50,
    ecutrho = 500,
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
   $K $K 1  0 0 0
EOF

pw.x < ${NAME}_${K}.in > ${NAME}_${K}.out
echo ${NAME}_${K}
grep ! ${NAME}_${K}.out

done