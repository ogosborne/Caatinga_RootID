#!/bin/bash
BASEDIR=$(pwd)
indir=STACKS_data/process_radtags/
mkdir ${BASEDIR}/STACKS_data/ustacks
outdir=STACKS_data/ustacks/US_M
# leaf samples
for M in $(seq 2 2 8) ; do
   mkdir ${outdir}${M}
   n=1
   for i in L_10 L_11 L_1 L_12 L_13 L_14 L_15 L_16 L_17 L_18 L_19 L_20 L_21 L_22 L_23 L_24 L_25 L_26 L_27 L_28 L_29 L_30 L_31 L_3 L_32 L_33 L_34 L_35 L_36 L_37 L_38 L_39 L_40 L_41 L_4 L_42 L_43 L_5 L_6 L_7 L_8 L_9 L_RTPa L_RTPb L_RTPc L_RTPd L_RTPe ; do
      ustacks -t gzfastq -f ${BASEDIR}/${indir}${i}.cat.fq -i ${n} -o ${BASEDIR}/${outdir}${M} -p 32 -d -m 5 -M ${M} -H --disable-gapped
      n=$((n+1))
  done
done
# root samples
for M in 2 4 6 8 ; do
   n=49
   for i in R_1A_05 R_1A_10 R_1A_20 R_1A_50 R_1B_05 R_1B_10 R_1B_20 R_1B_50 R_1C_05 R_1C_10 R_1C_20 R_1C_50 R_1D_05 R_1D_10 R_1D_20 R_1D_50 R_1E_05 R_1E_10 R_1E_20 R_1E_50 R_2A_05 R_2A_10 R_2A_20 R_2A_50 R_2B_05 R_2B_10 R_2B_20 R_2B_50 R_2C_05 R_2C_10 R_2C_20 R_2C_50 R_2D_05 R_2D_10 R_2D_20 R_2D_50 R_2E_05 R_2E_10 R_2E_20 R_2E_50 R_3A_05 R_3A_10 R_3A_20 R_3A_50 R_3B_05 R_3B_10 R_3B_20 R_3B_50 R_3C_05 R_3C_10 R_3C_20 R_3C_50 R_3D_05 R_3D_10 R_3D_20 R_3D_50 R_3E_05 R_3E_10 R_3E_20 R_3E_50 R_4A_05 R_4A_10 R_4A_20 R_4A_50 R_4B_05 R_4B_10 R_4B_20 R_4B_50 R_4C_05 R_4C_10 R_4C_20 R_4C_50 R_4D_05 R_4D_10 R_4D_20 R_4D_50 R_4E_05 R_4E_10 R_4E_20 R_4E_50 R_5A_05 R_5A_10 R_5A_20 R_5A_50 R_5B_05 R_5B_10 R_5B_20 R_5B_50 R_5C_05 R_5C_10 R_5C_20 R_5C_50 R_5D_05 R_5D_10 R_5D_20 R_5D_50 R_5E_05 R_5E_10 R_5E_20 R_5E_50 ; do
      ustacks -t gzfastq -f ${BASEDIR}/${indir}${i}.cat.fq -i ${n} -o ${BASEDIR}/${outdir}${M} -p 32 -d -m 1 -M ${M} -H --disable-gapped
      n=$((n+1))
  done
done


