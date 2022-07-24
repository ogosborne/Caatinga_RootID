#!/bin/bash
BASEDIR=$(pwd)
mkdir ${BASEDIR}/STACKS_data/process_radtags
indir=seq_data/
outdir=STACKS_data/process_radtags/

for i in L-10 L-11 L-12 L-13 L-14 L-15 L-16 L-17 L-18 L-19 L-1 L-20 L-21 L-22 L-23 L-24 L-25 L-26 L-27 L-28 L-29 L-30 L-31 L-32 L-33 L-34 L-35 L-36 L-37 L-38 L-39 L-3 L-40 L-41 L-42 L-43 L-4 L-5 L-6 L-7 L-8 L-9 L-RTPa L-RTPb L-RTPc L-RTPd L-RTPe R-1A-05 R-1A-10 R-1A-20 R-1A-50 R-1B-05 R-1B-10 R-1B-20 R-1B-50 R-1C-05 R-1C-10 R-1C-20 R-1C-50 R-1D-05 R-1D-10 R-1D-20 R-1D-50 R-1E-05 R-1E-10 R-1E-20 R-1E-50 R-2A-05 R-2A-10 R-2A-20 R-2A-50 R-2B-05 R-2B-10 R-2B-20 R-2B-50 R-2C-05 R-2C-10 R-2C-20 R-2C-50 R-2D-05 R-2D-10 R-2D-20 R-2D-50 R-2E-05 R-2E-10 R-2E-20 R-2E-50 R-3A-05 R-3A-10 R-3A-20 R-3A-50 R-3B-05 R-3B-10 R-3B-20 R-3B-50 R-3C-05 R-3C-10 R-3C-20 R-3C-50 R-3D-05 R-3D-10 R-3D-20 R-3D-50 R-3E-05 R-3E-10 R-3E-20 R-3E-50 R-4A-05 R-4A-10 R-4A-20 R-4A-50 R-4B-05 R-4B-10 R-4B-20 R-4B-50 R-4C-05 R-4C-10 R-4C-20 R-4C-50 R-4D-05 R-4D-10 R-4D-20 R-4D-50 R-4E-05 R-4E-10 R-4E-20 R-4E-50 R-5A-05 R-5A-10 R-5A-20 R-5A-50 R-5B-05 R-5B-10 R-5B-20 R-5B-50 R-5C-05 R-5C-10 R-5C-20 R-5C-50 R-5D-05 R-5D-10 R-5D-20 R-5D-50 R-5E-05 R-5E-10 R-5E-20 R-5E-50 ; do
    process_radtags -1 ${BASEDIR}/${indir}${i}_R1.fastq.gz -2 ${BASEDIR}/${indir}${i}_R2.fastq.gz  --renz_1 pstI --renz_2 apeKI -q -E phred33 -s 25 -o ${BASEDIR}/${outdir} -y gzfastq -t 142 --len-limit 142 -c
done

