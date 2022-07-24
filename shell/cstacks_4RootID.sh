#!/bin/bash
BASEDIR=$(pwd)
indir=STACKS_data/ustacks/
mkdir ${BASEDIR}/STACKS_data/cstacks
catdir=STACKS_data/cstacks/
popmap=data/popmap_leaf.txt

for M in $(seq 2 2 8) ; do
   # get n values
   minus1=$(echo "$M - 1" | bc)
   plus1=$(echo "$M + 1" | bc)
   # make outdir
   mkdir ${BASEDIR}/${catdir}US_M${M}
   for n in ${minus1} ${M} ${plus1} ; do
   	  # set up temp dir 
   	  mkdir ${BASEDIR}/${catdir}US_M${M}/CS_n${n}
      cd ${BASEDIR}/${catdir}US_M${M}/CS_n${n}
      mkdir temp ; cd temp
      cp -l ${BASEDIR}/${indir}US_M${M}/*.* .
      # run cstacks
      cstacks -n ${n} --disable-gapped -p 32 -P {BASEDIR}/${catdir}US_M${M}/CS_n${n}/temp/ -M ${BASEDIR}/${popmap}
      # keep cats, remove rest
      mv catalog.* .. ; cd ..
      rm -r temp
      cd ${BASEDIR}
   done
done
