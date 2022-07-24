#!/bin/bash
BASEDIR=$(pwd)
cstacksdir=STACKS_data/cstacks/
ustacksdir=STACKS_data/ustacks/
mkdir ${BASEDIR}/STACKS_data/sstacks
sstacksdir=STACKS_data/sstacks/
popmap=data/popmap_


for M in $(seq 2 2 8) ; do
   # get n values
   minus1=$(echo "$M - 1" | bc)
   plus1=$(echo "$M + 1" | bc)
   # make outdir
   mkdir ${BASEDIR}/${sstacksdir}US_M${M}
   for n in ${minus1} ${M} ${plus1} ; do
      # make outdir
      mkdir ${BASEDIR}/${sstacksdir}US_M${M}/CS_n${n}
      for t in leaf root ; do
         # make outdir
         mkdir ${BASEDIR}/${sstacksdir}US_M${M}/CS_n${n}/${t}
         cd ${BASEDIR}/${sstacksdir}US_M${M}/CS_n${n}/${t}
         # set up temp dir
         mkdir temp ; cd temp
         cp -l ${BASEDIR}/${cstacksdir}US_M${M}/CS_n${n}/*.gz .
         cp -l ${BASEDIR}/${ustacksdir}US_M${M}/*.* .
         # run sstacks
         sstacks -p 32 -P ${BASEDIR}/${sstacksdir}US_M${M}/CS_n${n}/${t}/temp/ -M ${BASEDIR}/${popmap}${t}.txt --disable-gapped
         # move tsv.gz and remove temp
         cd ..
         mv temp/*.matches.tsv.gz .
         rm -r temp
         cd ${BASEDIR}
      done
   done
done



