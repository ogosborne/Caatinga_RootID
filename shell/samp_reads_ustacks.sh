# paths
indir=SIM_data/sim/
outdir=SIM_data/subsamp/
usdir=SIM_data/ustacks/
mkdir $outdir
mkdir $usdir
# iterators
nreads=$(seq 100000 100000 1000000)
spp=(Ao  At  Eg  Ha  Ls  Pt)
# ustacks options
usM=6
usm=5
usp=32
# loop
for sp in ${spp[*]} ; do
    # for each n reads
    for n in $nreads ; do
    	# log for random seed
    	echo "species individual n_reads seed" > ${sp}_${n}.log
    	rand=$(echo $RANDOM)
    	for ind in $(seq 0 1 9) ; do
            # write seed to log
            echo "$sp $ind $n $rand" >> ${sp}_${n}.log
            # subsample
            $seqtk sample -s $rand ${indir}${sp}_${ind}.1.fa.gz $n > ${outdir}${sp}_${ind}_${n}.1.fa
            $seqtk sample -s $rand ${indir}${sp}_${ind}.2.fa.gz $n > ${outdir}${sp}_${ind}_${n}.2.fa
            # concatenate
            paste -d "" ${outdir}${sp}_${ind}_${n}.1.fa ${outdir}${sp}_${ind}_${n}.2.fa > ${outdir}${sp}_${ind}_${n}.fa
            # run ustacks
            ustacks -t fasta -f ${outdir}${sp}_${ind}_${n}.fa -i $rand -o ${usdir} -p ${usp} -d -m ${usm} -M ${usM} -H --disable-gapped
        done
    done
done
