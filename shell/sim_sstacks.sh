# paths
usdir=SIM_data/ustacks/
csdir=SIM_data/cstacks/
outdir=SIM_data/sstacks/
mkdir $outdir
# iterators
nreads=$(seq 100000 100000 1000000)
# for each n reads
for n in $nreads ; do
    # set up outdir and temp dir
    mkdir ${outdir}/N${n} ; cd ${outdir}/N${n}
    mkdir temp ; cd temp
    # get ustacks data
    cp -l ${usdir}/*_${n}.* .
    # get catalogues
    cp -l ${csdir}/N${n}/*.tsv .
    # run stacks
    sstacks -p 32 -P . -M ${csdir}N${n}.popmap --disable-gapped
    # move tsv.gz and remove temp
    cd ..
    mv temp/*.matches.tsv .
    gzip *.matches.tsv
    rm -r temp
done