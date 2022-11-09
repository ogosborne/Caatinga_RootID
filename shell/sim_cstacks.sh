# paths
usdir=SIM_data/ustacks/
csdir=SIM_data/cstacks/
mkdir $csdir
# iterators
nreads=$(seq 100000 100000 1000000)
spp=(Ao  At  Eg  Ha  Ls  Pt)
# cstacks options
csn=7
csp=32
# make popmaps
# for each n reads
for n in $nreads ; do
	# initiate file
	echo -n > ${csdir}N${n}.popmap
	# for each species
    for sp in ${spp[*]} ; do
        # for each individual
        for ind in $(seq 0 1 9) ; do
            echo -e "${sp}_${ind}_${n}\t1" >> ${csdir}N${n}.popmap
        done
    done
done
# run cstacks
for n in $nreads ; do
    mkdir ${csdir}/N${n}
    cd ${csdir}/N${n}
    mkdir temp; cd temp
    cp -l ${usdir}/*_${n}.* .
    cstacks -n ${csn} --disable-gapped -p ${csp} -P . -M ${csdir}N${n}.popmap
    mv catalog.* .. ; cd ..
    rm -r temp
done
