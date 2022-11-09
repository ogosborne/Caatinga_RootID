# First download and unzip the genomes from phytozome
# make chromosome lists
for i in $(ls *.fa | sed 's/.fa//g' ) ; do 
   grep "^>" ${i}.fa | sed 's/>//g' | cut -f1 -d" " > ${i}.chrom.list
done
# replace ambiguity chracters with N
for i in $(ls *.fa | sed 's/.fa//g' ) ; do 
   sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D]/N/g' ${i}.fa > ${i}.faN 
   mv ${i}.faN ${i}.fa
done
# simulate populations
mkdir SIM_data
popdir=SIM_data/populations/
mkdir $popdir

###### options
Ne=1000

###### species specific options
sp=Aofficinalis_498_Aspof.V1
sp_short=Ao
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10

###### species specific options
sp=Athaliana_447_TAIR10
sp_short=At
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10

###### species specific options
sp=Egrandis_297_v2.0
sp_short=Eg
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10

###### species specific options
sp=Hannuus_494_r1.0
sp_short=Ha
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10

###### species specific options
sp=Lsativa_467_v8
sp_short=Ls
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10

###### species specific options
sp=Ptrichocarpa_533_v4.0
sp_short=Pt
######
mkdir ${popdir}${sp_short}
radinitio --make-population \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${popdir}${sp_short} \
--n-pops 1 --pop-eff-size ${Ne} --n-seq-indv 10
