popdir=SIM_data/populations/
outdir=SIM_data/sim/
mkdir $outdir

#general options
npops=1
nind=10
enz1=pstI
#enz2 apeKI isn't available, use default
readlen=150
libtype=ddRAD
pcrcyc=10

###### species specific options
cov=19
sp=Aofficinalis_498_Aspof.V1
sp_short=Ao
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}

###### species specific options
cov=222
sp=Athaliana_447_TAIR10
sp_short=At
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}

###### species specific options
cov=37
sp=Egrandis_297_v2.0
sp_short=Eg
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}

###### species specific options
cov=9
sp=Hannuus_494_r1.0
sp_short=Ha
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}

###### species specific options
cov=15
sp=Lsativa_467_v8
sp_short=Ls
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}

###### species specific options
cov=99
sp=Ptrichocarpa_533_v4.0
sp_short=Pt
######
mkdir ${outdir}${sp_short}
radinitio --make-library-seq \
--genome ${sp}.fa \
--chromosomes ${sp}.chrom.list \
--out-dir ${outdir}${sp_short} \
--make-pop-sim-dir ${popdir}${sp_short} \
--library-type ${libtype} \
--pcr-cycles ${pcrcyc} \
--enz ${enz1} \
--coverage ${cov} \
--read-length ${readlen}