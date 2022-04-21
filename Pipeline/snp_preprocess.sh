1. cut

conda activate py36
wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/1.rawdata
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/2.cutdata
mkdir $nwd
cd $wd
for i in *R1.fq.gz
do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 2 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 50 -q 20,20 -o $nwd/${i%.fq*}_cut.fq -p $nwd/${t%.fq*}_cut.fq $i $t > $nwd/${i%%.*}_cut.log &
done

2. SNP align
```
#!/bin/sh
THREAD=6
tempfifo=$$.fifo
trap 'exec 1000>&-;exec 1000<&+;rm -rf $tempfifo;exit 0' 2
mkfifo $tempfifo
exec 1000<>$tempfifo
for ((i=0;i<$THREAD;i++))
do
echo >&1000
done

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/2.cutdata
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/2.cutdata/3.snp_align
mkdir -p $nwd
cd $wd
file=`ls *R1_cut.fq`
for i in $file
do
t=${i/R1/R2}
echo $i
echo $t
read -u 1000
{
bowtie2 -p 7 -x /home1/share/snpsplit/PWK_Phj_single_strain/PWK_Phj_bowtie2_index/PWK_Phj_N_masked  -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.R1*}.sam > $nwd/${i%%.R1*}.log 2>&1
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo
```

nohup sh ~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/bowtie_snp_align.sh &

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/3.snp_align
cd $wd
for i in *sam
do
nohup samtools sort -@ 2 -o ${i%.*}.sorted.bam $i &
done

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/3.snp_align
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/4.snp_proper
mkdir -p $nwd
cd $wd
for i in *sorted.bam
do
nohup samtools view -bf 0x2 -@ 2 -q 20 $i -o $nwd/${i%sort*}proper.bam &
done

wait
wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/4.snp_proper
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/5.snp_unique
cd $wd
mkdir -p $nwd
for i in *.bam;do
nohup sambamba markdup -r -t 2 -p $i $nwd/${i%proper*}duprm.bam > $nwd/${i%proper*}duprm.log &
done

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/5.snp_unique
nwd=$wd/../6.snp_merged
mkdir $nwd
cd $wd
for i in *rep1.duprm.bam
do
name=${i/_rep1.duprm.bam/}
ls ${name}*bam
nohup samtools merge $nwd/${name}.bam `ls ${name}*bam` -@ 2 &
done

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/6.snp_merged
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/7.snp_split
mkdir -p $nwd
cd $wd
for i in *bam
do
nohup SNPsplit  --snp_file /home1/share/snpsplit/PWK_Phj_single_strain/all_SNPs_PWK_PhJ_GRCm38.txt.gz $i --paired -o $nwd/${i%.*}  --conflicting  > ${i%.*}_snpsplit.log &
done

wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/7.snp_split
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort
cd $wd
mkdir $nwd
for i in *
do
cd $i
for j in *genome*bam
do
nohup samtools sort -@ 2 -o $nwd/${j%.*}.sorted.bam $j &
done
cd ..
done

wait
wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort
cd $wd
for i in *bam
do
nohup samtools index $i &
done

wait
conda activate py36
wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort
nwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/9.snp_bw
mkdir $nwd
cd $wd
for i in *genome1*bam
do
j=${i/genome1/genome2}
bamCoverage  -p 3 -b $i -o $nwd/${i%%.*}_c57.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > $nwd/${i%%.*}_c57.bw.log &
bamCoverage  -p 3 -b $j -o $nwd/${j%%.*}_pwk.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > $nwd/${j%%.*}_pwk.bw.log &
done

3. SNP count

######### create a bam version without ChrXY
wd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort/Tead4
cd $wd
mkdir without_chrxy
for i in *.bam
do
nohup samtools view -@ 2 -o ./without_chrxy/${i/.bam/_chr1-19.bam} $i `seq 1 19 | sed 's/^/chr/'` &
done

wait
cd $wd/without_chrxy
for i in *bam
do
samtools index -@ 2 $i &
done

conda activate py36
bwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort/Tead4/without_chrxy
wd=/home1/gyang/work_space/4.ProjectET/analysis/figs3
cd $bwd
nohup multiBamSummary bins --bamfiles `ls *bam` -e  -p 20 -bs 1000 --smartLabels -o $wd/temp.npz --outRawCounts $wd/project_et_t2_allele.counts &

bwd=~/work_space/4.ProjectET/analysis/standardize/1.snp_raw_fastq/8.snp_split_sort/Tfap2c/without_chrxy
wd=/home1/gyang/work_space/4.ProjectET/analysis/figs3
cd $bwd
nohup multiBamSummary bins --bamfiles `ls *bam` -e  -p 20 -bs 1000 --smartLabels -o $wd/temp.npz --outRawCounts $wd/project_et_t1_allele.counts &


sed -i "s/'//g" *counts
sed -i "s/#//g" *counts
sed -i "s/.sorted_chr1-19//g" *counts

4. Extract ASP motif loci

cd ~/work_space/4.ProjectET/analysis/figs3/AS_stat
mwd=/home1/gyang/work_space/4.ProjectET/analysis/motif_denovo
for i in *Tfap2c*.peak.bed
do
findMotifsGenome.pl $i mm10 ABP_motif -size given  -p 20 -find $mwd/Tfap2c.motif 1>ABP_motif/${i/peak.bed/peak.motif.bed}  2>ABP_motif/${i/peak.bed/peak.motif.bed}.log
done

mwd=/home1/gyang/work_space/4.ProjectET/analysis/motif_denovo
for i in *Tead4*.peak.bed
do
findMotifsGenome.pl $i mm10 ABP_motif -size given  -p 3 -find $mwd/Tead4.motif 1>ABP_motif/${i/peak.bed/peak.motif.bed}  2>ABP_motif/${i/peak.bed/peak.motif.bed}.log
done

conda activate py36
cd ~/work_space/4.ProjectET/analysis/figs3/AS_stat/ABP_motif
for i in *peak.motif.bed
do
echo $i
Rscript ~/software/call_abs_homer_motif_locus.R ../${i/peak.motif.bed/peak.bed} $i ${i/peak.motif.bed/peak.motif.abs.bed}
done

------- 20220308

cd /home1/gyang/work_space/4.ProjectET/analysis/figs3/tf_snp/2.ABP_related_genes_10k
mwd=/home1/gyang/work_space/4.ProjectET/analysis/motif_denovo
mkdir ABP_motif
findMotifsGenome.pl allstage_ABP10k.bed mm10 ABP_motif -size given  -p 20 -find $mwd/Tfap2c.motif 1>ABP_motif/allstage_Tfap2c_ABP10k.motif.bed  2>ABP_motif/allstage_Tfap2c_ABP10k.motif.bed.log

conda activate py36
Rscript ~/software/call_abs_homer_motif_locus.R ../allstage_Tfap2c_ABP10k.bed allstage_Tfap2c_ABP10k.motif.bed allstage_Tfap2c_ABP10k.abs.bed

cd /home1/gyang/work_space/4.ProjectET/analysis/figs3/tf_snp/2.ABP_related_genes_10k
mwd=/home1/gyang/work_space/4.ProjectET/analysis/motif_denovo
mkdir ABP_motif
findMotifsGenome.pl allstage_Tead4_ABP10k.bed mm10 ABP_motif -size given  -p 20 -find $mwd/Tead4.motif 1>ABP_motif/allstage_Tead4_ABP10k.motif.bed  2>ABP_motif/allstage_Tead4_ABP10k.motif.bed.log

conda activate py36
Rscript ~/software/call_abs_homer_motif_locus.R ../allstage_Tead4_ABP10k.bed allstage_Tead4_ABP10k.motif.bed allstage_Tead4_ABP10k.abs.bed

 conda activate py36
 bwd=~/work_space/4.ProjectET/analysis/figs3/build/epi_resource/k4me3
 pwd=/home1/gyang/work_space/4.ProjectET/analysis/figs3/tf_snp/3.ABP10k_ABW_motif_genomeanno
 wd=/home1/gyang/work_space/4.ProjectET/analysis/figs3/tf_snp/3.ABP10k_ABW_motif_genomeanno
 cd $bwd
 computeMatrix reference-point -S *.bw -R ${pwd}/ABP10k_ABW_motif.bed  --skipZeros\
                              --beforeRegionStartLength 1500 \
                            --afterRegionStartLength 1500  -o $wd/AAM_K4me3.mat.gz  --binSize 50 -p 40
cd $wd
plotProfile -m AAM_K4me3.mat.gz  -out AAM_K4me3.pdf

------- 20220318 Human
cd ~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
mwd=/home1/homer/motifs
mkdir AP2_motif
findMotifsGenome.pl NCB2018_lftover_hg38.bed hg38 AP2_motif -size given  -p 20 -find $mwd/ap2gamma.motif 1>AP2_motif/NCB2018_lftover_hg38.motif.bed  2>AP2_motif/NCB2018_lftover_hg38.motif.bed.log

conda activate py36
Rscript ~/software/call_abs_homer_motif_locus.R ../NCB2018_lftover_hg38.bed NCB2018_lftover_hg38.motif.bed  NCB2018_lftover_hg38.abs.bed

cd ~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
bedtools intersect -a NCB2018_lftover_hg38.abs.bed -b ~/work_space/Bio-Resource/eqtl/all.bed > ../AP2_motif_eqtl.bed

 wd=~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
 cd $wd
 mkdir PeakAnno
 for i in *bed
 do
annotatePeaks.pl $i hg38  1>PeakAnno/${i%_peak*}.peakanno  2>&1 &
done

conda activate py36
bwd=~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
pwd=~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
wd=~/work_space/4.ProjectET/analysis/figs3/tf_snp/4.hg38_TFAP2C
cd $bwd
computeMatrix reference-point -S *.bw -R ${pwd}/AP2_motif_eqtl_inter.bed --beforeRegionStartLength 2000 \
 --skipZeros   --afterRegionStartLength 2000 \
  -o $wd/AP2_ATAC_eqtl.mat.gz  --binSize 50 -p 10

cd $wd
plotProfile -m AP2_ATAC_eqtl.mat.gz   -out AP2_ATAC_eqtl2.pdf  --startLabel motif  --perGroup


Backcross_RNA_Scattergroup_anno_redname.csv


~/work_space/1.Mouse_Acetylation/buffet/resource_set/all_mm10_5most_genes_tss_tes_genename.tab

awk 'BEGIN{OFS="\t"}NR==FNR{ a[$0] } NR>FNR{if($5 in a){print $1,$2,$3}}' Backcross_RNA_Scattergroup_anno_redname.csv ~/work_space/1.Mouse_Acetylation/buffet/resource_set/all_mm10_5most_genes_tss_tes_genename.tab  | sort -k1V -k2n > red_genebody.bed
awk 'BEGIN{OFS="\t"}NR==FNR{ a[$0] } NR>FNR{if($5 in a){print $1,$2,$3}}' Backcross_RNA_Scattergroup_anno_bluename.csv ~/work_space/1.Mouse_Acetylation/buffet/resource_set/all_mm10_5most_genes_tss_tes_genename.tab | sort -k1V -k2n > blue_genebody.bed
awk 'BEGIN{OFS="\t"}NR==FNR{ a[$0] } NR>FNR{if($5 in a){print $1,$2,$3}}' Backcross_RNA_Scattergroup_anno_greyname.csv ~/work_space/1.Mouse_Acetylation/buffet/resource_set/all_mm10_5most_genes_tss_tes_genename.tab | sort -k1V -k2n > grey_genebody.bed

j figs3
for  i in *body.bed
do
nohup bedtools intersect -a $i -b ~/work_space/4.ProjectET/analysis/fig3/build/loop/8Cloop_raw.bedpe -wa -wb| cut -f 4,5,6 > ${i}.enhancer.bed &
done

for i in *enhancer.bed
do
nohup bedtools intersect -a $i -b /home1/share/snpsplit/PWK_Phj_single_strain/pwk_snp.bed -wa -wb | bedtools groupby -c 4 -o count > ${i}.pwk.snpcount &
done

for i in *enhancer.bed
do
nohup bedtools intersect -a $i -b /home1/share/snpsplit/129S1_Sv_single_strain/S129_snp.bed -wa -wb | bedtools groupby -c 4 -o count > ${i}.S129.snpcount &
done

for i in *enhancer.bed
do
nohup bedtools intersect -a $i -b /home1/share/snpsplit/CAST_EiJ_single_strain/CAST_snp.bed -wa -wb | bedtools groupby -c 4 -o count > ${i}.CAST.snpcount &
done

# j figs3
# for i in *body.bed
# do
# nohup bedtools intersect -a $i -b /home1/share/snpsplit/PWK_Phj_single_strain/pwk_snp.bed -wa -wb | bedtools groupby -c 4 -o count > ${i}.snpcount &
# done
