1.

```
wait
wd=~/work_space/4.ProjectET/batch30_cr/1.rawdata
cd $wd
mkdir fastqc
for i in *fq.gz
do
nohup fastqc $i -o ./fastqc &
done

```

2.

```
conda activate py36
wd=~/work_space/4.ProjectET/batch30_cr/1.rawdata
nwd=~/work_space/4.ProjectET/batch30_cr/2.cutdata
mkdir $nwd
cd $wd
for i in *R1.fq.gz
do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 1 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 50 -q 20,20 -o $nwd/${i%.fq*}_cut.fq -p $nwd/${t%.fq*}_cut.fq $i $t > $nwd/${i%%.*}_cut.log &
done

```
```
wd=~/work_space/4.ProjectET/batch22_cr/2.cutdata
cd $wd
mkdir fastqc
for i in *fq
do
nohup fastqc $i -o ./fastqc &
done

```

3. align to mm10 & sacCer3

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

wd=~/work_space/4.ProjectET/batch30_cr/2.cutdata
nwd=~/work_space/4.ProjectET/batch30_cr/3.align
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
bowtie2 -p 10 -x /home1/share/bowtie2_index/mm10  --local --very-sensitive-local -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.*}.sam > $nwd/${i%%.*}.log 2>&1
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo




echo -e 'sample\ttotal_paired(million)\tproperly_paired(million)\tratio\tmore_1_time(million)\tratio\talign_ratio_all' >  bowtie2_align_summary.tab
for i in *.log; do
echo $i
key=${i%%.*}
awk -v k=$key 'BEGIN{FS="[ \t();%]+";OFS="\t"} NR==5 {t=$1} NR==8{a=$2;b=$3} NR==9{c=$2;d=$3} NR==10{f=$1} END{print k,t/1e6,a/1e6,b"%",c/1e6,d"%",f"%"}' $i >> bowtie2_align_summary.tab
done

```

4. sam2bam

```
wait
wd=~/work_space/4.ProjectET/batch30_cr/3.align
cd $wd
for i in *sam; do nohup samtools sort -@ 5 -o ${i%.*}.sorted.bam $i &
done

wait
wd=~/work_space/4.ProjectET/batch30_cr/3.align
nwd=~/work_space/4.ProjectET/batch30_cr/4.proper
mkdir -p $nwd
cd $wd
for i in *sorted.bam; do
nohup samtools view -bf 0x2 -@ 5 -q 20 $i -o $nwd/${i%sort*}proper.bam &
done


```
5. (optional) Remove dup


```
wait
wd=~/work_space/4.ProjectET/batch30_cr/4.proper
nwd=~/work_space/4.ProjectET/batch30_cr/5.unique
cd $wd
mkdir $nwd
for i in *.proper.bam;do
nohup sambamba markdup -r -t 4 -p $i $nwd/${i%proper*}duprm.bam > $nwd/${i%proper*}duprm.log &
done
```


macs

wait
conda activate work
wd=~/work_space/4.ProjectET/batch30_cr/5.unique
nwd=~/work_space/4.ProjectET/batch30_cr/6.peak
mkdir $nwd
cd $wd
for i in *bam
do
nohup macs2 callpeak -t $i  -f BAMPE -g mm -n ${i%%.*} -B -q 0.05 --keep-dup 1 --outdir $nwd > $nwd/${i%%.*}.log &
done

7. Peak anno


wait
wd=~/work_space/4.ProjectET/batch30_cr/6.peak
nwd=~/work_space/4.ProjectET/batch30_cr/7.Peak_anno
mkdir $nwd
cd $wd
for i in *narrowPeak
do
nohup findMotifsGenome.pl $i mm10 $nwd/${i%_*}_peakAnalysis -size given -p 7  -len 8 > $nwd/${i%_*}_peakAnalysis.log &
done

a. SNP align
```
#!/bin/sh
THREAD=4
tempfifo=$$.fifo
trap 'exec 1000>&-;exec 1000<&+;rm -rf $tempfifo;exit 0' 2
mkfifo $tempfifo
exec 1000<>$tempfifo
for ((i=0;i<$THREAD;i++))
do
echo >&1000
done

wd=~/work_space/4.ProjectET/batch28_cr/2.cutdata
nwd=~/work_space/4.ProjectET/batch28_cr/a.snp_align
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
bowtie2 -p 8 -x /home1/share/snpsplit/PWK_Phj_single_strain/PWK_Phj_bowtie2_index/PWK_Phj_N_masked  -t --no-mixed --no-discordant --no-unal -1 $i -2 $t -S $nwd/${i%%.R1*}.sam > $nwd/${i%%.R1*}.log 2>&1
echo >&1000
}&
done
wait
echo "All Done!"
rm -rf $tempfifo
```

nohup sh /home1/gyang/work_space/4.ProjectET/batch28_cr/bowtie_align_snp.sh &

wd=~/work_space/4.ProjectET/batch28_cr/a.snp_align
cd $wd
for i in *sam
do
nohup samtools sort -@ 4 -o ${i%.*}.sorted.bam $i &
done

wd=~/work_space/4.ProjectET/batch14_cr/a.snp_align
nwd=~/work_space/4.ProjectET/batch14_cr/b.snp_proper
mkdir -p $nwd
cd $wd
for i in *sorted.bam
do
nohup samtools view -bf 0x2 -@ 3 -q 20 $i -o $nwd/${i%sort*}proper.bam &
done

wd=~/work_space/4.ProjectET/batch14_cr/b.snp_proper
nwd=~/work_space/4.ProjectET/batch14_cr/c.snp_unique
cd $wd
mkdir -p $nwd
for i in *.bam;do
nohup sambamba markdup -r -t 4 -p $i $nwd/${i%proper*}duprm.bam > $nwd/${i%proper*}duprm.log &
done

wd=~/work_space/4.ProjectET/batch14_cr/c.snp_unique
nwd=~/work_space/4.ProjectET/batch14_cr/d.snp_split
mkdir -p $nwd
cd $wd
for i in *duprm.bam
do
nohup SNPsplit  --snp_file /home1/share/snpsplit/PWK_Phj_single_strain/all_SNPs_PWK_PhJ_GRCm38.txt.gz $i --paired -o $nwd/${i%.duprm*}  --conflicting  > ${i%.duprm*}_snpsplit.log &
done

wd=~/work_space/4.ProjectET/batch14_cr/d.snp_split
nwd=~/work_space/4.ProjectET/batch14_cr/e.snp_split_sort
cd $wd
mkdir $nwd
for i in *
do
cd $i
for j in *genome*bam
do
nohup samtools sort -@ 3 -o $nwd/${j%.*}.sorted.bam $j &
done
cd ..
done

wd=~/work_space/4.ProjectET/batch14_cr/e.snp_split_sort
cd $wd
for i in *bam
do
nohup samtools index $i &
done

conda activate py36
wd=~/work_space/4.ProjectET/batch16_cr/e.snp_split_sort
nwd=~/work_space/4.ProjectET/batch16_cr/f.snp_bw
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

5.0

```
source activate py36
export TMPDIR=/home1/gyang/work_space/4.ProjectET/batch3_cr/tmp
mkdir $TMPDIR
wd=~/work_space/4.ProjectET/batch1_cr/5.unique
nwd=~/work_space/4.ProjectET/batch1_cr/5.1len_40_120
mkdir -p $nwd
cd $wd
for i in *.bam
do
name=${i%%.*}
echo $name
nohup alignmentSieve -p 2 -b $i -o $nwd/${name}.sieved.bam \
--minFragmentLength 40  --maxFragmentLength 120 --filterMetrics $nwd/${i%%.*}_sieve.log &
done

source activate work
wd=~/work_space/4.ProjectET/batch3_cr/5.1len_40_120
cd $wd
for i in *bam
do
nohup samtools index $i &
done

```
5.1 bw

```
wait
conda activate py36
wd=~/work_space/4.ProjectET/batch30_cr/5.unique
nwd=~/work_space/4.ProjectET/batch30_cr/z.bw
mkdir -p $nwd
cd $wd
for i in *.bam
do
name=${i%%.*}
echo $name
nohup bamCoverage  -p 3 -b $i -o $nwd/${name}.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > $nwd/${name}.bw.log &
done
```

 sample qc Correlation

```
wd=/home1/gyang/work_space/4.ProjectET/batch5_cr/z.bw
nwd=/home1/gyang/work_space/4.ProjectET/batch5_cr/sample_qc_and_comparison
mkdir $nwd
cd $wd
nohup multiBigwigSummary BED-file -b  `ls *bw`  --outRawCounts $nwd/all_batch5_promoter.tab -o $nwd/temp.npz --BED $rwd/5most_promoter_mm10.tab -p 48 &
```














A. PCA

```
conda activate py36
wd=/home1/gyang/work_space/4.ProjectET/batch19_cr/z.bw
nwd=/home1/gyang/work_space/4.ProjectET/batch19_cr/z.pca
mkdir $nwd
cd $wd
multiBigwigSummary bins  -b *bw  \
 -out $nwd/all-batch19.npz --outRawCounts $nwd/all-batch19.tab --smartLabels \
--numberOfProcessors 48 &

cd $nwd
nohup plotPCA -in all-batch19.npz --transpose --outFileNameData all-batch19.tab -T all-batch19_PCA  -o all-batch19.pdf --plotHeight 5 --plotWidth 10 &

nohup plotCorrelation -in batch6.npz  --whatToPlot heatmap -c pearson --removeOutliers --plotNumbers -o batch6_cor.png &
```
