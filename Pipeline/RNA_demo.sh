Fastqc

wd=~/work_space/4.ProjectET/batch34_rna/1.rawdata
cd $wd
mkdir fastqc
for i in *fq.gz ;do
nohup fastqc $i -o ./fastqc &
done

multiqc -n "batch17_rna_raw" ./ &

Cutadapt

conda activate py36
wd=~/work_space/4.ProjectET/batch34_rna/1.rawdata
nwd=~/work_space/4.ProjectET/batch34_rna/2.cutdata
mkdir $nwd
cd $wd
for i in *R1.fq.gz
do
t=${i/R1/R2}
echo $i
echo $t
nohup cutadapt -j 3 -a AGATCGGAAGAGC  -A AGATCGGAAGAGC --trim-n -m 80 -q 20,20 -o $nwd/${i%.fq*}_cut.fq -p $nwd/${t%.fq*}_cut.fq $i $t > $nwd/${i%%.*}_cut.log &
done


wd=/home1/gyang/work_space/Public_data/7.zhangyi_nt_rna/2.cutdata
cd $wd
mkdir fastqc
for i in *fq ;do nohup fastqc $i -o ./fastqc & done


3. align
wait
wd=~/work_space/4.ProjectET/batch34_rna/2.cutdata
nwd=~/work_space/4.ProjectET/batch34_rna/3.align
mkdir $nwd
cd $wd
for i in *.R1_cut.fq
do
echo $i
t=${i/R1_cut.fq/R2_cut.fq}
echo $t
nohup hisat2 -p 10 --dta-cufflinks --no-discordant  -t -x /home1/share/hisat2_index/mm10/genome -1 $i -2 $t -S $nwd/${i%%.*}.sam  > $nwd/${i%%.*}.log &
echo ${i%%.*}.sam
done


4.sam2bam
wait
wd=~/work_space/4.ProjectET/batch34_rna/3.align
cd $wd
mkdir flagstat
for i in *.sam
do
echo $i
nohup samtools flagstat -@ 1 $i > ./flagstat/${i/.sam/.flagstat} &
done

wait
for i in *sam
do
nohup samtools sort -@ 5 -o ${i%.*}.sorted.bam $i &
done

wait
wd=~/work_space/4.ProjectET/batch34_rna/3.align
nwd=~/work_space/4.ProjectET/batch34_rna/4.proper
mkdir $nwd
cd $wd
for i in *sorted.bam
do
nohup samtools view -bf 0x2 -q 20 -@ 4 $i -o $nwd/${i%sort*}proper.bam  &
done

wait
cd $nwd
for i in *bam
do
nohup samtools index -@ 2 $i &
done

4.1 merge rep
wd=~/work_space/1.Mouse_Acetylation/batch34_rna/4.proper
cd $wd
for i in *rep1.proper.bam
do
echo ${i%_rep1*}
j=${i%_rep1*}
nohup samtools merge ${j}_merged.bam `ls *$j*` -@ 5 &
done

5.stringtie
wait
wd=~/work_space/4.ProjectET/batch34_rna/4.proper
nwd=~/work_space/4.ProjectET/batch34_rna/5.stringtie
cd $wd
for i in *proper.bam; do
   outfolder=${i%.proper*}
   mkdir -p $nwd/$outfolder
   nohup stringtie -p 4 $i -b $nwd/$outfolder -e -G /home1/share/gtf/mm10.gtf -o $nwd/$outfolder/$outfolder_stringtie.gtf > $nwd/$outfolder.log &
done

6.

featureCounts_dir=~/work_space/4.ProjectET/batch17_rna/5.featurecounts
ann_dir=/home1/share/gtf
hisat2_bam_dir=~/work_space/4.ProjectET/batch17_rna/4.proper

mkdir $featureCounts_dir
echo $hisat2_bam_dir
cd $hisat2_bam_dir
thread=2
for i in *.bam
do
echo $i
nohup featureCounts -T $thread -M -p -B -C -t exon -g gene_id -R BAM -a $ann_dir/mm10.gtf -o $featureCounts_dir/${i%%.bam}.featureCounts_geneCounts $i >$featureCounts_dir/${i%%.bam}.featureCounts.log &
done


6.bigwig
conda activate py36
wd=~/work_space/4.ProjectET/batch34_rna/4.proper
nwd=~/work_space/4.ProjectET/batch34_rna/z.bw
mkdir -p $nwd
cd $wd
for i in *.bam
do
name=${i%.proper*}
echo $name
nohup bamCoverage  -p 8 -b $i -o $nwd/${name}.bw \
--normalizeUsing RPKM  --centerReads --binSize 25 -ignore chrM > $nwd/${name}.bw.log &
done

################### To quantify repeats

#### merge fastq files
cd ~/work_space/4.ProjectET/batch34_rna/2.cutdata
mkdir ../2.1cut_merge_fastq
for i in *_rep1.R1_cut.fq
do
j=${i%_rep1*}
echo $j
ls *${j}*R1*
cat `ls *${j}*R1*` > ../2.1cut_merge_fastq/${j}_merge.R1.cut.fq &
cat `ls *${j}*R2*` > ../2.1cut_merge_fastq/${j}_merge.R2.cut.fq &
done

conda activate squire
wd=~/work_space/4.ProjectET/batch34_rna/2.1cut_merge_fastq
cd $wd
for i in *R1*fq
do
j=${i/R1/R2}
nohup squire Map -1 $i -2 $j -o $wd/../3.1squire_map -f ~/software/squire_index/squire_fetch -r 140 -n ${i%_*} -b mm10 -p 30 &
done

cd $wd/../3.1squire_map
for i in *.bam
do
nohup squire Count -m ~/work_space/4.ProjectET/batch34_rna/3.1squire_map -c /home1/gyang/software/squire_index/squire_clean -o ~/work_space/4.ProjectET/batch34_rna/3.2squire_count -f ~/software/squire_index/squire_fetch -r 140 -n ${i%.*} -b mm10 -p 30 -t  /home1/gyang/software/squire_index/squire_clean  &
done


############################ FeatureCounts

featureCounts_dir=~/work_space/1.Mouse_Acetylation/batch22/4.1featureCounts
ann_dir=/home1/share/gtf
hisat2_bam_dir=~/work_space/1.Mouse_Acetylation/batch22/4.proper

mkdir $featureCounts_dir
echo $hisat2_bam_dir
cd $hisat2_bam_dir
thread=2
for i in *.bam
do
echo $i
nohup featureCounts -T $thread -M -p -B -C -t exon -g gene_id -R BAM -a $ann_dir/mm10.gtf -o $featureCounts_dir/${i%%.bam}.featureCounts_geneCounts $i >$featureCounts_dir/${i%%.bam}.featureCounts.log &
done

rm merge_geneCounts_featureCounts_hisat2.tab
file=merge_geneCounts_featureCounts_hisat2.tab

for i in *.featureCounts_geneCounts
do
k=${i%%.featureCounts_geneCounts}
echo $k
if [ ! -f "$file" ]; then
  cut -f 1,6,7 $i | sed '1,2d' | sort -k1,1 | sed '1i #id\tlength\t'$k > $file
else
  cut -f 1,7 $i | sed '1,2d' | sort -k1,1 | sed '1i #id\t'$k > temp
  join -t $'\t' $file temp > temp2
#  cut -f 2 temp | paste -d $'\t' $file - > temp2
  mv temp2 $file
fi
wc -l $file
done
rm temp


#### merge fastq files
cd ~/work_space/4.ProjectET/batch15_rna/2.cutdata
mkdir ../2.1cut_merge_fastq
for i in *_rep1.R1_cut.fq
do
j=${i%rep1*}
echo $j
ls *${j}*R1*
cat `ls *${j}*R1*` > ../2.1cut_merge_fastq/${j}_merge.R1.cut.fq &
cat `ls *${j}*R2*` > ../2.1cut_merge_fastq/${j}_merge.R2.cut.fq &
done

conda activate squire
wd=~/work_space/4.ProjectET/batch15_rna/2.1cut_merge_fastq
cd $wd
for i in *R1*fq
do
j=${i/R1/R2}
nohup squire Map -1 $i -2 $j -o $wd/../3.1squire_map -f ~/software/squire_index/squire_fetch -r 140 -n ${i%_*} -b mm10 -p 10 &
done

mkdir ~/work_space/4.ProjectET/batch15_rna/3.2squire_count
cd $wd/../3.1squire_map
for i in *.bam
do
echo $i
nohup squire Count -m  ~/work_space/4.ProjectET/batch15_rna/3.1squire_map -c /home1/gyang/software/squire_index/squire_clean -o ~/work_space/4.ProjectET/batch15_rna/3.2squire_count -f ~/software/squire_index/squire_fetch -r 140 -n ${i%.*} -b mm10 -p 10 -t  /home1/gyang/software/squire_index/squire_clean  &
done
