Segerated pre-rRNA processing in nucleolar subdomains coordinates with celluar demands
======================================================================================
Editing date: 20250706  
Author: Yu-Yao Zhang  
Email: zhangyuyao24@m.fudan.edu.cn

1.RNA-seq profiling of H9-derived arcuate neurons at differentiation day 0 (D0) and day 30 (D30) to analyze pre-rRNA processing factor expression dynamics
-------------------------------------
## 1.1 Requirements
------------------
```
Trimmomatics v0.32
Bowtie v1.1.2
HISAT2 v2.2.1
samtools v1.6
featureCounts v2.0.1
```
## 1.2 Trim adapter
```
parallel -j 3 'java -jar /picb/rnomics1/xuew/software/Trimmomatic-0.32/trimmomatic-0.32.jar PE -threads 16 -phred33  /path/to/{}_R1.fastq.gz /path/to/{}_R2.fastq.gz --baseout {}.qc.fq.gz ILLUMINACLIP:/path/to/Trimmomatic-0.32/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30' ::: $sample_list
```
### 1.3 Remove rRNA
map reads to human 45S rRNA (NCBI:NR_046235.3) to remove rRNA reads
```
parallel -j 3 'bowtie2 -x /path/to/bowtie_index  -1 /path/to/{}.qc_1P.fq.gz -2 /path/to/{}.qc_2P.fq.gz  -m 1 -k 1 -v 2 -S -p 12 -S {}.sam --un-conc-gz {}.fq.gz' ::: $sample_list 
```
### 1.4 Map to genome
map reads excluding rRNA reads to hg38 genome 
```
parallel -j 3 'hisat2 -x /path/to/hisat2_index --no-softclip --rna-strandness RF --score-min L,-16,0 --mp 7,7 --rfg 0,7 --rdg 0,7 --max-seeds 20 -k 10 -t -p 12 -1 /path/to/{}.fq.1.gz -2 /path/to/{}.fq.2.gz -S {}.bam' ::: $sample_list
parallel -j 3 'samtools view -bS  {} | samtools sort -@ 15 - -O BAM -o {= s:.sam::; =}.sorted.bam' ::: *.sam
parallel -j 3 'samtools index {}' ::: *.bam
parallel 'bamtools stats -in {} > {}.log' ::: *.bam
paste <(grep "Mapped reads" *log |cut -f1 -d":") <(grep "Mapped reads" *.log|cut -f3 -d":"|sed -E -e 's/\s+//g' -e 's/\(.*\)//'|xargs -I {} echo {}/2 |bc)|sed 's/\..*\t/\t/'    > sample_mappable_stats.txt
```
### 1.5 Compute read count matrix
quantify from bam files based on gencode v41 gene annotation
```
featureCounts -s 2 -p --fraction -O -T 16 -a /path/to/gencode.v41.annotation.gtf -t exon -g gene_id -o linear-RNA_Count.txt /path/to/*.bam
export mappable_frags=`seq 1 6|xargs|sed 's/ /,/g'`
tail -n +3 linear-RNA_Count.txt | perl -alne '@total=split(/,/,$ENV{'mappable_frags'}); @r=();for $i (0..$#total){$s=1;push @r, $F[$i+6]*$s};$,="\t";print ($F[0],@r)' > Counts.txt
sed -i -E "1iEnsemblID\t$sample_list" Counts.txt
```
### 1.6 Get differentially expressed pre-rRNA processing factors
performed in R v4.3
```
source('DESeq2.R')

count <-
read.table('/path/to/Counts.txt',header=T)

depth <-
read.table('/path/to/sample_mappable_stats.txt',col.names=c('sample','depth'))%>%
mutate(sample=str_replace(sample,pattern = "-",replacement = "."))

cpm <-
count%>%
gather(key = 'sample',value = 'count',-1)%>%
inner_join(depth)%>%
mutate(cpm=count*1000000/depth)%>%
select(1,2,5)%>%
spread(key = sample,value = cpm)

highexpressed_cpm <-
cpm[apply(cpm[, 2:7], 1, function(x) any(x > 1)), ]

highexpressed_count <-
count[count$EnsemblID %in% highexpressed_cpm$EnsemblID,]

highexpressed_count <-
highexpressed_count%>%
mutate(EnsemblID=str_remove(EnsemblID,'\\..*'))

highexpressed_count$Gene.name <-
mapIds(org.Hs.eg.db, keys = highexpressed_count$EnsemblID, 
       column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

highexpressed_count <-
highexpressed_count%>%
mutate(Gene.name = ifelse(is.na(Gene.name), EnsemblID, Gene.name))

sample <-
get_coldata(highexpressed_count[-1],c('time','rep'))

DEG <-
Count2DEG(highexpressed_count,col = 'Gene.name')%>%
DEG(colData = sample,formula = ~time)%>%
DESeq2fillNA

pre_rRNA_processing_gene<-
read.table('/path/to/pre_rRNA_processing_factors.txt')

pre-rRNA_processing_factors_DEG <-
left_join(pre_rRNA_processing_gene[1],DEG,by=c('Gene.name'='Row.names'))
```
