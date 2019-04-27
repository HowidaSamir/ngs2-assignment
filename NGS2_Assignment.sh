#!/bin/bash
#Explore the sample names
mkdir -p ~/home/ngs-01/workdir/ngs2-assignment-data && cd ~/home/ngs-01/workdir/ngs2-assignment-data/
for ((i=1;i<=5;i++)) ;
do
ls -tral ~/workdir/ngs2-assignment-data/shuffled_SRR8797509_*.part_001.part_001.fastq.gz
done

for ((i=1;i<=5;i++)) ;
do
ls -tral ~/workdir/ngs2-assignment-data/SRR8797509_*.part_001.part_001.fastq.gz
done

#Add and align all reads
#install STAR
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.7.0f.tar.gz
tar -xzf 2.7.0f.tar.gz
cd STAR-2.7.0f
# Compile under Linux
cd STAR/source
make STAR

#the shuffled
for R1 in ~/workdir/ngs2-assignment-data/shuffled_SRR8797509_*.part_001.part_001.fastq.gz;do
SM=$(basename "$R1" | cut -d"_" -f1)
LB=$(basename "$R1" | cut -d"_" -f1,2)
PL="Illumina"
RGID=$(zcat "$R1" | head -n1 | sed 's/:/_/g' |cut -d "_" -f1,2,3,4)
PU="$RGID"."$LB"
echo -e "@RG\tID:$RGID\tSM:$SM\tPL:$PL\tLB:$LB\tPU:$PU"
R2=$(echo "$R1" | sed 's/_R1_/_R2_/')
echo "$R1" "$R2"

mkdir -p ~/home/ngs-01/workdir/ngs2-assignment-data/STAR_align/STARIndex && cd ~/home/ngs-01/workdir/ngs2-assignment-data/STAR_align/STARIndex/
STAR --runThreadN {5} --runMode genomeGenerate --genomeDir /home/ngs-01/workdir/ngs2-assignment-data/STAR_align/STARIndex/ --genomeFastaFiles /home/ngs-01/workdir/ngs2-assignment-data/gencode.v29.pc_transcripts.chr22.simplified.fa --sjdbGTFfile /home/ngs-01/workdir/ngs2-assignment-data/ --sjdbOverhang {100 - 1}
STAR --runMode alignReads --outSAMtype BAM Unsorted --readFilesCommand zcat --genomeDir /home/ngs-01/workdir/ngs2-assignment-data/STAR_align/STARIndex --outFileNamePrefix {shuffled_SRR8797509_*.part_001.part_001}  --readFilesIn  /home/ngs-01/workdir/ngs2-assignment-data/shuffled_SRR8797509_1.part_001.part_001.fastq.gz /home/ngs-01/workdir/ngs2-assignment-data/shuffled_SRR8797509_2.part_001.part_001.fastq.gz
done

#P.S: the BAM files can be sorted directly through the STAR by this command: --outSAMtype BAM SortedByCoordinate

#Generate & sort BAM files:
for samfile in *.sam;do
  sample=${samfile%.sam}
  samtools view -hbo $sample.bam $samfile
  samtools sort $sample.bam -o $sample.sorted.bam
done

#Merge replicate with Picard tools:
# Install Picard tools
conda install -c bioconda picard 
picard_path=$CONDA_PREFIX/share/picard-2.19.0-0

# merge the replicates
java  -Xmx2g -jar $picard_path/picard.jar MergeSamFiles I=BD143_TGACCA_L005.sorted.bam I=BD143_TGACCA_L006.sorted.bam OUTPUT=BD143_TGACCA_merged.sorted.bam

# check for the changes in the header
samtools view -H BD143_TGACCA_L005.sorted.bam
samtools view -H BD143_TGACCA_L006.sorted.bam
samtools view -H BD143_TGACCA_merged.sorted.bam

# remove the individual replicates
rm BD143_TGACCA_L00*.sorted.bam

#mapping QC
for bamFile in *.sorted.bam;do
  output=${bamFile%.sorted.bam}
  samtools depth $bamFile | awk '{{sum+=$3}} END {{print "Average = ",sum/NR}}' > $output.cov
  samtools flagstat $bamFile > $output.stat
done

#Mark duplicate
for sample in *.sorted.bam;do
  name=${sample%.sorted.bam}
  java  -Xmx2g -jar $picard_path/picard.jar MarkDuplicates INPUT=$sample OUTPUT=$name.dedup.bam METRICS_FILE=$name.metrics.txt;
done

#Install GATK
conda install -c bioconda gatk4 

#indexing
# samples
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}
  java -Xmx2g -jar $picard_path/picard.jar BuildBamIndex VALIDATION_STRINGENCY=LENIENT INPUT=$sample
done

# Reference
ln -s ~/home/ngs-01/workdir/ngs2-assignment-data/gencode.v29.pc_transcripts.chr22.simplified.fa .
java -Xmx2g -jar $picard_path/picard.jar CreateSequenceDictionary R=gencode.v29.pc_transcripts.chr22.simplified.fa O=gencode.v29.pc_transcripts.chr22.simplified.dict
samtools faidx gencode.v29.pc_transcripts.chr22.simplified.fa

#Download known varinats
# Download known polymorphic sites
wget 'ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/hg18/1000G_phase1.snps.high_confidence.hg18.vcf.gz' -O 1000G_phase1.snps.high_confidence.hg18.vcf.gz

# Select variants on chr22 and correct chr name
gunzip 1000G_phase1.snps.high_confidence.hg18.vcf.gz
grep "^#" 1000G_phase1.snps.high_confidence.hg18.vcf > snps_hg18.vcf
grep "^22" 1000G_phase1.snps.high_confidence.hg18.vcf | sed 's/^22/chr22/' >> snps_hg18.vcf
gatk IndexFeatureFile -F snps_hg18.vcf

#Recalibrate Bases BQSR
for sample in *.dedup.bam;do
  name=${sample%.dedup.bam}

  gatk --java-options "-Xmx2G" BaseRecalibrator \
-R gencode.v29.pc_transcripts.chr22.simplified.fa -I $sample --known-sites snps_hg18.vcf \
-O $name.report

  gatk --java-options "-Xmx2G" ApplyBQSR \
-R dog_chr5.fa -I $sample -bqsr $name.report \
-O $name.bqsr.bam --add-output-sam-program-record --emit-original-quals
done

#Joint variant calling using HaplotypeCaller
## assess genotype likelihood per-sample
for sample in *.bqsr.bam;do
  name=${sample%.bqsr.bam}

  gatk --java-options "-Xmx2G" HaplotypeCaller \
  -R gencode.v29.pc_transcripts.chr22.simplified.fa -I $sample \
  --emit-ref-confidence GVCF \
  --pcr-indel-model NONE \
  -O $name.gvcf
done

## combine samples
gatk --java-options "-Xmx2G" CombineGVCFs \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V BD143_TGACCA_merged.gvcf \
-V BD174_CAGATC_L005.gvcf \
-V BD225_TAGCTT_L007.gvcf \
-O raw_variants.gvcf

## Joint Genotyping
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
-O raw_variants.vcf

## annotated output
gatk --java-options "-Xmx60G" GenotypeGVCFs \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V raw_variants.gvcf \
--max-alternate-alleles 6 \
--dbsnp snps_hg18.vcf \
-O raw_variants_ann.vcf

## check how many variant got annotated
grep -v "^#" raw_variants_ann.vcf | awk '{print $3}' | grep "^rs" | wc -l

#VCF statitics
#index the VCF file
conda install -c bioconda tabix
bgzip -c raw_variants_ann.vcf > raw_variants_ann.vcf.gz
tabix -p vcf raw_variants_ann.vcf.gz

#Calculate statistics about the vcf
conda install -c bioconda rtg-tools
rtg vcfstats raw_variants_ann.vcf.gz > stats.txt

#Split SNPs and indels
gatk --java-options "-Xmx2G" SelectVariants \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V raw_variants_ann.vcf \
--select-type-to-include SNP \
-O raw_variants_ann_SNP.vcf

gatk --java-options "-Xmx2G" SelectVariants \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V raw_variants_ann.vcf \
--select-type-to-include INDEL \
-O raw_variants_ann_INDEL.vcf

#Assess the different filters in both known and novel
for var in "SNP" "INDEL";do
 input="raw_variants_ann_"$var".vcf"
 for filter in "QD" "MQ" "MQRankSum" "FS" "SOR" "ReadPosRankSum" "AN" "DP" "InbreedingCoeff";do
  filterValues=$var.$filter
  awk -v k="$filter=" '!/#/{n=split($8,a,";"); for(i=1;i<=n;i++) if(a[i]~"^"k) {sub(k,$3" ",a[i]); print a[i]}}' $input > $filterValues
  grep -v "^\." $filterValues > known.$var.$filter
  grep "^\." $filterValues > novel.$var.$filter
done; done

mkdir filters && cd filters
mv ../{*.SNP.*,SNP.*,*.INDEL.*,INDEL.*} .

#to generate figures
wget https://raw.githubusercontent.com/dib-lab/dogSeq/master/scripts/densityCurves.R
sudo Rscript -e "install.packages('ggplot2', contriburl=contrib.url('http://cran.r-project.org/'))"
for f in SNP.* INDEL.*;do
  Rscript densityCurves.R "$f"
done

#Calc the DP threathols
cat SNP.DP INDEL.DP | awk '{sum+= $2; sumsq+= ($2)^2} END { print sum/NR, sqrt((sumsq-sum^2/NR)/NR), sum/NR + 5*sqrt((sumsq-sum^2/NR)/NR) }' 


#SNP Variant filteration

cd ~/home/ngs-01/workdir/ngs2-assignment-data/
gatk --java-options "-Xmx2G" VariantFiltration \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
-V raw_variants_ann_SNP.vcf \
--filter-name "snpQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "snpMQ" \
--filter-expression "vc.hasAttribute('MQ') && MQ < 40.0" \
--filter-name "snpMQRankSum" \
--filter-expression "vc.hasAttribute('MQRankSum') && MQRankSum < -12.5" \
--filter-name "snpFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 60.0" \
--filter-name "snpSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 4.0" \
--filter-name "snpReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -8.0" \
--filter-name "snpDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_SNP_clean.vcf

#to check the filtered record
grep -v "^#" raw_variants_ann_SNP_clean.vcf | awk '{if($7!="PASS")print $0}'


#INDEL Variant filteration
cd ~/home/ngs-01/workdir/ngs2-assignment-data/
gatk --java-options "-Xmx2G" VariantFiltration \
-R gencode.v29.pc_transcripts.chr22.simplified.fa \
V raw_variants_ann_SNP.vcf \
--filter-name "indelQD" \
--filter-expression "vc.hasAttribute('QD') && QD < 2.0" \
--filter-name "indelFS" \
--filter-expression "vc.hasAttribute('FS') && FS > 200.0" \
--filter-name "indelSOR" \
--filter-expression "vc.hasAttribute('SOR') && SOR > 10.0" \
--filter-name "indelReadPosRankSum" \
--filter-expression "vc.hasAttribute('ReadPosRankSum') && ReadPosRankSum < -20.0" \
--filter-name "indelInbreedingCoeff" \
--filter-expression "vc.hasAttribute('InbreedingCoeff') && InbreedingCoeff < -0.8" \
--filter-name "indelDP" \
--filter-expression "vc.hasAttribute('DP') && DP > 3105" \
-O raw_variants_ann_INDEL_clean.vcf
