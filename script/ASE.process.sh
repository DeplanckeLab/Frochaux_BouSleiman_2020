# Rootdir as parameter of this script
rootdir=$1;

# Here path for the droso genome
genomedir=/data/genome/drosophila_melanogaster/dm3;
inputGTF=$genomedir/genesBDGP5.25.2014-05-23.gtf;
inputFASTA=$genomedir/dm3.Wolb.fa;
barcodefile=$rootdir/fastq/barcodes.txt;
picardir=/software/picard-2.17.8;

# Creating output directories
mkdir -p $rootdir/BAM/NXT0222/
mkdir -p $rootdir/BAM/NXT0254/
mkdir -p $rootdir/Results/QC/NXT0222/
mkdir -p $rootdir/Results/QC/NXT0254/
mkdir -p $rootdir/TMP/

# 0- Prepare data
gatk IndexFeatureFile -F $genomedir/dgrp2.vcf

# 1- Separate Alignment of 2 BRB-seq libraries
STAR --runMode alignReads --twopassMode Basic --outSAMmapqUnique 60 --runThreadN 8 --genomeDir $genomedir --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $rootdir/BAM/NXT0222/ --readFilesIn $rootdir/fastq/NXT0222/MF_VB_eQTL_S1_R2_001.fastq.gz
STAR --runMode alignReads --twopassMode Basic --outSAMmapqUnique 60 --runThreadN 8 --genomeDir $genomedir --outFilterMultimapNmax 1 --readFilesCommand zcat --outSAMtype BAM Unsorted --outFileNamePrefix $rootdir/BAM/NXT0254/ --readFilesIn $rootdir/fastq/NXT0254/MF_VB_eQTL_S1_R2_001.fastq.gz

# 2- Separate Annotation of BAM
# NXT0222
java -jar /software/BRBseqTools.1.3.jar AnnotateBAM -f $rootdir/fastq/NXT0222/MF_VB_eQTL_S1_R1_001.fastq.gz -b $rootdir/BAM/NXT0222/Aligned.out.bam -o $rootdir/BAM/NXT0222/ -c $barcodefile -gtf $inputGTF -p B -rg -lib NXT0222
# NXT0254
java -jar /software/BRBseqTools.1.3.jar AnnotateBAM -f $rootdir/fastq/NXT0254/MF_VB_eQTL_S1_R1_001.fastq.gz -b $rootdir/BAM/NXT0254/Aligned.out.bam -o $rootdir/BAM/NXT0254/ -c $barcodefile -gtf $inputGTF -p BU -UMI 15 -rg -lib NXT0254

# 3- Merging BAMs & sorting
samtools merge $rootdir/BAM/round_robin.merged.annotated.bam $rootdir/BAM/NXT0222/round_robin.annotated.sorted.bam $rootdir/BAM/NXT0254/round_robin.annotated.sorted.bam
samtools sort $rootdir/BAM/round_robin.merged.annotated.bam -o $rootdir/BAM/round_robin.merged.sorted.bam

# 4- Mark Duplicates
java -Djava.io.tmpdir=$rootdir/TMP/ -jar $picardir/picard.jar MarkDuplicates I=$rootdir/BAM/round_robin.merged.sorted.bam O=$rootdir/BAM/round_robin.merged.marked.dups.bam BARCODE_TAG=BC READ_ONE_BARCODE_TAG=BX CREATE_INDEX=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT M=$rootdir/BAM/round_robin.merged.marked.dups.metrics.txt

# 5- Split N CIGAR
gatk SplitNCigarReads --TMP_DIR $rootdir/TMP/ -R $inputFASTA -I $rootdir/BAM/round_robin.merged.marked.dups.bam -O $rootdir/BAM/round_robin.merged.split.N.cigars.bam

# 6- Base Recalibration
gatk BaseRecalibrator --TMP_DIR $rootdir/TMP/ --known-sites $genomedir/dgrp2.vcf -R $inputFASTA -I $rootdir/BAM/round_robin.merged.split.N.cigars.bam -O $rootdir/BAM/round_robin.merged.recalibration_report.table
gatk ApplyBQSR --TMP_DIR $rootdir/TMP/ -I $rootdir/BAM/round_robin.merged.split.N.cigars.bam -O $rootdir/BAM/round_robin.merged.recalibrated.bam -R $inputFASTA --bqsr-recal-file $rootdir/BAM/round_robin.merged.recalibration_report.table

# 7- Quality Control
fastqc --outdir=$rootdir/Results/QC/NXT0222/ $rootdir/fastq/NXT0222/MF_VB_eQTL_S1_R1_001.fastq.gz
fastqc --outdir=$rootdir/Results/QC/NXT0222/ $rootdir/fastq/NXT0222/MF_VB_eQTL_S1_R2_001.fastq.gz
fastqc --outdir=$rootdir/Results/QC/NXT0254/ $rootdir/fastq/NXT0254/MF_VB_eQTL_S1_R1_001.fastq.gz
fastqc --outdir=$rootdir/Results/QC/NXT0254/ $rootdir/fastq/NXT0254/MF_VB_eQTL_S1_R2_001.fastq.gz
fastqc --outdir=$rootdir/Results/QC/ $rootdir/BAM/round_robin.merged.recalibrated.bam

# 8- Demultiplexing (in case it is needed)
mkdir -p $rootdir/BAM/Demult
gatk SplitReads -I $rootdir/BAM/round_robin.merged.recalibrated.bam -O $rootdir/BAM/Demult --split-sample

# 9- ASE Calculation
mkdir -p $rootdir/Results/ASE/

cross=(1x2_UC_1 2x3_UC_1 3x4_UC_1 4x5_UC_1 5x6_UC_1 6x7_UC_1 7x8_UC_1 8x9_UC_1 9x10_UC_1 10x11_UC_1 11x12_UC_1 12x13_UC_1 1x2_Pe_1 2x3_Pe_1 3x4_Pe_1 4x5_Pe_1 5x6_Pe_1 6x7_Pe_1 7x8_Pe_1 8x9_Pe_1 9x10_Pe_1 10x11_Pe_1 11x12_Pe_1 12x13_Pe_1 13x14_UC_1 14x15_UC_1 17x18_UC_1 18x19_UC_1 19x20_UC_1 20x1_UC_1 13x14_Pe_1 14x15_Pe_1 17x18_Pe_1 18x19_Pe_1 19x20_Pe_1 20x1_Pe_1 1x2_UC_2 2x3_UC_2 3x4_UC_2 4x5_UC_2 5x6_UC_2 6x7_UC_2 7x8_UC_2 8x9_UC_2 9x10_UC_2 10x11_UC_2 11x12_UC_2 12x13_UC_2 1x2_Pe_2 2x3_Pe_2 3x4_Pe_2 4x5_Pe_2 5x6_Pe_2 6x7_Pe_2 7x8_Pe_2 8x9_Pe_2 9x10_Pe_2 10x11_Pe_2 11x12_Pe_2 12x13_Pe_2 13x14_UC_2 14x15_UC_2 17x18_UC_2 18x19_UC_2 19x20_UC_2 20x1_UC_2 13x14_Pe_2 14x15_Pe_2 17x18_Pe_2 18x19_Pe_2 19x20_Pe_2 20x1_Pe_2)
sample1=(line_217 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_217 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_890 line_129 line_229 line_639 line_142 line_721 line_890 line_129 line_229 line_639 line_142 line_217 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_217 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_890 line_129 line_229 line_639 line_142 line_721 line_890 line_129 line_229 line_639 line_142)
sample2=(line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_890 line_332 line_229 line_639 line_142 line_217 line_890 line_332 line_229 line_639 line_142 line_217 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_386 line_426 line_287 line_802 line_208 line_313 line_707 line_379 line_486 line_535 line_738 line_721 line_890 line_332 line_229 line_639 line_142 line_217 line_890 line_332 line_229 line_639 line_142 line_217)

for i in {0..71}
do
	java -jar /software/ASECounter.0.1.jar -vcf $genomedir/dm3/dgrp2.vcf -bam $rootdir/BAM/Demult/round_robin.merged.recalibrated.${cross[i]}.bam -gtf $inputGTF -o $rootdir/Results/ASE/ -s1 ${sample1[i]} -s2 ${sample2[i]} -cross ${cross[i]}
done

#java -jar /software/ASECounter.0.1.jar -vcf $genomedir/dm3/dgrp2.vcf -bam $rootdir/BAM/Demult/round_robin.merged.recalibrated.1x2_Pe_1.bam -gtf $inputGTF -o $rootdir/Results/ASE/ -s1 line_217 -s2 line_386 -cross "1x2_Pe_ 1"
#java -jar /software/ASECounter.0.1.jar -vcf $genomedir/dm3/dgrp2.vcf -bam $rootdir/BAM/Demult/round_robin.merged.recalibrated.1x2_UC_2.bam -gtf $inputGTF -o $rootdir/Results/ASE/ -s1 line_217 -s2 line_386 -cross "1x2_UC_ 2"
