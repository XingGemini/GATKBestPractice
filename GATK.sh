#!/bin/bash
BASEDIR=$(dirname "$0")
CONFIG=$BASEDIR/config

source $CONFIG

CHROM=( chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY )

STDOUT='/dev/stdout'
STDIN='/dev/stdin'
TMPDIR=$PWD/tmp

HG19REF=$REFERENCEDIR/ucsc.hg19.fasta
HG19REFBWT=$REFERENCEDIR/ucsc.hg19.fasta.bwt
HG19REFFAI=$REFERENCEDIR/ucsc.hg19.fasta.fai
HG19REFDICT=$REFERENCEDIR/ucsc.hg19.dict
ENRICHMENTDIR=$REFERENCEDIR/hg19/Enrichment/SureSelectV5'
ENRICHMENTBED=$ENRICHMENTDIR/S04380110_Regions.bed

KNOWNINDEL1000GVCF=$REFERENCEDIR/1000G_phase1.indels.hg19.sites.vcf
KNOWNINDELMILLSVCF=$REFERENCEDIR/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf
KNOWNDBSNPVCF=$REFERENCEDIR/dbsnp_138.hg19.vcf

FASTQ1=$1
FASTQ2=$2
SAMPLEID=$3
STEP=$4
FASTQBAM=$TMPDIR/$SAMPLEID.raw.bam
MARKADAPTERBAM=$TMPDIR/$SAMPLEID.markIlluminaAdapter.bam
MARKADAPTERRPT=$TMPDIR/$SAMPLEID\_markIlluminaAdapter.metrics.txt
MAPPEDBAM=$TMPDIR/$SAMPLEID.mapped.bam
SORTEDBAM=$TMPDIR/$SAMPLEID.mapped.sorted.bam
MARKDUPBAM=$TMPDIR/$SAMPLEID.mapped.sorted.markdup.bam
MARKDUPRPT=$TMPDIR/$SAMPLEID.mapped.sorted.markdup.metrics.txt
#REALNEDBAM=$TMPDIR/$SAMPLEID.mapped.sorted.markdup.realn.bam
FINALBAM=$PWD/$SAMPLEID.mapped.sorted.markdup.realn.recal.bam
GVCFGZ=$PWD/$SAMPLEID.g.vcf.gz
RAWVARIANTVCFGZ=$TMPDIR/$SAMPLEID.raw.variant.vcf.gz
RAWINDELVCFGZ=$TMPDIR/$SAMPLEID.raw.indel.vcf.gz
RAWSNPVCFGZ=$TMPDIR/$SAMPLEID.raw.snp.vcf.gz
FILTEREDVARIANTVCFGZ=$PWD/$SAMPLEID.filtered.variant.vcf.gz
FILTEREDINDELVCFGZ=$TMPDIR/$SAMPLEID.filtered.indel.vcf.gz
FILTEREDSNPVCFGZ=$TMPDIR/$SAMPLEID.filtered.snp.vcf.gz

USAGE="$0 <FASTQ1> <FASTQ2> <SAMPLEID> <STEP>  where STEP can be step number or \'all\'"

if [ -z "$1" ]; then
     echo $USAGE
    exit
fi

if [ ! -d "$TMPDIR" ]; then
    mkdir $TMPDIR
fi

if [ $STEP == 0 -o $STEP == 'all' ]; then
    echo 'Step 0: Fastq to Sam'
    $PICARD FastqToSam \
        FASTQ=$FASTQ1 \
        FASTQ2=$FASTQ2 \
        OUTPUT=$FASTQBAM \
        READ_GROUP_NAME=$SAMPLEID \
        SAMPLE_NAME=$SAMPLEID \
        LIBRARY_NAME=$SAMPLEID \
        PLATFORM=ILLUMINA \
        PLATFORM_UNIT=FlowCell1_Lane1_1 \
        SEQUENCING_CENTER=MACROGEN_MD \
        RUN_DATE=$(date +%Y%m%d) \
        &> Step0.log
fi

if [ $STEP == 1 -o $STEP == 'all' ]; then
    echo 'Step 1: Mark ILLUMINA adpaters'
    $PICARD MarkIlluminaAdapters \
        I=$FASTQBAM \
        O=$MARKADAPTERBAM \
        M=$MARKADAPTERRPT \
        TMP_DIR=$TMPDIR \
        &> Step1.log
fi

if [ $STEP == 2 -o $STEP == 'all' ]; then
    echo 'Step 2: Mapping to Reference'
    if [ ! -f $HG19REFBWT ]; then
        echo 'Step 2.0.1  prepare bwa reference directory'
        $BWA index -a bwtsw $HG19REF
    fi

    if [ ! -f $HG19REFFAI ]; then
        echo 'Step 2.0.2 prepare faidx for reference fastq'
        $SAMTOOLS faidx $HG19REF
    fi

    if [ ! -f $HG19REFDICT ]; then
        echo 'Step 2.0.3 prepare picard dict for reference fastq'
        $PICARD CreateSequenceDictionary \
            R=$HG19REF \
            O=$HG19REFDICT
    fi

    echo 'Step 2.1 Mapping to reference'
    set -o pipefail
    $PICARD SamToFastq \
        I=$MARKADAPTERBAM \
        FASTQ=$STDOUT \
        CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 INTERLEAVE=true NON_PF=true \
        TMP_DIR=$TMPDIR | \
    #$BWA mem -M -t 32 -p $HG19REF $STDIN | \
    $BWA mem -M -t 24 -p $HG19REF $STDIN | \

    $PICARD MergeBamAlignment \
        ALIGNED_BAM=$STDIN \
        UNMAPPED_BAM=$FASTQBAM \
        OUTPUT=$MAPPEDBAM \
        R=$HG19REF  CREATE_INDEX=true ADD_MATE_CIGAR=true \
        CLIP_ADAPTERS=false CLIP_OVERLAPPING_READS=true \
        INCLUDE_SECONDARY_ALIGNMENTS=true MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS \
        TMP_DIR=$TMPDIR \
        &> Step2.log
fi

if [ $STEP == 3 -o $STEP == 'all' ]; then
    echo 'Step 3: MarkDuplicates'

    echo 'Step 3.1: Sort aligned bam'
    $PICARD SortSam \
        I=$MAPPEDBAM \
        O=$SORTEDBAM \
        SORT_ORDER=coordinate \
        &> Step3.1.log

    echo 'Step 3.2: Mark Duplicates'
    $PICARD MarkDuplicates \
        I=$SORTEDBAM \
        O=$MARKDUPBAM \
        METRICS_FILE=$MARKADAPTERRPT \
        &> Step3.2.log

    echo 'Step 3.3: Build Bam Index'
    $PICARD BuildBamIndex \
        I=$MARKDUPBAM \
        &> Step3.3.log
fi

if [ $STEP == 4 -o $STEP == 'all' ]; then
    echo 'Step 4: Realgin indel'

    for i in "${CHROM[@]}"
    do
        if [ ! -f $ENRICHMENTDIR/$i\.bed ];
        then
            echo $i

            echo 'grep "^$i\t" $ENRICHMENTBED |cut -f 1-3 > $ENRICHMENTDIR\/$i\.bed  '
            grep -P "^$i\t" $ENRICHMENTBED |cut -f 1-3 > $ENRICHMENTDIR\/$i\.bed  # grep split bed file in chromosome files and keep only the first 3 columns
        fi
    done

    echo 'Step 4: Logging'  > Step4.log
    for i in "${CHROM[@]}"
    do
        echo "Building Target List for $i"
        $GATK -T RealignerTargetCreator \
            -R $HG19REF \
            -I $MARKDUPBAM \
            -L $ENRICHMENTDIR/$i\.bed \
            -known $KNOWNINDELMILLSVCF \
            -known $KNOWNINDEL1000GVCF \
            -o $TMPDIR/$i\.realignment_targets.list \
            &> tmp.log
        cat tmp.log >> Step4.log

        echo "Realgin indels for $i"
        $GATK -T IndelRealigner \
            -R $HG19REF \
            -I $MARKDUPBAM \
            -L $ENRICHMENTDIR/$i\.bed \
            -known $KNOWNINDELMILLSVCF \
            -known $KNOWNINDEL1000GVCF \
            -targetIntervals $TMPDIR/$i\.realignment_targets.list \
            -o $TMPDIR/$i\.mapped.sorted.markdup.realn.bam \
            &> tmp.log
        cat tmp.log >> Step4.log
    done
fi

if [ $STEP == 5 -o $STEP == 'all' ]; then
    echo "Step5 Recalibration Base Quality"

    echo 'Step 5: Logging'  > Step5.log

    for i in "${CHROM[@]}"
    do
        echo "5.1 Build recalibration base quality table for $i"
        $GATK -T BaseRecalibrator \
            -R $HG19REF \
            -I $TMPDIR/$i\.mapped.sorted.markdup.realn.bam \
            -knownSites $KNOWNINDELMILLSVCF \
            -knownSites $KNOWNINDEL1000GVCF \
            -knownSites $KNOWNDBSNPVCF \
            -o $TMPDIR/$i\.mapped.sorted.markdup.realn.recal.table \
            &> tmp.log
        cat tmp.log >> Step5.log
    done

    for i in "${CHROM[@]}"
    do
        echo "5.2 Recalibration base quality for $i"
        $GATK -T PrintReads \
            -R $HG19REF \
            -I $TMPDIR/$i\.mapped.sorted.markdup.realn.bam \
            -BQSR $TMPDIR/$i\.mapped.sorted.markdup.realn.recal.table \
            -o $TMPDIR/$i\.mapped.sorted.markdup.realn.recal.bam
        cat tmp.log >> Step5.log
    done

    BAMLIST=''
    for i in "${CHROM[@]}"
    do
        BAMLIST+="  $TMPDIR/$i.mapped.sorted.markdup.realn.recal.bam"
    done

    echo '5.3 Concadinate all individual Bam files'
    echo "    $SAMTOOLS cat -o $FINALBAM $BAMLIST >> Step5.log"
    $SAMTOOLS cat -o $FINALBAM $BAMLIST >> Step5.log

    echo 'Step 5.4: Build Bam Index'
    $PICARD BuildBamIndex \
        I=$FINALBAM \
        &> tmp.log
    cat tmp.log >> Step5.log

    #rm -f $BAMLIST
fi

if [ $STEP == 6 -o $STEP == 'all' ]; then
    echo "Step6: Variant Calling"

    echo "Step6: Variant Calling Logging" > Step6.log

    VCFLIST=''
    for i in "${CHROM[@]}"
    do
        echo "6.1 Haplotype Calling for $i"
        CHRGVCFGZ=$TMPDIR/$SAMPLEID.$i\.g.vcf.gz
        VCFLIST+="-V $CHRGVCFGZ "
        $GATK -T HaplotypeCaller \
            -R $HG19REF \
            -I $FINALBAM \
            -L $ENRICHMENTDIR/$i\.bed \
            --genotyping_mode DISCOVERY \
            --emitRefConfidence GVCF \
            --variant_index_type LINEAR --variant_index_parameter 128000 \
            -stand_emit_conf 10 \
            -stand_call_conf 30 \
            -o $CHRGVCFGZ \
            &> tmp.log
        cat tmp.log >> Step6.log
    done

    echo 'combine gvcf files'
    echo '/home/xing/jre1.8.0_91/bin/java -cp /home/xing/data/Tools/GATK/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
        -R $HG19REF \
        -assumeSorted \
        -o $RAWVARIANTVCFGZ \
        $VCFLIST \
        &> tmp.log'
    $GATKCP org.broadinstitute.gatk.tools.CatVariants \
        -R $HG19REF \
        -assumeSorted \
        -out $GVCFGZ \
        $VCFLIST \
        &> tmp.log
    cat tmp.log >> Step6.log

fi

if [ $STEP == 7 -o $STEP == 'all' ]; then
    echo "Step7: Genotype Calling"

    echo "Step7: Genotype Calling Log" > Step7.log

    VCFLIST=''
    for i in "${CHROM[@]}"
    do
        echo "7.1 Genotype Calling for $i"
        CHRRAWVARIANTSVCFGZ=$TMPDIR/$SAMPLEID.$i\.raw.variants.vcf.gz
        VCFLIST+="-V $CHRRAWVARIANTSVCFGZ "
        $GATK -T GenotypeGVCFs \
            -R $HG19REF \
            -V $TMPDIR/$SAMPLEID.$i\.g.vcf.gz \
            -o $CHRRAWVARIANTSVCFGZ \
            &> tmp.log
        cat tmp.log >> Step7.log
    done

    echo 'combine raw vcf files'
    $GATKCP org.broadinstitute.gatk.tools.CatVariants \
        -R $HG19REF \
        -assumeSorted \
        -out $RAWVARIANTVCFGZ \
        $VCFLIST \
        &> tmp.log
    cat tmp.log >> Step7.log
fi

if [ $STEP == 8 -o $STEP == 'all' ]; then
    echo "Step8: Filter Indels"

    echo "Step8: Filter Indels " > Step8.log

    $GATK -T SelectVariants \
        -R $HG19REF \
        -V $RAWVARIANTVCFGZ \
        -selectType INDEL \
        -o $RAWINDELVCFGZ \
        &> tmp.log
    cat tmp.log >> Step8.log

    $GATK -T VariantFiltration \
        -R $HG19REF \
        -V $RAWINDELVCFGZ \
        -filter "QD<2.0 || FS>200.0 ||SOR > 10.0 || InbreedingCoeff< -0.8 || ReadPosRankSum< -20.0" \
        -filterName "MG_indel_filter" \
        -o $FILTEREDINDELVCFGZ \
        &> tmp.log
    cat tmp.log >> Step8.log
fi

if [ $STEP == 9 -o $STEP == 'all' ]; then
    echo "Step9: Filter SNPs"

    echo "Step9: Filter SNPs " > Step9.log

    $GATK -T SelectVariants \
        -R $HG19REF \
        -V $RAWVARIANTVCFGZ \
        -selectType SNP \
        -o $RAWSNPVCFGZ \
        &> tmp.log
    cat tmp.log >> Step9.log

    $GATK -T VariantFiltration \
        -R $HG19REF \
        -V $RAWSNPVCFGZ \
        -filter "QD < 2.0 || FS > 60.0 || SOR > 4.0 ||  MQ < 40.0 || MQRankSum< -12.5 || ReadPosRankSum< -8.0" \
        -filterName "MG_SNP_filter" \
        -o $FILTEREDSNPVCFGZ \
        &> tmp.log
    cat tmp.log >> Step9.log
fi

if [ $STEP == 10 -o $STEP == 'all' ]; then
    echo "Step10: Merge Variants"

    echo "Step10: Merge Variants" > Step10.log

    $GATK -T CombineVariants \
        -R $HG19REF \
        -V $FILTEREDINDELVCFGZ \
        -V $FILTEREDSNPVCFGZ \
        -o $FILTEREDVARIANTVCFGZ \
        -genotypeMergeOptions UNSORTED \
        &> tmp.log
    cat tmp.log >> Step10.log
fi
