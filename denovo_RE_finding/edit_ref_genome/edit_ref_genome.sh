#!/bin/bash
module load samtools/1.16.1
module load bedtools/2.31.0
module load BWA/0.7.17
module load picard/3.1.1

# rm mod.fa
# python edit_ref_genome.py 40
# samtools faidx mod.fa
# bwa index mod.fa


JA="${1}"
JA="$(readlink -f "$JA")"
JOBINDEX=2
STEP=`wc -l ${JA} | awk -F' ' '{print $1}'`
# STEP=100

STOP=$(( JOBINDEX + STEP + 1 ))
IDX=$((JOBINDEX))
while [[ $IDX -lt $STOP ]];
do
    echo "$IDX"
    PULL_LINE=$((IDX))
    QUIT_LINE=$((PULL_LINE + 1))
    PROCESS_LINE=($(sed -n "${PULL_LINE}p;${QUIT_LINE}q" "$JA"))
    PLATEKEY="${PROCESS_LINE[0]}"
    GENOME_BUILD="${PROCESS_LINE[1]}"
    BAMFILE="${PROCESS_LINE[2]}"
    samtools view -h -bS $BAMFILE chr17:50198604-50200000 > test.bam
    java -jar /tools/aws-workspace-apps/re_admin/source_code/picard/3.1.1/picard.jar \
        SamToFastq I=test.bam FASTQ=Pair1_picard.fastq     SECOND_END_FASTQ=Pair2_picard.fastq \
        UNPAIRED_FASTQ=UNpaired_picard.fastq     INCLUDE_NON_PF_READS=true     VALIDATION_STRINGENCY=LENIENT
    # wc -l Pair1_picard.fastq
    # bwa index mod.fa
    # reference="/genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa"
    reference="mod.fa"
    bwa mem -Y $reference Pair1_picard.fastq Pair2_picard.fastq > test_reads.sam
    samtools sort -o test_reads_std.bam test_reads.sam
    samtools index test_reads_std.bam
    reference="mod.fa"
    json="CC_AT/EH/COL1A1_mod_rc.json"
    /nas/weka.gel.zone/re_gecip/shared_allGeCIPs/AT_CC/congenital_heart_disease/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter --reads test_reads_std.bam \
            --reference ${reference} \
            --output-prefix test \
            --variant-catalog ${json} --threads 3 \
            -a path-aligner
    tail -n -1 test.vcf >> col1a1_outs
    GENE="COL1A1_rc"
    json="CC_AT/EH/${GENE}_mod_rc.json"

    bam_path="test_realigned.bam"
    # bam_path="ensure_not_solveable/known_genes/${bam}_known_genes_realigned.bam"
    python3 /re_gecip/neurology/Valentina/GraphAlignmentViewer/GraphAlignmentViewer_record.py \
        --read_align $bam_path \
        --variant_catalog $json \
        --locus_id COL1A1 \
        --reference_fasta mod.fa \
        --file_format v3 \
        --output_prefix test_${PLATEKEY}
    IDX=$((IDX + 1))
done

