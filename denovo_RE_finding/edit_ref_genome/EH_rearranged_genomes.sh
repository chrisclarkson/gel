#!/bin/bash
module load samtools/1.16.1
module load bedtools/2.31.0
module load BWA/0.7.17
module load picard/3.1.1
# source /resources/conda/miniconda3/bin/activate

# conda activate /re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/test_env/
# rm mod.fa
# python edit_ref_genome.py 40
# samtools faidx mod.fa
# bwa index mod.fa
cd /re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined
GENE="${1}"
JA="${2}"
JA="$(readlink -f "$JA")"

region="chr17:50198604-50200000"
JOBINDEX=${4:-$LSB_JOBINDEX}
STEP=`wc -l ${JA} | awk -F' ' '{print $1}'`
# STEP=100

OUTDIR="${3}"
mkdir "${OUTDIR}/${GENE}"

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
    echo $PLATEKEY

    # samtools view -h -bS $BAMFILE $region > ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_extract.bam
    # java -jar /tools/aws-workspace-apps/re_admin/source_code/picard/3.1.1/picard.jar \
    #     SamToFastq I=${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_extract.bam FASTQ=${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair1_picard.fastq     SECOND_END_FASTQ=${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair2_picard.fastq \
    #     UNPAIRED_FASTQ=${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_UNpaired_picard.fastq     INCLUDE_NON_PF_READS=true     VALIDATION_STRINGENCY=LENIENT
    # reference="/re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/mod.fa"
    # bwa mem -Y $reference ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair1_picard.fastq ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair2_picard.fastq > ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads.sam
    # samtools sort -o ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads_std.bam ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads.sam
    # samtools index ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads_std.bam
    json="${OUTDIR}/${GENE}_mod_rc.json"
    # /nas/weka.gel.zone/re_gecip/shared_allGeCIPs/AT_CC/congenital_heart_disease/ExpansionHunter-v5.0.0-linux_x86_64/bin/ExpansionHunter \
    #         --reads ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads_std.bam \
    #         --reference ${reference} \
    #         --output-prefix ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE} \
    #         --variant-catalog ${json} --threads 3 \
    #         -a path-aligner

    bam_path="${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_realigned.bam"
    # bam_path="ensure_not_solveable/known_genes/${bam}_known_genes_realigned.bam"
    python3 /re_gecip/neurology/Valentina/GraphAlignmentViewer/GraphAlignmentViewer_record.py \
        --read_align $bam_path \
        --variant_catalog $json \
        --locus_id ${GENE}_rc \
        --reference_fasta $reference \
        --file_format v3 \
        --output_prefix ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}

    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair1_picard.fastq
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_Pair2_picard.fastq
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_UNpaired_picard.fastq
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_extract.bam
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads.sam
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads_std.bam
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_extract.bam.bai
    rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_reads_std.bam.bai
    IDX=$((IDX + 1))
done

# gene="COL1A1"
# bsub -q short -P re_gecip_neurology -J "Array[1-23831:100]" -R rusage[mem=48000] -M 48000 -n 3 \
#     -o "${PWD}/jobs_rearranged/myJob.out%I" -e "${PWD}/jobs_rearranged/myJob.err%I" \
#     "${PWD}/EH_rearranged_genomes.sh $gene ${PWD}/gms_jobs_updated.tsv ${PWD}/EH_rearranged_genomes"
# bsub -q short -P re_gecip_neurology -J "Array[1-34794:100]" -R rusage[mem=48000] -M 48000 -n 3 \
#     -o "${PWD}/jobs_rearranged/myJob.out%I" -e "${PWD}/jobs_rearranged/myJob.err%I" \
#     "${PWD}/EH_rearranged_genomes.sh $gene ${PWD}/neuro_platekeys.txt ${PWD}/EH_rearranged_genomes"

