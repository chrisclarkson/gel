#!/bin/bash

# JA="all_patients_with_codes_collapsed.tsv"
# JA="$(readlink -f "$JA")"
# JOBINDEX=${1:-$LSB_JOBINDEX}
JA="${2}"
JA="$(readlink -f "$JA")"
JOBINDEX=${4:-$LSB_JOBINDEX}
# STEP=`wc -l ${JA} | awk -F' ' '{print $1}'`
STEP=100
GENE=${1}
OUTDIR="${3}"

STOP=$(( JOBINDEX + STEP ))
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

    if [ ${GENOME_BUILD} = "GRCh38" ]
    then
       reference="/genomes-bertha-dev/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa"
       json="${OUTDIR}/${GENE}.json"
    else
       reference="/genomes-bertha-dev/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
       json="${OUTDIR}/${GENE}_grch37.json"
    fi
    # rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf
    # nlines=`wc -l ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf | awk -F' ' '{print $1}'`
    # if [ $nlines -eq 9 ]
    # then
        # rm ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf
    # fi
    mkdir ${OUTDIR}/${GENE}
    if [ ! -f ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf ]
    then
    # /tools/aws-workspace-apps/re_admin/source_code/ExpansionHunter/3.2.2/ExpansionHunter-v3.2.2-linux_x86_64/bin/ExpansionHunter
    /re_gecip/neurology/Tucci_lab_tools/ExpansionHunter --reads ${BAMFILE} \
        --reference ${reference} \
        --output-prefix  ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE} \
        --variant-catalog ${json} --threads 3
    else
        echo "${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf exists"
    fi
    IDX=$((IDX + 1))
    # bam_path="${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}_realigned.bam"
    # # bam_path="ensure_not_solveable/known_genes/${bam}_known_genes_realigned.bam"
    # python3 /re_gecip/neurology/Valentina/GraphAlignmentViewer/GraphAlignmentViewer_record.py \
    #     --read_align $bam_path \
    #     --variant_catalog $json \
    #     --locus_id $GENE \
    #     --reference_fasta $reference \
    #     --file_format v3 \
    #     --output_prefix ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}
done
# done < /home/cclarkson/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/genes_to_explore_ataxia
# done
# for gene in DAB1 DACH1 APELA C1orf168; 
# do 
# /re_gecip/neurology/Tucci_lab_tools/liftover_json.sh EH/${gene}.json EH/${gene}_grch37.json; 
# done
# for gene in GIPC1
# do
# bsub -q short -P re_gecip_neurology -J "Array[1-49445:100]" -R rusage[mem=1000] -M 1000 -n 3 \
#     -o "${PWD}/jobs/myJob.out%I" -e "${PWD}/jobs/myJob.err%I" \
#     "${PWD}/eh.sh $gene ${PWD}/combined_jobs.tsv ${PWD}/EH"
# bsub -q short -P re_gecip_neurology -J "Array[1-15249:100]" -R rusage[mem=1000] -M 1000 -n 3 \
#     -o "${PWD}/jobs/myJob.out%I" -e "${PWD}/jobs/myJob.err%I" \
#     "${PWD}/eh.sh $gene ${PWD}/cancer_paths.txt ${PWD}/EH"
# done
# ./eh.sh THAP11 ./job_array.txt ./EH_output
