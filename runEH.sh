#!/bin/bash

GENE=${1}
JA="${2}"
JA="$(readlink -f "$JA")"
OUTDIR="${3}"
JOBINDEX=${4:-$LSB_JOBINDEX}
STEP=100

STOP=$(( JOBINDEX + STEP ))
IDX=$((JOBINDEX))
while [[ $IDX -lt $STOP ]];
do
    echo "$IDX"
    PULL_LINE=$((IDX))
    QUIT_LINE=$((PULL_LINE + 1))
    PROCESS_LINE=($(sed -n "$PULL_LINE"p;"$QUIT_LINE"q "$JA"))
    PLATEKEY="${PROCESS_LINE[0]}"
    GENOME_BUILD="${PROCESS_LINE[1]}"
    BAMFILE="${PROCESS_LINE[2]}"
    SEX="${PROCESS_LINE[4]}"
    echo $PLATEKEY
    echo $GENOME_BUILD
    if [ ${GENOME_BUILD} = "GRCh38" ]
    then
        reference="/genomes-bertha-dev/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa"
        json="${OUTDIR}/${GENE}.json"
    else
        reference="/genomes-bertha-dev/resources/genomeref/Illumina/Homo_sapiens/Ensembl/GRCh37/Sequence/WholeGenomeFasta/genome.fa"
        json="${OUTDIR}/${GENE}_grch37.json"
    fi
    mkdir ${OUTDIR}/${GENE}
    if [ ! -f ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf ]
    then
        ExpansionHunter --reads ${BAMFILE} \
            --reference $reference \
            --output-prefix ${OUTDIR}/${GENE}/${PLATEKEY}_${GENE} \
            --variant-catalog ${json} --threads 3 --sex ${SEX}
    else
        echo "${OUTDIR}/${GENE}/${PLATEKEY}_${GENE}.vcf exists"
    fi
    IDX=$((IDX + 1))
done

# for gene in BCLAF3_revised
#
# bsub -q short -P re_gecip_neurology -J "Array[1-20262:100]" -R rusage[mem=48000] -M 48000 -n 3 \
# -o "${PWD}/jobs/myJob.out%I" -e "${PWD}/jobs/myJob.err%I" \
# "${PWD}/eh.sh $gene ${PWD}/100k_paths_annotated ${PWD}/EH"
# done
