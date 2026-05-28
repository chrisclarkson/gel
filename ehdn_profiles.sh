#!/bin/bash

module load ExpansionHunterDenovo/0.9.0
# export REF_PATH="/tools/aws-workspace-apps/re_admin/source_code/samtools_ref_cache/ref/cache/%2s/%2s/%s:http://www.ebi.ac.uk/ena/cram/md5/%s"
# export REF_CACHE="/tools/aws-workspace-apps/re_admin/source_code/samtools_ref_cache/ref/cache/%2s/%2s/%s"


JA="$1"

JA="$(readlink -f "$JA")"
JOBINDEX=${1:-$LSB_JOBINDEX}
OUT_DIR="$2"
# JOBINDEX=50
# STEP=$(echo $LSB_BATCH_JID | cut -d'[' -f 2 | cut -d ']' -f1)
STEP=10
STOP=$(( JOBINDEX + STEP ))
IDX=$((JOBINDEX))
# IDX=1

while [[ $IDX -lt $STOP ]];
do
    echo $IDX
    PULL_LINE=$((IDX))
    QUIT_LINE=$((PULL_LINE + 1))
    PROCESS_LINE=($(sed -n "$PULL_LINE"p;"$QUIT_LINE"q "$JA"))
    echo ${PROCESS_LINE[2]}
    # if [ ! -f ${PROCESS_LINE[1]}/str_profiles/${PROCESS_LINE[0]}.str_profile.json ]
    # if [ -f GMS/non_neuros/str_profiles/${PROCESS_LINE[0]}.str_profile.json ]
    # then
    # samtools view -h -b -T /genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa ${PROCESS_LINE[2]} > crams_to_bams/${PROCESS_LINE[0]}.bam
    ExpansionHunterDenovo profile --reads ${PROCESS_LINE[2]} \
        --reference /genomes/resources/genomeref/Illumina/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa \
        --output-prefix ${OUT_DIR}/${PROCESS_LINE[0]} --min-anchor-mapq 50 --max-irr-mapq 40
    # else
    #     echo "${PROCESS_LINE[1]}/str_profiles/${PROCESS_LINE[0]}.str_profile.json exists"
    # fi
    IDX=$((IDX + 1))
done

# bsub -q short -P re_gecip_neurology -J "Array[1-N:10]" -R rusage[mem=48000] -M 48000 -n 3 \
# -o "${PWD}/jobs/myJob.out%I" -e "${PWD}/jobs/myJob.err%I" \
# "${PWD}/ehdn_profiles.sh ${PWD}/paths ${PWD}/EH"
