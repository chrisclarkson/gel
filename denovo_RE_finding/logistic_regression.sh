#!/bin/bash
echo $PWD
source /resources/conda/miniconda3/bin/activate

conda activate /home/cclarkson/test_env

cd /re_gecip/shared_allGeCIPs/FA_AT/EHdenovo
JA="gms_and_100k_no_relations_europeans_casecontrol_locus_min5_annotated.tsv"
# JA="unfinished_jobs.txt"
JA="$(readlink -f "$JA")"
# JOBINDEX=${LSB_JOBINDEX}
JOBINDEX=${1:-$LSB_JOBINDEX}
# JOBINDEX=$1
PULL_LINE=$((JOBINDEX + 1))
QUIT_LINE=$((PULL_LINE + 1))
PROCESS_LINE=($(sed -n "${PULL_LINE}p;${QUIT_LINE}q" "$JA"))
echo $PWD
STEP=10
STOP=$(( JOBINDEX + STEP ))
IDX=$((JOBINDEX))
while [[ $IDX -lt $STOP ]];
do
    echo "$IDX"
    PULL_LINE=$((IDX))
    QUIT_LINE=$((PULL_LINE + 1))
    PROCESS_LINE=($(sed -n "${PULL_LINE}p;${QUIT_LINE}q" "$JA"))
    chrom=${PROCESS_LINE[0]}
    start=${PROCESS_LINE[1]}
    end=${PROCESS_LINE[2]}
    motif=${PROCESS_LINE[3]}
    gene=${PROCESS_LINE[4]}
    if [ ! -f ${PWD}/logistic_regression_output/${chrom}:${start}-${end}_casecontrol_100k.tsv ]
    then
        python /re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/analysis_scripts/motif_count_compare/motif_count_compare/main_logistic_regression.py \
                --motif $motif \
                --region ${chrom}:${start}-${end} \
                --output ${PWD}/logistic_regression_output/${chrom}:${start}-${end}_${gene}_${motif}_casecontrol_gms \
                --manifest manifest_gms_no_rels --export_ranked_genomes True
        python /re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/analysis_scripts/motif_count_compare/motif_count_compare/main_logistic_regression.py \
                --motif $motif \
                --region ${chrom}:${start}-${end} \
                --output ${PWD}/logistic_regression_output/${chrom}:${start}-${end}_${gene}_${motif}_casecontrol_100k \
                --manifest manifest_100k_no_rels --export_ranked_genomes True
    else
    echo "${PWD}/logistic_regression_output/${chrom}:${start}-${end}_casecontrol.tsv exists"
    fi
    IDX=$((IDX + 1))
done

# bsub -q short -P re_gecip_neurology -J "Array[1-23152:10]" -R rusage[mem=48000] -M 48000 -n 6 -o $PWD/jobs/job%I.out -e $PWD/jobs/job%I.err $PWD/logistic_regression.sh
# a=list.files(pattern='tsv')
# data=do.call('rbind',lapply(a,function(x){out=read.table(x,header=T,stringsAsFactors=F);return(c(x,out[,2]))}))
# data[,12]=as.numeric(data[,12])
# data=data[order(data[,12],decreasing=T),]
# head(data)
# head(data[grep(pattern='100k',data[,1]),])
# colnames(data)=c('locus','case_genomes','control_genomes','stat','pvalue','stat_normalised','pvalue_normalised','top_z_score','top_z_score_normalised','lm_slope','lm_pvalue','specificity_x_tp','specificity','sensitivity','total_normalised_reads')
# data_gms=data[grep(pattern='gms',data$locus),]
# data_100k=data[grep(pattern='100k',data$locus),]
# data_gms$locus=gsub('gms.tsv','',data_gms$locus)
# data_100k$locus=gsub('100k.tsv','',data_100k$locus)
# data_merge=merge(data_100k,data_gms,by='locus')
# data_merge=data_merge[data_merge$pvalue_normalised.x<0.05,]
# data_merge=data_merge[data_merge$pvalue_normalised.y<0.05,]
# data_merge$score=data_merge$specificity_x_tp.x+data_merge$specificity_x_tp.y
# data_merge=data_merge[order(data_merge$score,decreasing=T),]
# write.table(data_merge[,c('locus','pvalue_normalised.x','pvalue_normalised.y','specificity_x_tp.x','specificity_x_tp.y')])



