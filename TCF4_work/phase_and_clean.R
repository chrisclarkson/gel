source('phase_and_clean_and_plot_tools.R')

data=unique(data.frame(fread('TCF4_sample_0-3_eh_both_motifs.tsv',sep='\t'), stringsAsFactors=F))
data_0=clean_data(data,4,c('first_second_gt','GAA1','GAG1','AAG','GAG2','CAG1','CAC','CAG2'))
data_phased=phase_keep_read_counts(data_0,
list(c('GAA1','GAG1','AAG','GAG2'),c('CAC','CAG1','CAG2')),
c('EH_GAG','EH_CAG'),multiples=NULL,'CAG')

# Filtering based on matching integer values
# data_phased[which(as.integer(data_phased$CAG_sum_a1)==as.integer(data_phased$EH_CAG_A1)),c('GT1','CAG_sum_a1')]=c(NA,NA)
# data_phased[which(as.integer(data_phased$CAG_sum_a2)==as.integer(data_phased$EH_CAG_A2)),c('GT2','CAG_sum_a2')]=c(NA,NA)
# data_phased[which(as.integer(data_phased$GAG_sum_a1)==as.integer(data_phased$EH_GAG_A1)),c('GT1','GAG_sum_a1')]=c(NA,NA)
# data_phased[which(as.integer(data_phased$GAG_sum_a2)==as.integer(data_phased$EH_GAG_A2)),c('GT2','GAG_sum_a2')]=c(NA,NA)

data_phased$bam_file=basename(gsub('.expansionhunter_realigned.bam','',data_phased$bam_file))

data_0_2=unique(data.frame(fread('TCF4_sample_0-2.tsv',sep='\t'), stringsAsFactors=F))
data_0_2$bam_file=basename(gsub('.expansionhunter_realigned.bam','',data_0_2$bam_file))
data_0_2=clean_data(data_0_2,4,c('first_second_gt'))
data_0_2=format_v0_data(data_0_2)

data_alleles_imputed=merge_and_phase(data_phased,data_0_2,c('GT1','GT2','CAG_sum_a1','CAG_sum_a2','EH_CAG_A1','EH_CAG_A2'),'TCF4_phased.txt')
# nrow(data_alleles_imputed)

data_alleles_imputed=data_alleles_imputed[!(data_alleles_imputed$GT1_imputed=='imputed' & data_alleles_imputed$GT2_imputed=='imputed'),]
# nrow(data_alleles_imputed)

# g_tcf4=melt_columns_and_plot_structure_across_ancestries_and_undone(
# data.frame(data_alleles_imputed, stringsAsFactors=F),
# c('GT1','GT2'), NULL, c('1','0','1','2'),
# file_suffix='.expansionhunter_realigned.bam',
# get_rid=NULL, facet_ancestries=T, get_rid_of_numbers=F,
# name='TCF4', full_list_of_platekeys=data_0_2$bam_file
# )


revc=function(x){
    opposite_strand=function(x){
        x=gsub('A','t',x);x=gsub('G','c',x);x=gsub('T','a',x);x=gsub('C','g',x)
        x=toupper(x)
        return(x)
    }
    reverse_string=function(x){
        return(sapply(lapply(strsplit(x, NULL), rev), paste, collapse=""))
    }
    x=opposite_strand(x)
    x=reverse_string(x)
    return(x)
}

revc_cols=revc(c('GAA1','GAG1','AAG','AGG2'))
revc_data=data_alleles_imputed[,c('GAA1_A1','GAG1_A1','AAG_A1','AGG2_A1','GAA1_A2','GAG1_A2','AAG_A2','AGG2_A2')]
colnames(revc_data)=c(paste0(revc_cols,'_A1'),paste0(revc_cols,'_A2'))

revc_data$bam_file=data_alleles_imputed$bam_file
revc_data$GT1=data_alleles_imputed$GT1
revc_data$GT2=data_alleles_imputed$GT2
revc_data$new_structure_a1=''
revc_data$new_structure_a2=''

for( r in revc_cols){
    revc_data$new_structure_a1=paste0(revc_data$new_structure_a1,r,'x',as.integer(revc_data[,paste0(r,'_A1')]),'|')
    revc_data$new_structure_a2=paste0(revc_data$new_structure_a2,r,'x',as.integer(revc_data[,paste0(r,'_A2')]),'|')
}

revc_data$new_structure_a1[grep(pattern='CAC',data_alleles_imputed$GT1)]=paste(revc_data$new_structure_a1[grep(pattern='CAC',data_alleles_imputed$GT1)],'GTG interruption present')
revc_data$new_structure_a2[grep(pattern='CAC',data_alleles_imputed$GT2)]=paste(revc_data$new_structure_a2[grep(pattern='CAC',data_alleles_imputed$GT2)],'GTG interruption present')
print(head(revc_data))
revc_data=revc_data[,c('bam_file','new_structure_a1','new_structure_a1','Q1_1','Q1_2')]
colnames(revc_data)=c('bam_file','new_structure_a1','new_structure_a1','EH_CAG_A1','EH_CAG_A2')
write.table(revc_data,'alleles_imputed_and_formatted.tsv',sep='\t',row.names=F,quote=F)


