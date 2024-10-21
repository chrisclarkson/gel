library(data.table)
library(ggplot2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
data=fread(args[1],sep='\t',fill=T)

data$platekey=gsub("/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg38/","",data$bam_file)
data$platekey=gsub("/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg37/","",data$platekey)
data$platekey=basename(data$platekey)
data$platekey=gsub("_HTT_realigned.bam","",data$platekey)
data$platekey=gsub("_realigned.bam","",data$platekey)
genotypes=sapply(data$first_second_gt,function(x){
    split=strsplit(x,"/")[[1]]
    if(length(split)==1){out=NA}else{out=split[2]}
    return(out)
})

foo <- data.frame(do.call('rbind', strsplit(as.character(data$first_second_gt),'/',fixed=TRUE)))

read_counts=do.call('rbind',lapply(data$read_counts,function(x){
    split=strsplit(x,"/")[[1]]
    if(length(split)==1){out=c(split,NA)}else{out=split}
    return(out)
}))

out=cbind(data$platekey,foo)
colnames(out)=c('platekey','RC_GT1','RC_GT2')
write.table(out,args[2],sep='\t',row.names=F,quote=F)
