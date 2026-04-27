library(data.table)
library(ggplot2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
data=fread(args[1],sep='\t',fill=T)

data$platekey=gsub("/weka/re_gecip/.snapshots/2023070701b01/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg38/","",data$bam_file)
data$platekey=gsub("/weka/re_gecip/.snapshots/2023070701b01/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg37/","",data$platekey)
data$platekey=gsub("_HTT_realigned.bam","",data$platekey)

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

cag_stretch_counts=t(apply(data,1,function(x){
    split=as.numeric(strsplit(x[4],"/")[[1]])
    read_counts=as.numeric(strsplit(x[3],"/")[[1]])
    if(length(split)==1){
        out=c(read_counts,NA)
    }else{
        idx=c(which.min(split),which.max(split))
        out=read_counts[idx]
    }
    return(out)
}))
colnames(cag_stretch_counts)=c('A1_read_count','A2_read_count')

alleles=t(apply(data,1,function(x){
    split=as.numeric(strsplit(x[4],"/")[[1]])
    alleles=strsplit(x[2],"/")[[1]]
    if(length(split)==1){
        out=c(alleles,NA)
    }else{
        idx=c(which.min(split),which.max(split))
        out=alleles[idx]
    }
    return(out)
}))
colnames(alleles)=c('A1','A2')

cag_stretches=do.call('rbind',lapply(1:nrow(data),function(i){
    x=data$CAG[i]
    read_counts=as.numeric(strsplit(data$read_counts[i],"/")[[1]])
    split=strsplit(x,"/")[[1]]
    ratio=read_counts[which.min(as.numeric(split))]/read_counts[which.max(as.numeric(split))]
    alleles=strsplit(data$first_second_gt[i],"/")[[1]]
    split=split[!is.na(split)]
    alleles=alleles[!is.na(alleles)]
    if(length(split)==1){
        out=c(split,NA,ratio,alleles,NA)
    }else{

        if(split[1]!=split[2]){
        out=c(min(as.numeric(split),na.rm=T),
            max(as.numeric(split),na.rm=T),ratio,
            alleles[which.min(split)],
            alleles[which.max(split)])
        }else{
            out=c(min(as.numeric(split),na.rm=T),
            max(as.numeric(split),na.rm=T),ratio,
            alleles[1],
            alleles[2])
        }
    }
    return(out)
}))
# cag_stretches=cag_stretches[,-ncol(cag_stretches)]
cag_stretches=data.frame(cag_stretches,stringsAsFactors=F)
# print('here')
# print(dim(data))
# print(dim(cag_stretches))
# print(cag_stretches[10720,])
colnames(cag_stretches)=c('X1','X2','ratio','A1','A2')
cag_stretches$cat="RC"

cag_eh_stretches=do.call('rbind',lapply(data$EH_rep_no,function(x){
    split=strsplit(x,"/")[[1]]
    if(length(split)==1){out=c(split,NA)}else{out=c(min(as.numeric(split)),max(as.numeric(split)))}
    return(out)
}))
cag_eh_stretches=data.frame(cag_eh_stretches,stringsAsFactors=F)
cag_eh_stretches$cat="EH"
idx=which(as.integer(cag_stretches$X1)!=as.integer(cag_eh_stretches$X1))
cag_stretches_sub=cag_stretches[idx,]
cag_eh_stretches_sub=cag_eh_stretches[idx,]

idx2=which(as.integer(cag_stretches_sub$X1)==as.integer(cag_eh_stretches_sub$X2))
idx3=which(as.integer(cag_stretches_sub$X2)==as.integer(cag_eh_stretches_sub$X1))
cag_stretches_sub$X2[idx2]=cag_stretches_sub$X1[idx2]
# cag_stretches_sub$X1[idx2]=0
cag_stretches_sub$A2[idx2]=cag_stretches_sub$A1[idx2]
# cag_stretches_sub$A1[idx2]=NA
cag_stretches_sub$X1[idx3]=cag_stretches_sub$X2[idx3]
# cag_stretches_sub$X1[idx3]=0
cag_stretches_sub$A1[idx3]=cag_stretches_sub$A2[idx3]
# cag_stretches_sub$A1[idx3]=NA

# cag_stretches[idx,]=cag_stretches_sub


out=cbind(data$platekey,cag_stretches,cag_eh_stretches)
colnames(out)=c('bam_file','RC_A1','RC_A2','ratio','RC_GT1','RC_GT2','cat','EH_A1','EH_A2')
out=out[,c('bam_file','EH_A1','EH_A2','RC_A1','RC_A2','RC_GT1','RC_GT2')]

write.table(out,args[2],sep='\t',row.names=F,quote=F)

