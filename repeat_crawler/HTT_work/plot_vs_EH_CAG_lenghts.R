library(data.table)
library(ggplot2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)
# data=fread("most_likely_gts_parents_last2_excluded_eh_annotated.tsv",sep='\t',fill=T)
# data=fread("most_likely_gts_parents_cag_included_cleaned_eh_annotated.tsv",sep='\t',fill=T)
data=fread(args[1],sep='\t',fill=T)

# data$platekey=gsub("/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg38/","",data$bam_file)
# data$platekey=gsub("/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg37/","",data$platekey)
# data$platekey=gsub("_HTT_realigned.bam","",data$platekey)
# data$eh_read_counts=info$Spanning

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
    # print(read_counts)
    # if(!all(read_counts<2) | all(is.na(read_counts))){
    #     split[read_counts<2]=NA
    #     alleles[read_counts<2]=NA
    # }
    split=split[!is.na(split)]
    alleles=alleles[!is.na(alleles)]
    if(length(split)==1){
        out=c(split,NA,ratio,alleles,NA)
    }else{
        out=c(min(as.numeric(split),na.rm=T),
            max(as.numeric(split),na.rm=T),ratio,
            alleles[which.min(split)],
            alleles[which.max(split)])
    }
    return(out)
}))
# cag_stretches=cag_stretches[,-ncol(cag_stretches)]
cag_stretches=data.frame(cag_stretches,stringsAsFactors=F)
colnames(cag_stretches)=c('X1','X2','ratio','A1','A2')
cag_stretches$cat="RC"
cag_eh_stretches=do.call('rbind',lapply(data$EH_rep_no,function(x){
    split=strsplit(x,"/")[[1]]
    if(length(split)==1){out=c(split,NA)}else{out=c(min(as.numeric(split)),max(as.numeric(split)))}
    return(out)
}))
cag_eh_stretches=data.frame(cag_eh_stretches,stringsAsFactors=F)
cag_eh_stretches$cat="EH"
idx=which(as.integer(cag_stretches$X1)!=as.integer(cag_eh_stretches$X1) & data$amend=='normal')
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

cag_stretches[idx,]=cag_stretches_sub

homo_hetero=ifelse(cag_eh_stretches$X1==cag_eh_stretches$X2,"homo","hetero")
cags=rbind(cbind(data$bam_file,cag_stretches$X1,cag_eh_stretches$X1,cag_stretch_counts[,1],'A1',homo_hetero,cag_stretches$X3,cag_stretches$A1),
           cbind(data$bam_file,cag_stretches$X2,cag_eh_stretches$X2,cag_stretch_counts[,2],'A2','hetero',cag_stretches$X3,cag_stretches$A2)
)
cags=data.frame(cags,stringsAsFactors=F)
colnames(cags)=c("bam_file","RC","EH",'read_count',"cat","homo_hetero",'allele')
cags$RC=as.numeric(as.character(cags$RC))
cags$EH=as.numeric(as.character(cags$EH))
cags$read_count=as.numeric(as.character(cags$read_count))

# cags$ratio=as.numeric(as.character(cags$ratio))
# cags=cags[cags$allele!='(CAG|CCT)',]
cags=cags[!is.na(cags$allele),]
counts=cags%>%group_by(RC,EH,cat,homo_hetero,allele)%>%
    count()

g=ggplot(counts,aes(x=RC,y=EH,size=n,col=allele))+
    geom_point()+
    geom_abline()+xlim(0,35)+ylim(0,35)+theme(axis.title.x=element_text(size=20),axis.title.y=element_text(size=20))+
    facet_grid(allele~cat)

png(args[2],height=1020,width=1440)
print(g)
dev.off()

