library(data.table)
library(ggplot2)
library(dplyr)
args = commandArgs(trailingOnly=TRUE)

data=fread(args[1],sep='\t',fill=T)

data$platekey=gsub("/weka/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg38/","",data$bam_file)
data$platekey=gsub("/weka/re_gecip/shared_allGeCIPs/AD_VGD/HTT_bigjson_allconfig/Final_hg37/","",data$platekey)
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


data$first_second_gt[which(read_counts[,1]=="1")]=foo$X2[which(read_counts[,1]=="1")]
data$first_second_gt[which(read_counts[,2]=="1")]=foo$X1[which(read_counts[,2]=="1")]

# data$read_counts[which(read_counts[,1]=="1")]=read_counts[which(read_counts[,1]=="1"),2]
# data$read_counts[which(read_counts[,2]=="1")]=read_counts[which(read_counts[,2]=="1"),1]


# data$first_second_gt[which(read_counts[,1]=="2")]=foo$X2[which(read_counts[,1]=="2")]
# data$first_second_gt[which(read_counts[,2]=="2")]=foo$X1[which(read_counts[,2]=="2")]

# data$read_counts[which(read_counts[,1]=="2")]=read_counts[which(read_counts[,1]=="2"),2]
# data$read_counts[which(read_counts[,2]=="2")]=read_counts[which(read_counts[,2]=="2"),1]
# data[twos_idx,]=data_twos


ratios=as.integer(read_counts[,1])/as.integer(read_counts[,2])
data$first_second_gt[which(ratios>4)]=foo$X1[which(ratios>4)]
data$read_counts[which(ratios>4)]=read_counts[which(ratios>4),1]

tab=table(data$first_second_gt)
counts=data.frame()
for(i in names(tab)){
    count=tab[names(tab)==i]
    if(grepl(pattern='/',i)){
        split=strsplit(i,"/")[[1]]
        for(s in split){
            if(s%in%counts$GT){counts$count[counts$GT==s]=as.numeric(counts$count[counts$GT==s])+count}else{counts=rbind(counts,cbind(GT=s,count=count))}
        }
    }else{
        s=i
        if(s%in%counts$GT){counts$count[counts$GT==s]=as.numeric(counts$count[counts$GT==s])+count*2}else{counts=rbind(counts,cbind(GT=s,count=as.numeric(count)*2))}
    }
    counts$count=as.numeric(counts$count)
}
total=nrow(data)*2
counts$freq=counts$count/total
counts=counts[order(counts$count,decreasing=T),]
write.table(counts,paste0('alleles_counted_',args[2]),sep="\t",quote=F,row.names=F)

write.table(data,args[2],sep="\t",row.names=F,col.names=T,quote=F)

