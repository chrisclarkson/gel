library(data.table)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
f="gms_and_100k_no_relations_casecontrol_locus_min5_annotated.tsv"
f="gms_and_100k_no_relations_europeans_casecontrol_locus_min5_annotated.tsv"
data=fread(f,sep='\t',header=T)
data=data.frame(data,stringsAsFactors=F)

manifest=read.table('manifest_merge_gms',sep='\t',header=F,stringsAsFactors=F)
cases=manifest[manifest$V2=='case',]
controls=manifest[manifest$V2=='control',]

patient_contained=function(x,cases){
    counts=lapply(x,function(x){
        pats=unlist(strsplit(x,','))
        pats=gsub("\\:\\d+.\\d+$","",pats)
        count=sum(pats%in%cases)/length(pats)
        return(count)
    })
    return(unlist(counts))
}

read_count_mean=function(x){
    counts=lapply(x,function(x){
        pats=unlist(strsplit(x,','))
        pats=gsub("^.*\\:","",pats)
        count=mean(as.numeric(pats))
        return(count)
    })
    return(unlist(counts))
}

patient_numbers=function(x){
    counts=lapply(x,function(x){
        pats=unlist(strsplit(x,','))
        return(length(pats))
    })
    return(unlist(counts))
}

get_patients=function(x){
    patients=lapply(x,function(x){
        pats=unlist(strsplit(x,','))
        pats=gsub("\\:\\d+$","",pats)
        return(pats)
    })
    return(unlist(patients))
}

counts=patient_contained(data$counts,cases$V1)
read_counts=read_count_mean(data$counts)
lengths=patient_numbers(data$counts)

counts=data.frame(percentage_pats=counts*100,lengths=lengths,read_counts,stringsAsFactors=F)
data=data.frame(data,counts,stringsAsFactors=F)
data[(data$percentage_pats>90 & data$lengths>=3),colnames(data)!='counts']

counts_summary=counts%>%group_by(percentage_pats,lengths)%>%count()
counts_summary$percentage_pats=round(counts_summary$percentage_pats,1)
data$percentage_pats=round(data$percentage_pats,1)
datag=makeGRangesFromDataFrame(data[,c('contig','start','end','percentage_pats','lengths')],
    seqnames.field='contig',start.field='start',end.field='end',keep.extra.columns=T)
datag$percentage_pats=round(datag$percentage_pats)
ataxia_genes=makeGRangesFromDataFrame(
    read.table('../Chris_Clarkson/clinical_data/known_genes.bed',sep='\t',header=F,stringsAsFactors=F),
    seqnames.field='V1',start.field='V2',end.field='V3',keep.extra.columns=T)
enrich=function(i){
    datag=datag[(datag$percentage_pats==counts_summary$percentage_pats[i] & datag$lengths==counts_summary$lengths[i]),]
    return(length(findOverlaps(datag,ataxia_genes)@from)
    )
}
overlaps=unlist(lapply(1:nrow(counts_summary),enrich))
counts_summary$overlaps=overlaps

pdf('perc_cases_vs_number_case_control_gms_100k_coloured.pdf',width=8)
ggplot(counts_summary,aes(x=percentage_pats,y=lengths,col=as.character(overlaps),size=n))+
    geom_point()+
    scale_y_log10()+
    scale_color_manual(values=c("0"="#000000","1"="#6e0000","2"="#a00000","3"="#db0000","4"="#ff0000","5"="orange","6"="yellow"))+
    xlab('% people holding a given repeat that were cases')+
    ylab('Number of people with repeat')+
    theme(axis.title=element_text(size=20))+labs(col="Number of overlaps\nwith known RE loci",size='Number of loci')+ 
    guides(color = guide_legend(override.aes = list(size=5)))
dev.off()
head(data[((!grepl(pattern='dist',data$gene) & nchar(data$motif)%in%c(3,4,5,6,9))),colnames(data)!='counts'])
