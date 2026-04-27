library(data.table)
library(ggplot2)
library(grid)
library(gridExtra)
library(ggsignif)
library(ggpubr)
library(cowplot)
library(dplyr)
library(scales)

count_alleles=function(data){
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
    return(counts)
}


clean_data=function(data,thresh,cols){
    ratios_idx=apply(data,1,function(x){
        split=as.numeric(strsplit(x[4],"/")[[1]])
        read_counts=as.numeric(strsplit(x['read_counts'],"/")[[1]])
        if(length(split)==1){
            read_counts=c(read_counts,NA)
        }
        if(!is.na(read_counts[2])){
        if((read_counts[1]/read_counts[2])>=thresh){
            out=1
        }else if((read_counts[2]/read_counts[1])>=thresh){
            out=2
        }else{
            out=0
        }
        }else{out=1}
        return(out)
    })
    data_sub=data[ratios_idx!=0,]
    ratios_idx_sub=ratios_idx[ratios_idx!=0]
    to=nrow(data_sub)
    for(col in cols){
        replacement=unlist(lapply(1:to,function(x){
            favorite=ratios_idx_sub[x]
            x=data_sub[x,col]
            return(strsplit(x,"/")[[1]][favorite])
        }))
        data[ratios_idx!=0,col]=replacement
    }
    return(data)
}

get_counts=function(data,cols){
    to=nrow(data)*2
    counts_mat=1:to
    for(col in c(cols)){
        counts=do.call('rbind',lapply(data[,col],function(x){
            split=strsplit(x,"/")[[1]]
            if(length(split)==1){out=c(split,NA)}else{out=split}
            return(out)
        }))
        if(col=='first_second_gt'){
            counts=c(counts[,1],counts[,2])
        }else{
            counts=as.numeric(c(counts[,1],counts[,2]))
        }
        counts_mat=cbind(counts_mat,counts)
    }
    counts_mat=counts_mat[,-1]
    colnames(counts_mat)=cols
    counts_mat=data.frame(counts_mat,stringsAsFactors=F)
    counts_mat$bam_file=c(data$bam_file,data$bam_file)
    return(counts_mat)
}

split_values <- function(x) {
  if (grepl("/", x)) {
    parts <- strsplit(x, "/")[[1]]
  } else {
    parts=c(x, NA) # If no "/", fill with NAs
  }
  return(parts)
}

rearrange_cols=function(data,col1,col2,idx){
    tmp=cbind(data[,c(col1,col2)],idx)
    tmp=t(apply(tmp,1,function(x){x[c(col1,col2)]=x[c(col1,col2)][as.numeric(c(x[["RC_idx1"]],x[["RC_idx2"]]))];return(x)}))
    return(tmp[,c(col1,col2)])
}
phase=function(data,rc_cols,eh_cols,multiples=NA,phase_according_to){
    multiple_dataframe=data.frame()
    if(!is.null(multiples)){
        for(m in multiples){
            m=unlist(strsplit(m,':'))
            m1=m[1]
            m2=m[2]
            multiple_dataframe=rbind(multiple_dataframe,cbind(m1,m2))
        }
        multiple_dataframe=data.frame(multiple_dataframe,stringsAsFactors=F)
        multiple_dataframe$m2=as.numeric(multiple_dataframe$m2)
    }
    out_final=1:nrow(data)
    for(i in 1:length(rc_cols)){
        cols_out=c()
        out=1:nrow(data)
        cols=rc_cols[i][[1]]
        for(col in cols){
            new_cols=do.call('rbind',lapply(data[,col],split_values))
            out=cbind(out,new_cols)
            cols_out=c(cols_out,paste0(col,'_A1'),paste0(col,'_A2'))
        }
        if(i>1){
            out=data.frame(out[,-1],stringsAsFactors=F)
            colnames(out)=cols_out
        }else{
            out=data.frame(out,stringsAsFactors=F)
            colnames(out)[-1]=cols_out
        }
        if(!is.null(multiples)){
            for(m in multiple_dataframe$m1){
                if(sum(colnames(out)==paste0(m,'_A1'))==1){
                    out[,colnames(out)==paste0(m,'_A1')]=as.numeric(out[,colnames(out)==paste0(m,'_A1')])*multiple_dataframe$m2[multiple_dataframe$m1==m]
                    out[,colnames(out)==paste0(m,'_A2')]=as.numeric(out[,colnames(out)==paste0(m,'_A2')])*multiple_dataframe$m2[multiple_dataframe$m1==m]
                }
            }}
        if(length(cols)>1){
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a1')]=colSums(apply(out[,paste0(cols,'_A1')],1,as.numeric))
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a2')]=colSums(apply(out[,paste0(cols,'_A2')],1,as.numeric))
        }else{
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a1')]=out[,paste0(cols,'_A1')]
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a2')]=out[,paste0(cols,'_A2')]
        }
        eh_lengths=apply(do.call('rbind',sapply(data[,eh_cols[i]],function(x){x=strsplit(x,"/");if(length(unlist(x))!=2){x=list(c(NA,NA))};return(x)})),2,as.numeric)
        colnames(eh_lengths)=c(paste0(eh_cols[i],'_A1'),paste0(eh_cols[i],'_A2'))
        if(!('GT1' %in% colnames(out_final))){
            repeat_structures=do.call('rbind',sapply(data[,'first_second_gt'],function(x){return(strsplit(x,"/"))}))
            colnames(repeat_structures)=c('GT1','GT2')
        }
        if(i==1){
            out_final=cbind(cbind(out_final,out),repeat_structures,eh_lengths)
        }else{
            out_final=cbind(cbind(out_final,out),eh_lengths)
        }
    }
    out_final=data.frame(out_final[,!(colnames(out_final)%in%c('out_final','out','out.1','out.2','GT1.1','GT2.1'))],stringsAsFactors=F)
    RC_idx1=apply(out_final[,c(paste0(phase_according_to,'_sum_a1'),paste0(phase_according_to,'_sum_a2'))],1,which.min)
    RC_idx2=apply(out_final[,c(paste0(phase_according_to,'_sum_a1'),paste0(phase_according_to,'_sum_a2'))],1,which.max)
    idxs=cbind(RC_idx1,RC_idx2)
    EH_data=out_final[,grep(pattern='^EH_',colnames(out_final))]
    sum_data=out_final[,grep(pattern='_sum_',colnames(out_final))]
    data_for_phasing=out_final[,-grep(pattern='^EH_',colnames(out_final))]
    data_for_phasing=data_for_phasing[,-grep(pattern='_sum_',colnames(data_for_phasing))]
    structure_data=data_for_phasing[,c('GT1','GT2')]
    data_for_phasing=data_for_phasing[,!(colnames(data_for_phasing) %in% c('GT1','GT2'))]
    phase_cols=unique(gsub("_A\\d","",colnames(data_for_phasing)))
    for(col in phase_cols){
        data_for_phasing[,c(paste0(col,'_A1'),paste0(col,'_A2'))]=rearrange_cols(data_for_phasing,paste0(col,'_A1'),paste0(col,'_A2'),idxs)
    }
    structure_data=rearrange_cols(structure_data,'GT1','GT2',idxs)
    for(col in eh_cols){
        EH_data[,paste0(col,'_A1')]=apply(EH_data[,c(paste0(col,'_A1'),paste0(col,'_A2'))],1,min)
        EH_data[,paste0(col,'_A2')]=apply(EH_data[,c(paste0(col,'_A1'),paste0(col,'_A2'))],1,max)
        sum_data[,c(paste0(gsub('EH_','',col),'_sum_a1'),paste0(gsub('EH_','',col),'_sum_a2'))]=rearrange_cols(sum_data,paste0(gsub('EH_','',col),'_sum_a1'),paste0(gsub('EH_','',col),'_sum_a2'),idxs)
    }
    
    # out_final$EH_A1=apply(out_final[,c('EH1','EH2')],1,min)
    # out_final$EH_A2=apply(out_final[,c('EH1','EH2')],1,max)
    # return(out[,c('RC_A1','RC_A2','EH_A1','EH_A2')])
    
    out_final=data.frame(cbind(bam_file=data$bam_file,data_for_phasing,structure_data,sum_data,EH_data),stringsAsFactors=F)
    return(out_final)
}


phase_keep_read_counts=function(data,rc_cols,eh_cols,multiples=NA,phase_according_to){
    print(phase_according_to)
    read_counts=do.call('rbind',lapply(data$read_counts,split_values))
    read_counts=data.frame(read_counts)
    read_counts$X1=as.numeric(read_counts$X1)
    read_counts$X2=as.numeric(read_counts$X2)
    colnames(read_counts)=c('read_count1','read_count2')
    # print(head(read_counts))
    multiple_dataframe=data.frame()
    if(!is.null(multiples)){
        for(m in multiples){
            m=unlist(strsplit(m,':'))
            m1=m[1]
            m2=m[2]
            multiple_dataframe=rbind(multiple_dataframe,cbind(m1,m2))
        }
        multiple_dataframe=data.frame(multiple_dataframe,stringsAsFactors=F)
        multiple_dataframe$m2=as.numeric(multiple_dataframe$m2)
    }
    out_final=1:nrow(data)
    for(i in 1:length(rc_cols)){
        cols_out=c()
        out=1:nrow(data)
        cols=rc_cols[i][[1]]
        for(col in cols){
            new_cols=do.call('rbind',lapply(data[,col],split_values))
            out=cbind(out,new_cols)
            cols_out=c(cols_out,paste0(col,'_A1'),paste0(col,'_A2'))
        }
        if(i>1){
            out=data.frame(out[,-1],stringsAsFactors=F)
            colnames(out)=cols_out
        }else{
            out=data.frame(out,stringsAsFactors=F)
            colnames(out)[-1]=cols_out
        }
        if(!is.null(multiples)){
            for(m in multiple_dataframe$m1){
                if(sum(colnames(out)==paste0(m,'_A1'))==1){
                    out[,colnames(out)==paste0(m,'_A1')]=as.numeric(out[,colnames(out)==paste0(m,'_A1')])*multiple_dataframe$m2[multiple_dataframe$m1==m]
                    out[,colnames(out)==paste0(m,'_A2')]=as.numeric(out[,colnames(out)==paste0(m,'_A2')])*multiple_dataframe$m2[multiple_dataframe$m1==m]
                }
            }}
        if(length(cols)>1){
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a1')]=colSums(apply(out[,paste0(cols,'_A1')],1,as.numeric))
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a2')]=colSums(apply(out[,paste0(cols,'_A2')],1,as.numeric))
        }else{
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a1')]=out[,paste0(cols,'_A1')]
            out[,paste0(gsub('EH_','',eh_cols[i]),'_sum_a2')]=out[,paste0(cols,'_A2')]
        }
        eh_lengths=apply(do.call('rbind',sapply(data[,eh_cols[i]],function(x){x=strsplit(x,"/");if(length(unlist(x))!=2){x=list(c(NA,NA))};return(x)})),2,as.numeric)
        colnames(eh_lengths)=c(paste0(eh_cols[i],'_A1'),paste0(eh_cols[i],'_A2'))
        if(!('GT1' %in% colnames(out_final))){
            repeat_structures=do.call('rbind',sapply(data[,'first_second_gt'],function(x){return(strsplit(x,"/"))}))
            colnames(repeat_structures)=c('GT1','GT2')
        }
        if(i==1){
            out_final=cbind(cbind(out_final,out),repeat_structures,eh_lengths)
        }else{
            out_final=cbind(cbind(out_final,out),eh_lengths)
        }
    }
    out_final=data.frame(out_final[,!(colnames(out_final)%in%c('out_final','out','out.1','out.2','GT1.1','GT2.1'))],stringsAsFactors=F)
    RC_idx1=apply(out_final[,c(paste0(phase_according_to,'_sum_a1'),paste0(phase_according_to,'_sum_a2'))],1,which.min)
    RC_idx2=apply(out_final[,c(paste0(phase_according_to,'_sum_a1'),paste0(phase_according_to,'_sum_a2'))],1,which.max)
    idxs=cbind(RC_idx1,RC_idx2)
    EH_data=out_final[,grep(pattern='^EH_',colnames(out_final))]
    sum_data=out_final[,grep(pattern='_sum_',colnames(out_final))]
    data_for_phasing=out_final[,-grep(pattern='^EH_',colnames(out_final))]
    data_for_phasing=data_for_phasing[,-grep(pattern='_sum_',colnames(data_for_phasing))]
    structure_data=data_for_phasing[,c('GT1','GT2')]
    data_for_phasing=data_for_phasing[,!(colnames(data_for_phasing) %in% c('GT1','GT2'))]
    phase_cols=unique(gsub("_A\\d","",colnames(data_for_phasing)))
    for(col in phase_cols){
        data_for_phasing[,c(paste0(col,'_A1'),paste0(col,'_A2'))]=rearrange_cols(data_for_phasing,paste0(col,'_A1'),paste0(col,'_A2'),idxs)
    }
    structure_data=rearrange_cols(structure_data,'GT1','GT2',idxs)
    for(col in eh_cols){
        EH_data[,paste0(col,'_A1')]=apply(EH_data[,c(paste0(col,'_A1'),paste0(col,'_A2'))],1,min)
        EH_data[,paste0(col,'_A2')]=apply(EH_data[,c(paste0(col,'_A1'),paste0(col,'_A2'))],1,max)
        sum_data[,c(paste0(gsub('EH_','',col),'_sum_a1'),paste0(gsub('EH_','',col),'_sum_a2'))]=rearrange_cols(sum_data,paste0(gsub('EH_','',col),'_sum_a1'),paste0(gsub('EH_','',col),'_sum_a2'),idxs)
    }
    # print('oi')
    read_counts[,c('read_count1','read_count2')]=rearrange_cols(read_counts,'read_count1','read_count2',idxs)
    out_final=data.frame(cbind(bam_file=data$bam_file,data_for_phasing,structure_data,sum_data,EH_data,read_counts),stringsAsFactors=F)
    return(out_final)
}

plot_ratios=function(data_phased,rc_cols,eh_cols,thresh){
    # idx=ifelse(data_phased$read_count1>data_phased$read_count2,0,1)
    cols=c(rc_cols,eh_cols)
    for(col in cols){
        data_phased=data_phased[!is.na(data_phased[,col]),]
    }
    numbers=data.frame(do.call('rbind',lapply(1:nrow(data_phased),function(i){
        # print(data_phased[i,])
        if(data_phased$read_count2[i]>data_phased$read_count2[i]){
            ratio=data_phased$read_count1[i]/data_phased$read_count2[i]
            read_count=data_phased$read_count1[i]
            if(data_phased[i,rc_cols[2]]==data_phased[i,eh_cols[2]]){
                correct_or_not='correct'
            }else{
                correct_or_not='incorrect'
            }
        }else{
            ratio=data_phased$read_count2[i]/data_phased$read_count1[i]
            read_count=data_phased$read_count1[i]
            if(data_phased[i,rc_cols[1]]==data_phased[i,eh_cols[1]]){
                correct_or_not='correct'
            }else{
                correct_or_not='incorrect'
            }
        }
        return(c(read_count,ratio,correct_or_not))
    })),stringsAsFactors=F)
    colnames(numbers)=c('read_count','ratio','correct_or_not')
    # print(head(numbers))
    g_ratios=ggplot(numbers,aes(x=as.numeric(as.character(ratio)),y=as.numeric(read_count),col=correct_or_not))+
        geom_point()+geom_vline(xintercept=thresh,linetype='dashed')+
        theme(plot.title = element_text(hjust = 0.5))
    return(g_ratios)
}

switcharoo=function(data_phased,phase_to){
    idx1=apply(data_phased,1,function(x){
        as.numeric(x[paste0(phase_to,"_sum_a1")])%in%as.numeric(c(x[paste0("EH_",phase_to,"_A1")],x[paste0("EH_",phase_to,"_A2")]))
    })
    idx2=apply(data_phased,1,function(x){
        as.numeric(x[paste0(phase_to,"_sum_a2")])%in%as.numeric(c(x[paste0("EH_",phase_to,"_A1")],x[paste0("EH_",phase_to,"_A2")]))
    })
    data_phased$cond=paste0(idx1,'_',idx2)
    # data_phased=data_phased[data_phased$cond!='FALSE_FALSE',]
    true_trues=data_phased[data_phased$cond=='TRUE_TRUE',]
    true_falses=data_phased[data_phased$cond=='TRUE_FALSE',]
    false_trues=data_phased[data_phased$cond=='FALSE_TRUE',]
    true_falses[true_falses[,paste0(phase_to,"_sum_a1")]!=true_falses[,paste0("EH_",phase_to,"_A1")],c('GT1','GT2',paste0(phase_to,"_sum_a1"),paste0(phase_to,"_sum_a2"))]=true_falses[true_falses[,paste0(phase_to,"_sum_a1")]!=true_falses[,paste0("EH_",phase_to,"_A1")],c('GT2','GT1',paste0(phase_to,"_sum_a2"),paste0(phase_to,"_sum_a1"))]
    false_trues[false_trues[,paste0(phase_to,"_sum_a2")]!=false_trues[,paste0("EH_",phase_to,"_A2")],c('GT1','GT2',paste0(phase_to,"_sum_a1"),paste0(phase_to,"_sum_a2"))]=false_trues[false_trues[,paste0(phase_to,"_sum_a2")]!=false_trues[,paste0("EH_",phase_to,"_A2")],c('GT2','GT1',paste0(phase_to,"_sum_a2"),paste0(phase_to,"_sum_a1"))]
    return(rbind(true_trues,true_falses,false_trues))
}

melt_columns_and_plot_structure_heatmap=function(data,struct_cols,rc_cols,
                                    components_cols,
                                    component_names,eh_cols,
                                    get_rid=NULL,allele_n_threshold=100000,
                                    name,accepted_structures=NULL,sample_gts=NULL){
    out=data.frame()
    for(i in 1:length(struct_cols)){
        component_cols=components_cols[i][[1]]
        d_comps=data[,component_cols]
        colnames(d_comps)=c(component_names)
        if(!is.null(rc_cols)){
            
            d=cbind(data[,struct_cols[i]],data[,rc_cols[i]],d_comps,data[,eh_cols[i]])
        }else{
            d=cbind(data[,struct_cols[i]],d_comps,data[,eh_cols[i]])
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    if(!is.null(rc_cols)){
        colnames(out)=c('structure','rc_rep',component_names,'repeat_size')
        out=data.frame(out,stringsAsFactors=F)
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
    }else{
        colnames(out)=c('structure',component_names,'repeat_size')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    if(!is.null(accepted_structures)){
        out=out[out$structure%in%accepted_structures,]
    }
    out$structure_annotated=NA
    tab=table(out$structure)
    for(t in names(tab)){
        n=tab[names(tab)==t]
        if(n<allele_n_threshold){
            out=out[-which(out$structure==t),]
        }else{
            out$structure_annotated[out$structure==t]=paste0(t,' N=',n)
        }
    }
    f=paste0('structure_vs_eh_repsize_',name,'_heatmap.pdf')
    file.remove(f)
    pdf(f,height=20,width=10)
    if(!is.null(sample_gts)){
        for(gt in sample_gts){
            df_sub=out[out$structure==gt,]
            out=out[out$structure!=gt,]
            df_sub=df_sub[sample(round(nrow(df_sub)*0.2)),]
            df_sub$structure_annotated=paste(df_sub$structure_annotated,'(20% sampled)')
            out=rbind(out,df_sub)
        }
    }else{
        out=out[1:65535,]
    }
    structure_annotated=out$structure_annotated
    plot_out=apply(as.matrix(data.frame(out,stringsAsFactors=F)[,component_names]),2,as.numeric)
    colnames(plot_out)=component_names
    print(head(plot_out))
    col_fun<-colorRamp2(c(1,2,3,4,5,8,15,25,40,41), c(RColorBrewer::brewer.pal(9,'Reds'),'yellow'))
    draw(Heatmap(plot_out,
            name='Component repeat length',
            cluster_rows=T,
            # cluster_row_slices=F,
            # column_split=row_splits_filtered,
            # col=col_fun,
            # row_names_gp = gpar(fontsize=40),
            # column_names_gp = gpar(fontsize=10),
            row_title_rot=0,
            border='black',
            column_names_rot=70,
            column_title_rot=90,
            row_split=structure_annotated,
            cluster_columns=F,
            cluster_column_slices=F,
            # column_title_side = "bottom",
            show_row_names=F,
            show_column_names=T,
            gap=unit(5,'mm'),
            column_gap = unit(5, "mm"),
            # row_title_gp = gpar(fontsize = 40)
            heatmap_legend_param = list(title='Component repeat length',
                at = c(1,2,3,4,5,8,15,25,40,41), 
                labels = c('1','2','3','4','5','8','15','25','40','>40'),
                labels_gp=gpar(fontsize=40),
                title_gp=gpar(fontsize=40),
                legend_height=unit(10,'cm'),
                legend_width = unit(20, "cm")),
            column_title_gp = gpar(fontsize = 40)
            )
            # heatmap_legend_side="top", 
            # padding = unit(c(t = 10, r = 2, b = 2, l = 11), "cm")
            )

    dev.off()
    print(f)

}
melt_columns_and_plot_structure=function(data,struct_cols,rc_cols,eh_cols,get_rid=NULL,allele_n_threshold=0,name,do_comparisons=NULL){
    out=data.frame()
    for(i in 1:length(struct_cols)){
        if(!is.null(rc_cols)){
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,rc_cols[i]],data[,eh_cols[i]])
        }else{
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,eh_cols[i]])
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    if(!is.null(rc_cols)){
        colnames(out)=c('bam_file','structure','rc_rep','repeat_size')
        out=data.frame(out,stringsAsFactors=F)
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
    }else{
        colnames(out)=c('bam_file','structure','repeat_size')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    out$structure_annotated=NA
    tab=table(out$structure)
    # print(tab)
    for(t in names(tab)){
        n=tab[names(tab)==t]
        if(n<allele_n_threshold){
            out=out[-which(out$structure==t),]
        }else{
            out$structure_annotated[out$structure==t]=paste0(t,' N=',n)
        }        
    }
    print(head(out))
    g=ggplot(na.omit(out),aes(x=structure_annotated,y=as.numeric(as.character(repeat_size))))+
        geom_violin()+
        geom_boxplot(width=0.1)+
        theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab('Repeat structure')+
        ylab(paste('Expansion Hunter Repeat size'))
    if(!is.null(do_comparisons)){
        comparisons=na.omit(unique(out$structure_annotated))
        comparisons=t(combn(comparisons,2))
        combinations=list()
        for(i in 1:nrow(comparisons)){
            combinations[[i]]=c(comparisons[i,])
        }
        print(combinations)
        g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1)
    }
    f=paste0('structure_vs_eh_repsize_',name,'.png')
    file.remove(f)
    png(f,height=1400,width=1000)
    print(g)
    dev.off()
    print(f)
    # return(g)
    return(list(sizes=out,plot=g))
}

bind_ancestries=function(data, file_suffix='_realigned.bam'){
    data$bam_file=basename(data$bam_file)
    data$bam_file=gsub('EH_','',data$bam_file)
    data$bam_file=gsub(file_suffix,'',data$bam_file)
    population_meta=fread('/re_gecip/neurology/Kristina/POPULATION/list_34191_unrel_notNeuro_withSuperpopu.tsv',sep='\t',header=T)
    data=merge(data,population_meta,by.x=c('bam_file'),by.y=c('platekey'),all.x=T)
    data$superpopu[is.na(data$superpopu)]='Unknown'
    data$superpopu[data$superpopu=='AFR']='African'
    data$superpopu[data$superpopu=='EUR']='European'
    data$superpopu[data$superpopu=='AMR']='American'
    data$superpopu[data$superpopu=='EAS']='East Asian'
    data$superpopu[data$superpopu=='SAS']='South Asian'
    return(data)
}

melt_columns_and_plot_structure_across_ancestries=function(data,struct_cols,rc_cols,
                                        eh_cols,get_rid=NULL,allele_n_threshold=0,
                                        name,do_comparisons=NULL,file_suffix='_realigned.bam',
                                        percentage_threshold=0,allele_order=NULL){
    out=data.frame()
    data=bind_ancestries(data,file_suffix)
    for(i in 1:length(struct_cols)){
        if(!is.null(rc_cols)){
            d=cbind(data[,struct_cols[i]],data[,rc_cols[i]],data[,eh_cols[i]],data$superpopu)
        }else{
            d=cbind(data[,struct_cols[i]],data[,eh_cols[i]],data$superpopu)
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    if(!is.null(rc_cols)){
        colnames(out)=c('structure','rc_rep','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
    }else{
        colnames(out)=c('structure','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    percentage=data.frame(out,stringsAsFactors=F) %>%
        group_by(ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100,
            upperE=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100,
            lowerE=(-1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100)
    g_freq_stacked=ggplot(percentage[percentage$freq>percentage_threshold,],aes(x=ancestry,y=freq,fill=structure))+
        geom_bar(stat='identity')+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        xlab('Population')+
        ylab('Percentages')
    out$structure_annotated=NA
    out_new=data.frame()
    for(t in unique(out$structure)){
        out_sub=out[out$structure==t,]
        # if(length(unique(out_sub$ancestry))>1){
            for(ancestry in unique(out_sub$ancestry)){
                out_sub2=out_sub[out_sub$ancestry==ancestry,]
                n=nrow(out_sub2)
                if(n>0){
                    out_sub2=out_sub2[out_sub2$structure==t,]
                    # if(n<allele_n_threshold){
                    #     out_sub2=out_sub2[-which(out_sub2$structure==t),]
                    # }else{
                        out_sub2$structure_annotated[out_sub2$structure==t]=paste0(t,' N=',n)
                    # }
                    out_new=rbind(out_new,out_sub2)
                }
            }
        # }
    }
    out=out_new
    if(!is.null(percentage_threshold)){
        percentage=percentage[(percentage$freq)>=percentage_threshold,]
        out_new=data.frame()
        for(ancestry in unique(percentage$ancestry)){
            percentage_sub=percentage[percentage$ancestry==ancestry,]
            out_sub=out[out$ancestry==ancestry,]
            out_sub=out_sub[out_sub$structure%in%unique(percentage_sub$structure),]
            out_new=r69366
        out=out_new
    }
    g=ggplot(na.omit(out),aes(y=ancestry,x=as.numeric(as.character(repeat_size))))+
        geom_violin()+
        geom_boxplot()+
        theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        ylab('')+
        xlab(paste('Expansion Hunter Repeat size'))+
        facet_wrap(~structure,ncol=1,scales='free_y')+
        theme(axis.text.y=element_blank())
    percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    g_freq=ggplot(na.omit(percentage),aes(y=ancestry_annotated,x=freq))+
        geom_bar(stat='identity',col='black',fill='white')+
        geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
        ylab('Repeat structure')+
        xlab('Percentages')+
        facet_wrap(~factor(structure),ncol=1,scales='free')+theme_cowplot()
    g=ggarrange(g_freq,g,nrow=1,ncol=2,widths=c(1,2),align='h')
    if(!is.null(do_comparisons)){
        comparisons=na.omit(unique(out$structure))
        comparisons=t(combn(comparisons,2))
        combinations=list()
        for(i in 1:nrow(comparisons)){
            combinations[[i]]=c(comparisons[i,])
        }
        print(combinations)
        g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1)
    }
    f=paste0('structure_vs_eh_repsize_ancestries_',name,'.png')
    file.remove(f)
    png(f,height=1000,width=800)
    print(g)
    dev.off()
    print(f)
    return(list(sizes=out,proportions=percentage,plot=g,plot_proportions=g_freq_stacked))
}}

melt_columns_and_plot_structure_across_ancestries2=function(data,struct_cols,rc_cols,
                                        eh_cols,get_rid=NULL,allele_n_threshold=0,
                                        name,do_comparisons=NULL,file_suffix='_realigned.bam',
                                        percentage_threshold=0,allele_order=NULL,
                                        pathogenic_repeats=NULL,facet_ancestries=FALSE,
                                        title=NULL,get_rid_of_numbers=TRUE,phase_to=NULL){
    if(!is.null(phase_to)){
        data=switcharoo(data,phase_to)
    }
    out=data.frame()
    data=bind_ancestries(data,file_suffix)
    print(head(data))
    for(i in 1:length(struct_cols)){
        if(!is.null(rc_cols)){
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,rc_cols[i]],data[,eh_cols[i]],data$superpopu)
        }else{
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,eh_cols[i]],data$superpopu)
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    # print(head(out))
    if(!is.null(rc_cols)){
        colnames(out)=c('bam_file','structure','rc_rep','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
        # print('hereeee')
        # print(dim(out))
        # print(head(out[as.integer(out$rc_rep)!=as.integer(out$repeat_size),]))
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
        print(dim(out))
    }else{
        # print(nrow(out))
        colnames(out)=c('bam_file','structure','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    out=na.omit(out)
    percentage=data.frame(out,stringsAsFactors=F) %>%
        group_by(ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100,
            upperE=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100,
            lowerE=(-1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100)
    percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    if(get_rid_of_numbers){percentage$structure=gsub("(?<!x)\\d","",percentage$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        percentage$structure[((percentage$freq)<percentage_threshold) & !(percentage$structure%in%pathogenic_repeats)]='rare'
        percentage$structure[(percentage$structure%in%pathogenic_repeats)]=paste('Pathogenic:',percentage$structure[(percentage$structure%in%pathogenic_repeats)])
    }else{
        a=data.frame(percentage[(percentage$freq)>=percentage_threshold,])
        print(a[!(a$structure%in%allele_order),])
        percentage$structure[(percentage$freq)<percentage_threshold]='rare'
    }
    # percentage=na.omit(percentage)
    out$structure_annotated=NA
    out$ancestry_annotated=NA
    out_new=data.frame()
    for(t in unique(out$structure)){
        out_sub=out[out$structure==t,]
        # if(length(unique(out_sub$ancestry))>1){
            for(ancestry in unique(out_sub$ancestry)){
                out_sub2=out_sub[out_sub$ancestry==ancestry,]
                # if(ancestry!='rare'){
                    n=nrow(out_sub2)
                # }else{
                #     n=nrow(out[out$ancestry==ancestry,])
                # }
                if(n>0){
                    out_sub2=out_sub2[out_sub2$structure==t,]
                    if(n<allele_n_threshold){
                        out_sub2=out_sub2[-which(out_sub2$structure==t),]
                    }else{
                        out_sub2$structure_annotated[out_sub2$structure==t]=paste0(t,' N=',n)
                        out_sub2$ancestry_annotated=paste0(ancestry,' N=',n)
                    }
                    out_new=rbind(out_new,out_sub2)
                }
            }
        # }
    }
    out=out_new
    if(get_rid_of_numbers){out$structure=gsub("(?<!x)\\d","",out$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        pathogenic_sub=out[out$structure%in%pathogenic_repeats,]
        pathogenic_sub$structure=paste('Pathogenic:',pathogenic_sub$structure)
        pathogenic_sub$structure_annotated=paste('Pathogenic:',pathogenic_sub$structure_annotated)
        out=out[!(out$structure%in%pathogenic_repeats),]
        # print(nrow(out))
    }
    if(!is.null(percentage_threshold)){
        percentage2=percentage[percentage$freq>=percentage_threshold,]
        out_new=data.frame()
        for(ancestry in unique(percentage2$ancestry)){
            percentage_sub=percentage2[percentage2$ancestry==ancestry,]
            out_sub=out[out$ancestry==ancestry,]
            out_sub$structure[!(out_sub$structure%in%unique(percentage_sub$structure))]='rare'
            out_sub$ancestry_annotated[out_sub$structure=='rare']=paste0(ancestry,' N=',sum(out_sub$structure=='rare'))
            out_new=rbind(out_new,out_sub)
        }
        out=out_new
    }
    if(!is.null(pathogenic_repeats)){
        out=rbind(out,pathogenic_sub)
        print(nrow(out))
    }
    if(!is.null(allele_order)){
        if(get_rid_of_numbers){allele_order=gsub("(?<!x)\\d","",allele_order,perl=T)}
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }else{
        allele_order=unique(out$structure)
        # print('s')
        print(allele_order)
        if(('rare'%in%allele_order)){
            allele_order=allele_order[allele_order!='rare']
            allele_order=c(allele_order,'rare')
        }
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731')
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    if(facet_ancestries){
        g=ggplot(na.omit(out),aes(y=ancestry_annotated,x=as.numeric(as.character(repeat_size)),fill=structure,group=ancestry_annotated))+
            geom_violin()+
            geom_boxplot(width=0.1)+
            theme_cowplot()+
            theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))+
            facet_wrap(~structure,ncol=1,scales='free_y')
    }else{
        g=ggplot(na.omit(out),aes(y=factor(structure,levels=rev(allele_order)),x=as.numeric(as.character(repeat_size)),fill=structure))+
            geom_violin()+
            geom_boxplot(width=0.1)+
            theme_cowplot()+
            theme(legend.key.size=unit(2, "cm"),axis.title.y=element_text(size=20),axis.title.x=element_text(size=25),axis.text.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))
        for(a in allele_order){
            print(a)
            nalleles=nrow(out[as.character(out$structure)==a,])
            print(nalleles)
            if(nalleles>0){
                g=g+annotate(geom='text',y=a, x=round(max(as.numeric(as.character(out$repeat_size)))*0.9,0),
                            size=8, label= paste0('N=',nalleles),vjust=1.6,hjust='inwards')
            }
            # print(nalleles)
        }
    }
    g=g+fillScale
    if(!is.null(do_comparisons)){
        print('compare')
        if(!facet_ancestries){
            comparisons=na.omit(unique(as.character(out$structure)))
            comparisons=t(combn(comparisons,2))
            combinations=list()
            for(i in 1:nrow(comparisons)){
                combinations[[i]]=c(comparisons[i,])
            }
        }else{
            comparisons=na.omit(unique(as.character(out$ancestry_annotated)))
            comparisons=t(combn(comparisons,2))
            combinations=list()
            for(i in 1:nrow(comparisons)){
                combinations[[i]]=c(comparisons[i,])
            }
        }
        print(combinations)

        g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1,orientation='y')
        # g=g+stat_compare_means(method = "wilcox.test", ref.group='.all.')
    }
    # g=g+theme(legend.key.height=unit(0.2,'cm'),legend.key.width=unit(0.2,'cm'),legend.height=element_text(size=5),legend.text=element_text(size=5))
    # print(percentage[is.na(percentage$ancestry),])
    # print(percentage[is.na(percentage$structure),])
    # percentage$structure[is.na(percentage$structure)]='rare'
    g_freq_stacked=ggplot(percentage,aes(y=ancestry,x=freq,fill=structure))+
        geom_bar(stat='identity')+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        ylab('Population')+
        xlab('Percentages')+theme_cowplot()+
            theme(legend.position='none',axis.title.y=element_text(size=10),axis.title.x=element_text(size=25),axis.text.y = element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20))
    g_freq_stacked=g_freq_stacked+fillScale
    # percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    # g_freq=ggplot(na.omit(percentage),aes(y=ancestry_annotated,x=freq))+
    #     geom_bar(stat='identity',col='black',fill='white')+
    #     geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
    #     ylab('Repeat structure')+
    #     xlab('Percentages')+
    #     facet_wrap(~factor(structure),ncol=1,scales='free')+theme_cowplot()
    
    g1=ggarrange(g_freq_stacked+ggtitle(name)+theme(plot.title=element_text(size=20)),g,nrow=2,ncol=1,heights=c(1,3))
    # if(!is.null(title)){
    #     g=g+ggtitle(title)
    # }
    
    f=paste0('structure_vs_eh_repsize_ancestries_',name,'.pdf')
    file.remove(f)
    pdf(f,height=20,width=15)
    print(g1)
    dev.off()
    print(f)
    return(list(sizes=out,proportions=percentage,plot=g,plot_proportions=g_freq_stacked))
}


compare_cis_modifiers=function(data_new,out_name,facet_ancestries=F){
    data_new$category_annotated=data_new$category
    for(category in unique(data_new$category)){
        n=sum(data_new$category==category,na.rm=T)
        print(n)
        data_new$category_annotated[data_new$category==category]=paste0(data_new$category_annotated[data_new$category==category],' (N=',n,')')
    }
    print(head(data_new))
    # comparisons=na.omit(unique(data_new$category_annotated))
    # comparisons=t(combn(comparisons,2))
    # combinations=list()
    # for(i in 1:nrow(comparisons)){
    #     combinations[[i]]=c(comparisons[i,])
    # }
    g=ggplot(data_new,aes(x=category_annotated,y=as.numeric(repeat_size)))+
        geom_violin()+geom_boxplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Cis-modifiers')+ylab('Expansion Hunter repeat size')+geom_smooth(aes(group=1),method='lm')
    # g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1)
    mod=summary(lm(as.numeric(data_new$repeat_size)~as.numeric(factor(data_new$category_annotated))))
    print(mod)
    slope=round(mod$coefficients[2,1],1)
    pval=mod$coefficients[2,4]
    if(grepl(pattern='e',as.character(pval))){
        split=unlist(strsplit(as.character(pval),'e'))
        pval=paste0(round(as.numeric(split[1]),1),split[2])
    }
    g=g+annotate(geom='text',x=levels(factor(data_new$category_annotated))[1], y=max(as.numeric(as.character(data_new$repeat_size)))-5,size=6, 
            label=paste0('Slope=',slope,'\nP-value=',pval),vjust=0,hjust=0)
    f=paste0(out_name,'.png')
    print(f)
    png(f)
    print(g)
    dev.off()
    return(g)
}


plot_pairs=function(data_mat,rc_cols,rc_compare_cols,eh_cols,out){
    plot_list=list()
    counter=1
    for(i in 1:length(rc_cols)){
        cols=rc_cols[[i]]
        print(cols)
        if(!is.null(cols)){
            name=gsub("_A[0-9]+","",cols[1])
            name=gsub("_a[0-9]+","",name)
            print(name)
            # print(length(data_mat[,c(cols[1])]))
            rc_data=rbind(cbind(name,data_mat[,c(cols[1])]),cbind(name,data_mat[,c(cols[2])]))
            cols_eh=eh_cols[i][[1]]
            eh_data=c(data_mat[,cols_eh[1]],data_mat[,cols_eh[2]])
            cols_rc_compare=rc_compare_cols[[i]]
            # print('working?')
            # print(length(data_mat[,c(cols_rc_compare[1])]))
            rc_compare_data=c(data_mat[,c(cols_rc_compare[1])],data_mat[,c(cols_rc_compare[2])])
            # print(dim(rc_data))
            # print(length(eh_data))
            # print(length(rc_compare_data))
            data_plot=cbind(rc_data,eh_data,rc_compare_data)
            data_plot[,c(2,3,4)]=apply(data_plot[,c(2,3,4)],2,as.numeric)
            colnames(data_plot)=c('V1','V2','V3','V4')
            data_plot=data.frame(data_plot)
            data_plot=data_plot[data_plot$V3==data_plot$V4,]
            data_plot=na.omit(data_plot)
            # print(nrow(data_plot))
            counts=data.frame(table(data_plot$V2))
            counts$Var1=as.integer(as.character(counts$Var1))
            counts=counts[order(counts$Var1,decreasing=F),]
            if(nrow(counts)>2){
                more_than_two=TRUE
            }else{
                more_than_two=FALSE
            }
            if(more_than_two){
                g=ggplot(data_plot,aes(x=as.integer(V2),y=as.numeric(V3)))+
                    xlab(name)+ylab(gsub("EH_","Expansion Hunter ",gsub("_A[0-9]+","",eh_cols[i][[1]])))+
                    geom_violin(aes(group = cut_width(V2, 1)),scale = "width")+geom_point()+
                    geom_smooth(method='lm')
                mod=summary(lm(as.numeric(data_plot$V3)~as.numeric(as.character(data_plot$V2))))
                slope=round(mod$coefficients[2,1],1)
                pval=mod$coefficients[2,4]
                if(grepl(pattern='e',as.character(pval))){
                    split=unlist(strsplit(as.character(pval),'e'))
                    pval=paste0(round(as.numeric(split[1]),1),"10^",split[2])
                }else{
                    pval=round(pval,2)
                }
                g=g+annotate(geom='text',x=min(as.numeric(as.character(data_plot$V2))), y=max(as.numeric(as.character(data_plot$V3)))-5,size=6, 
                        label=paste0('Slope=',slope,'\nP-value=',pval),vjust=0,hjust=0)
            }else{
                g=ggplot(data_plot,aes(x=as.integer(V2),y=as.numeric(V3)))+
                    xlab(name)+ylab(gsub("EH_","Expansion Hunter ",gsub("_A[0-9]+","",eh_cols[i][[1]])))+
                    geom_violin(aes(group = cut_width(V2, 1)),scale = "width")+geom_point()
                mod=wilcox.test(as.numeric(as.character(data_plot$V3))[data_plot$V2==unique(data_plot$V2)[1]],
                                as.numeric(as.character(data_plot$V3))[data_plot$V2==unique(data_plot$V2)[2]])
                # print(mod)
                pval=mod$p.value
                if(grepl(pattern='e',as.character(pval))){
                    split=unlist(strsplit(as.character(pval),'e'))
                    pval=paste0(round(as.numeric(split[1]),1),"10^",split[2])
                }else{
                    pval=round(pval,2)
                }
                # print('here')
                # print(pval)
                x=min(as.numeric(as.character(data_plot$V2)))
                y=max(as.numeric(as.character(data_plot$V3)))-5
                # print(x)
                # print(y)
                g=g+annotate(geom='text',x=x,
                    y=y,size=6,
                    label=paste0('Wilcoxon Rank Sum\nP-value=',as.character(pval)),
                    vjust=0,hjust=0)
            }
        }else{
            g=ggplot()
        }
        g=g+
        theme(axis.text.x=element_text(size=30),axis.text.y=element_text(size=30),
              axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))
        if(more_than_two){
            # print(counts)
            to=nrow(counts)-1
            for(i in 1:to){
                difference=counts$Var1[i+1]-counts$Var1[i]
                # print(difference)
                if(difference>1){
                    difference=difference-1
                    for(j in difference){
                        counts=rbind(counts,c(counts$Var1[i]+j,0))
                    }
                }
            }
        }
        g_bar=ggplot(counts,aes(x=Var1,y=Freq))+
            geom_col()+
            scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
              labels = trans_format("log10", math_format(10^.x)))+
            xlab('')+
            ylab('Count')+
            theme_cowplot()+
            theme(axis.text.x=element_text(size=30),
            # axis.text.y=element_text(size=14),
                 axis.title.x=element_text(size=30),axis.title.y=element_text(size=30))
        # g=ggarrange(g_bar,g,ncol=1,nrow=2,align='v',heights=c(1,4))
        plot_list[[counter]]=g
        counter=counter+1
    }
    f=paste0('violin_plots_phased_',out,'.png')
    file.remove(f)
    print(f)
    png(f,width=1000,height=1000)
    do.call('grid.arrange',c(plot_list))
    dev.off()
    return(plot_list)
}


plot_pairs2=function(data_mat,x_axis_phasers,rc_cols,rc_compare_cols,eh_cols,out,thresh=0){
    plot_list=list()
    counter=1
    for(i in 1:length(rc_cols)){
        cols=rc_cols[[i]]
        print(cols)
        if(!is.null(cols)){
            name=gsub("_A[0-9]+","",cols[1])
            name=gsub("_a[0-9]+","",name)
	    name=gsub("_sum","",name)
            print(name)
            phaser_cols=x_axis_phasers[[i]]
            x_axis_phaser_data=c(data_mat[,c(phaser_cols[1])],data_mat[,c(phaser_cols[2])])
            rc_data=rbind(cbind(name,data_mat[,c(cols[1])]),cbind(name,data_mat[,c(cols[2])]))
            cols_eh=eh_cols[i][[1]]
            eh_data=c(data_mat[,cols_eh[1]],data_mat[,cols_eh[2]])
            cols_rc_compare=rc_compare_cols[[i]]
            rc_compare_data=c(data_mat[,c(cols_rc_compare[1])],data_mat[,c(cols_rc_compare[2])])
            data_plot=cbind(rc_data,eh_data,rc_compare_data,x_axis_phaser_data)
            data_plot[,c(2,3,4)]=apply(data_plot[,c(2,3,4)],2,as.numeric)
            colnames(data_plot)=c('V1','V2','V3','V4','V5')
            data_plot=data.frame(data_plot)
            print(head(data_plot))
            data_plot=data_plot[as.numeric(data_plot$V3)==as.numeric(data_plot$V4),]
            data_plot=data_plot[as.numeric(data_plot$V2)==as.numeric(data_plot$V5),]
            data_plot=na.omit(data_plot)
            print(head(data_plot))
            counts=data.frame(table(data_plot$V2))
            print(counts)
            counts$Var1=as.integer(as.character(counts$Var1))
            counts=counts[order(counts$Var1,decreasing=F),]
            counts=counts[as.numeric(counts$Freq)>thresh,]
            print(counts)
            data_plot=data_plot[data_plot$V2%in%counts$Var1,]
            if(nrow(counts)>2){
                more_than_two=TRUE
            }else{
                more_than_two=FALSE
            }
            data_plot$n_annotated=NA
            annotations_levels=c()
            for(n in counts$Var1){
                print(n)
                annotation=paste0(n,' N=',sum(data_plot$V2==n,na.rm=T))
                data_plot$n_annotated[data_plot$V2==n]=annotation
                annotations_levels=c(annotations_levels,annotation)
            }
            print('here')
            data_plot$n_annotated=factor(data_plot$n_annotated,levels=annotations_levels)
            print(levels(data_plot$n_annotated))
            if(more_than_two){
                g=ggplot(data_plot,aes(x=n_annotated,y=as.numeric(V3)))+
                    xlab(name)+ylab(gsub("EH_","Expansion Hunter ",gsub("_A[0-9]+","",eh_cols[i][[1]])))+
                    geom_violin(aes(group = cut_width(V2, 1)),scale = "width")+
                    geom_boxplot(aes(group = cut_width(V2, 1)),scale = "width")+
                    geom_smooth(aes(group=1),method='lm')
                mod=summary(lm(as.numeric(data_plot$V3)~as.numeric(as.character(data_plot$V2))))
                slope=round(mod$coefficients[2,1],1)
                pval=mod$coefficients[2,4]
                if(grepl(pattern='e',as.character(pval))){
                    split=unlist(strsplit(as.character(pval),'e'))
                    pval=paste0(round(as.numeric(split[1]),1),"10^",split[2])
                }else{
                    pval=round(pval,2)
                }
                g=g+annotate(geom='text',x=levels(data_plot$n_annotated)[1], y=max(as.numeric(as.character(data_plot$V3)))-5,size=6, 
                        label=paste0('Slope=',slope,'\nP-value=',pval),vjust=0,hjust=0)
            }else{
                g=ggplot(data_plot,aes(x=n_annotated,y=as.numeric(V3)))+
                    xlab(name)+ylab(gsub("EH_","Expansion Hunter ",gsub("_A[0-9]+","",eh_cols[i][[1]])))+
                    geom_violin(aes(group = cut_width(V2, 1)),scale = "width")+
                    geom_boxplot(aes(group = cut_width(V2, 1)),scale = "width")
                mod=wilcox.test(as.numeric(as.character(data_plot$V3))[data_plot$V2==unique(data_plot$V2)[1]],
                                as.numeric(as.character(data_plot$V3))[data_plot$V2==unique(data_plot$V2)[2]])
                pval=mod$p.value
                if(grepl(pattern='e',as.character(pval))){
                    split=unlist(strsplit(as.character(pval),'e'))
                    pval=paste0(round(as.numeric(split[1]),1),"10^",split[2])
                }else{
                    pval=round(pval,2)
                }
                x=min(as.numeric(as.character(data_plot$V2)))
                y=max(as.numeric(as.character(data_plot$V3)))-5
                g=g+annotate(geom='text',x=x,
                    y=y,size=6,
                    label=paste0('Wilcoxon Rank Sum\nP-value=',as.character(pval)),
                    vjust=0,hjust=0)
            }
        }else{
            g=ggplot()
        }
        g=g+theme_cowplot()+theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1, size=15),
                 axis.text.y=element_text(size=15),
                 axis.title.x=element_text(size=15),
                 axis.title.y=element_text(size=15))
        # g=ggarrange(g_bar,g,ncol=1,nrow=2,align='v',heights=c(1,4))
        plot_list[[counter]]=g
        counter=counter+1
    }
    f=paste0('violin_plots_phased_',out,'.png')
    file.remove(f)
    print(f)
    png(f,width=1000,height=1000)
    do.call('grid.arrange',c(plot_list))
    dev.off()
    return(plot_list)
}

format_v0_data=function(data,extra_cols=NULL){
    genotypes=sapply(data$first_second_gt,function(x){
        split=strsplit(x,"/")[[1]]
        if(length(split)==1){out=NA}else{out=split[2]}
        return(out)
    })
    foo <- data.frame(do.call('rbind', strsplit(as.character(data$first_second_gt),'/',fixed=TRUE)))
    out=cbind(data$bam_file,foo)
    colnames(out)=c('bam_file','GT1_v0','GT2_v0')
    if(!is.null(extra_cols)){
        for(col in extra_cols){
            foo <- data.frame(do.call('rbind', strsplit(as.character(data[,col]),'/',fixed=TRUE)))
            cols_out=c(colnames(out),paste0(col,'1'),paste0(col,'2'))
            out=cbind(out,foo)
            colnames(out)=cols_out
        }
    }
    return(out)
}

merge_and_phase=function(data_v2,data_v0,input_cols,out_file){
    translate_cols_v2=c('RC_GT1','RC_GT2','RC_A1','RC_A2','EH_A1','EH_A2')
    for(idx in 1:length(translate_cols_v2)){
        colnames(data_v2)[colnames(data_v2)==input_cols[idx]]=translate_cols_v2[idx]
    }
    translate_cols_v0=c('GT1_v0','GT2_v0')
    for(idx in 1:length(translate_cols_v0)){
        colnames(data_v0)[colnames(data_v0)==input_cols[idx]]=translate_cols_v0[idx]
    }
    eha1_rca1=data_v2[data_v2$EH_A1 == data_v2$RC_A1,]
    eha2_rca2=data_v2[data_v2$EH_A2==data_v2$RC_A2,]
    eha1_rca2=data_v2[data_v2$EH_A1 == data_v2$RC_A2,]
    eha2_rca1=data_v2[data_v2$EH_A2 == data_v2$RC_A1,]
    eha1_rca2=eha1_rca2[!(eha1_rca2$bam_file%in%eha1_rca1$bam_file),]
    eha1_rca2=eha1_rca2[!is.na(eha1_rca2$bam_file),]
    eha2_rca2=eha2_rca2[!(eha2_rca2$bam_file%in%eha1_rca1$bam_file) & !(eha2_rca2$bam_file%in%eha1_rca2$bam_file),]
    eha2_rca2=eha2_rca2[!is.na(eha2_rca2$bam_file),]
    eha2_rca1=eha2_rca1[!(eha2_rca1$bam_file%in%eha1_rca1$bam_file) & !(eha2_rca1$bam_file%in%eha1_rca2$bam_file) & !(eha2_rca1$bam_file %in% eha2_rca2$bam_file),]
    eha2_rca1=eha2_rca1[!is.na(eha2_rca1$bam_file),]
    eha1_rca1$Q1_1<-eha1_rca1$RC_A1
    eha1_rca1$GT1<-eha1_rca1$RC_GT1
    eha1_rca1$Q1_2<-eha1_rca1$EH_A2
    eha1_rca2$Q1_1<-eha1_rca2$RC_A2
    eha1_rca2$GT1<-eha1_rca2$RC_GT2
    eha1_rca2$Q1_2<-eha1_rca2$EH_A2
    eha2_rca1$Q1_1<-eha2_rca1$RC_A1
    eha2_rca1$GT1<-eha2_rca1$RC_GT1
    eha2_rca1$Q1_2<-eha2_rca1$EH_A1
    eha2_rca2$Q1_1<-eha2_rca2$RC_A2
    eha2_rca2$GT1<-eha2_rca2$RC_GT2
    eha2_rca2$Q1_2<-eha2_rca2$EH_A1
    data_v2_sorted<-rbind(eha1_rca1, eha1_rca2, eha2_rca1,eha2_rca2)
    data_v2_sorted$GT1_imputed[data_v2_sorted$EH_A1!=data_v2_sorted$RC_A1]='imputed'
    data_v2_sorted$GT2_imputed[data_v2_sorted$EH_A2!=data_v2_sorted$RC_A2]='imputed'
    data_v2_sorted$GT1_imputed[is.na(data_v2_sorted$GT1_imputed)]='not imputed'
    data_v2_sorted$GT2_imputed[is.na(data_v2_sorted$GT2_imputed)]='not imputed'
    # put the allele types into categories a1
    # add v0 genotypes 
    data_v2_v0<-merge(data_v2_sorted, data_v0,by = "bam_file")
    # print(head(data_v2_v0))
    # Checking which v0 genotypes match the phased genotype
    matching_g1=data_v2_v0[as.character(data_v2_v0$GT1)==as.character(data_v2_v0$GT1_v0),]
    matching_g2=data_v2_v0[as.character(data_v2_v0$GT1)==as.character(data_v2_v0$GT2_v0),]
    matching_g2=matching_g2[!(matching_g2$bam_file%in%matching_g1$bam_file),]
    # Assigning the second unphased genotype
    matching_g1$GT2=matching_g1$GT2_v0
    matching_g2$GT2=matching_g2$GT1_v0
    out=rbind(matching_g1,matching_g2)
    # print(paste('difference',nrow(data_v2_v0)-nrow(out)))
    unmatched<-data_v2_v0[!data_v2_v0$bam_file %in%out$bam_file,]
    write.table(unique(out), file=out_file,sep = "\t", quote = F, row.names = F)
    write.table(unmatched, file=paste0("unmatched_",out_file),sep = "\t", quote = F, row.names = F)
    return(out)
}


plot_percentage_grid=function(){
    data=data.frame(data,stringsAsFactors=F)
    population_meta=fread('/re_gecip/neurology/Kristina/POPULATION/list_34191_unrel_notNeuro_withSuperpopu.tsv',sep='\t',header=T)
    data=merge(data,population_meta,by.x=c('bam_file'),by.y=c('platekey'))
    plot_data=rbind(cbind(platekey=data$bam_file,category=data$GT1_category,population=data$superpopu,CAG_size=data$Q1_1,allele='A1'),
                    cbind(platekey=data$bam_file,category=data$GT2_category,population=data$superpopu,CAG_size=data$Q1_2,allele='A2'))
    plot_data=data.frame(plot_data,stringsAsFactors=F)
    plot_data=plot_data[!is.na(as.character(plot_data$category)),]
    percentages=data.frame()
    count_percentage_across_size=function(plot_data,population=NULL,return_leg=FALSE){
        plot_data=data.frame(plot_data,stringsAsFactors=F)
        if(population!='All'){
            plot_data=plot_data[plot_data$population==population,]
        }else{
            population='All'
        }
        percentages=data.frame()
        for(cat in unique(plot_data$category)){
            plot_data_sub=plot_data[plot_data$category==cat,]
            plot_data_sub=na.omit(plot_data_sub)
            if(nrow(plot_data_sub)>0){
                percentage=plot_data_sub %>%
                    group_by(CAG_size) %>%
                    summarise(n = n()) %>%
                    mutate(total=sum(n),freq = n / sum(n),
                        upperE=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n),
                        lowerE=(-1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n))
                percentage=data.frame(percentage,stringsAsFactors=F)
                percentage$category=cat
                percentage=percentage[!is.na(as.character(percentage$category)),]
                percentages=rbind(percentages,percentage)
            }
        }
        # write.table(percentages,paste0(population,'_percentages.tsv'),sep='\t',row.names=F,quote=F)
        percentages$population=population
        return(percentages)
    }
    for(p in c('All',"AFR","AMR","EAS","EUR","SAS")){
        print(p)
        percentages=rbind(percentages,count_percentage_across_size(plot_data,p))
    }
    percentages=na.omit(percentages)
    percentages$freq=percentages$freq*100
    percentages$population=factor(percentages$population,levels=c('All',"AFR","AMR","EAS","EUR","SAS"))
    percentages$category=factor(percentages$category,levels=c('canonical','Q2_dup','Q2_loss','partial_Q2_loss','Q2_P1_loss','P1_loss'))
    scale_function=function(x, scale,shift){
        return((x)*scale-shift)
    }
    inv_scale_function=function(x, scale,shift){
        return((x+shift)/scale)
    }
    grid_list=data.frame()
    for(p in c('All',"EUR","AFR","AMR","EAS","SAS")){
        percentages_sub=percentages[(percentages$population==p),]
        row_list=list()
        for(cat in c('canonical','Q2_dup','Q2_loss','partial_Q2_loss','Q2_P1_loss','P1_loss')){
            percentages_sub2=percentages_sub[(percentages_sub$category==cat),]
            # percentages_sub2$category=factor(as.character(percentages_sub2$category),levels=c('canonical','Q2_dup','Q2_loss','partial_Q2_loss','Q2_P1_loss','P1_loss'))
            # percentages_sub2$population=factor(as.character(percentages_sub2$population),levels=c('All',"AFR","AMR","EAS","EUR","SAS"))
            max_first=max(percentages_sub2$freq)
            min_first=min(percentages_sub2$freq)
            max_second=max(percentages_sub2$n)
            min_second=min(percentages_sub2$n)
            scale=(max_second-min_second)/(max_first-min_first)
            shift=min_first-min_second
            grid_list=rbind(grid_list,cbind(population=p,category=cat,scale,shift))
        }
        # row=ggarrange(row_list,nrow=1)
    }
    plot_correctly=function(p,cat,scale,shift,percentages){
        percentages_sub2=percentages[(percentages$population==p) & (percentages$cat==cat),]
        scale=as.numeric(scale)
        shift=as.numeric(shift)
        if(nrow(percentages_sub2)>0){
            if(cat=='canonical'){
                ylabel=paste0(p,'\nPercentage')
            }else{
                ylabel=''
            }
            if(cat=='P1_loss'){
                name='Count'
            }else{
                name=''
            }
            if(p=='All'){
                title=gsub('_',' ',cat)
            }else{
                title=''
            }
            if(p=='SAS'){
                xlabel='CAG stretch size'
            }else{
                xlabel=''
            }
            g_freq=ggplot(percentages_sub2,
                aes(x=as.integer(as.character(CAG_size))))+
                geom_bar(aes(y=freq),width=1,stat='identity')+xlim(6,40)+
                scale_y_continuous(sec.axis = sec_axis(~scale_function(.,scale,shift), name = name))+
                scale_color_manual(values=c("#E69F00"))+
                xlab(xlabel)+ylab(ylabel)+
                theme(legend.position='none',axis.text.y.right = element_text(color = "#E69F00"),axis.title.y.right = element_text(color = "#E69F00"))+
                ggtitle(title)
        }else{
            g_freq=ggplot()
        }
        return(g_freq)
    }
    grid_list2=apply(grid_list,1,function(x){plot_correctly(x[1],x[2],x[3],x[4],percentages)})
    f='grid3.pdf'
    file.remove(f)
    pdf(f,width=13,height=7)
    do.call(grid.arrange,c(grid_list2,ncol=6,nrow=6))
    dev.off()
}

plot_imputed_distibrution=function(data,name){
    data_plot=rbind(cbind(data$EH_A1,data$GT1_imputed),cbind(data$EH_A2,data$GT2_imputed))
    colnames(data_plot)=c('EH','imputed_or_not')
    data_plot=data.frame(data_plot,stringsAsFactors=F)
    data_plot$imputed_or_not[data_plot$imputed_or_not=='not imputed']='determined using spanning reads'
    g=ggplot(data_plot,aes(x=as.numeric(EH),fill=imputed_or_not,col=imputed_or_not))+
        geom_histogram()
    f=paste0(name,'_histogram_imputed_vs_not.png')
    file.remove(f)
    print(f)
    png(f)
    print(g)
    dev.off()
    return(g)
}


gene_diversity=function(data_new){
    n=length(unique(data_new$structure))
    counts=table(data_new$structure)
    probs=counts/n
    gd=(n/(1-n))*(1-sum(probs**2))
}


# gene_diversity=function(data_new){
#     counts=table(data_new$structure)
#     probs=counts/nrow(data_new)
#     print(sum(probs))
#     gd_df=data.frame()
#     for(uniq_length in unique(data_new$repeat_size)){
#         d_sub=data_new[data_new$repeat_size==uniq_length,]
#         n=nrow(d_sub)+1
#         probs_sub=probs[names(probs)%in%d_sub$structure]
#         gd=(1-sum(probs_sub**2))
#         gd_df=rbind(gd_df,cbind(uniq_length,gd,n))
#     }
#     colnames(gd_df)=c('repeat_size','gene_diversity','number')
#     return(gd_df[order(as.numeric(gd_df[,1]),decreasing=F),])
# }

gene_diversity=function(data_new){
    gd_df=data.frame()
    # for(ancestry in unique(data_new$ancestry)){
        for(uniq_length in unique(data_new$repeat_size)){
            d_sub=data_new[data_new$repeat_size==uniq_length,];
            n=nrow(d_sub)+1;
            counts=table(d_sub$structure);
            probs=counts/nrow(d_sub);
            probs_sub=probs[names(probs)%in%d_sub$structure];
            gd=(1-sum(probs_sub**2))*(n/(n-1));
            gd_df=rbind(gd_df,cbind(uniq_length,gd,n));
        }
    # }
    colnames(gd_df)=c('repeat_size','gene_diversity','number')
    gd_df=gd_df[order(as.numeric(gd_df[,1]),decreasing=F),]
    library(gtools)
    gd_df$Q=as.integer(quantcut(as.integer(gd_df$repeat_size),q=200))/2
    return(gd_df)
}



gene_diversity=function(data_new,facet_ancestries=F){
    # library(gtools);
    data_new=data_new[order(as.numeric(data_new$repeat_size),decreasing=F),];
    data_new$repeat_size=as.numeric(data_new$repeat_size)
    # data_new$Q=as.integer(quantcut(as.integer(data_new$repeat_size),q=200))/2;
    gd_df=data.frame();
    percentiles=c(0,0.25,0.50,0.75,0.90,1);
    to=length(percentiles)-1;
    quantiles=quantile(as.numeric(data_new$repeat_size),probs=percentiles);
    if(facet_ancestries){
    for(ancestry in unique(data_new$ancestry)){
        for(i in 1:to){
            percentile=quantiles[i];
            next_percentile=quantiles[i+1];
            d_sub=data_new[(data_new$ancestry==ancestry & data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            # d_sub=data_new[data_new$repeat_size==uniq_length,];
            n=nrow(d_sub)+1;
            if(n!=1){
            counts=table(d_sub$structure);
            probs=counts/nrow(d_sub);
            # probs_sub=probs[names(probs)%in%d_sub$structure];
            gd=(1-sum(probs**2))*(n/(n-1));}else{gd=0}
            most_common_structure=names(probs)[order(probs,decreasing=T)];
            most_common_structure=most_common_structure[1:3];
            gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),as.numeric(gsub('%','',names(quantiles)[i+1])),paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),gd,n,ancestry,paste(most_common_structure,collapse=';')));
        }
    }}else{
        for(i in 1:to){
            percentile=quantiles[i];
            next_percentile=quantiles[i+1];
            d_sub=data_new[(data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            # d_sub=data_new[data_new$repeat_size==uniq_length,];
            n=nrow(d_sub)+1;
            if(n!=1){
            counts=table(d_sub$structure);
            probs=counts/nrow(d_sub);
            # probs_sub=probs[names(probs)%in%d_sub$structure];
            gd=(1-sum(probs**2))*(n/(n-1));
            }else{gd=0};
            most_common_structure=names(probs)[order(probs,decreasing=T)];
            most_common_structure=most_common_structure[1:3];
            gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),
                as.numeric(gsub('%','',names(quantiles)[i+1])),paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),
                gd,n,'All',paste(most_common_structure,collapse=';')));
        }
    }
    colnames(gd_df)=c('qmin','qmax','repeat_sizeQ','gene_diversity','number','ancestry','most_common_structure')
    # gd_df=gd_df[order(as.numeric(gd_df[,1]),decreasing=F),]
    # gd_df[,1]=as.numeric(gd_df[,1])
    # gd_df$Q=with(gd_df,cut(repeat_size,breaks=quantile(gd$repeat_size,probs=c(0,0.25,0.5,0.75,0.9,0.95,1))))
    return(gd_df)
}

gene_diversity_across_structures=function(data_new,facet_ancestries=F){
    # library(gtools);
    data_new=data_new[order(as.numeric(data_new$repeat_size),decreasing=F),];
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    repeat_sizes=unique(na.omit(data_new$repeat_size));
    # data_new$Q=as.integer(quantcut(as.integer(data_new$repeat_size),q=200))/2;
    gd_df=data.frame();
    counts_df=data.frame();
    percentiles=c(0,0.25,0.50,0.75,0.90,1);
    to=length(percentiles)-1;
    quantiles=quantile(as.numeric(data_new$repeat_size),probs=percentiles);
    if(facet_ancestries){
    for(ancestry in unique(data_new$ancestry)){
        for(structure in unique(data_new$structure)){
            # percentile=quantiles[i];
            # next_percentile=quantiles[i+1];
            d_sub=data_new[(data_new$ancestry==ancestry & as.character(data_new$structure)==structure),];
            # d_sub=data_new[(data_new$ancestry==ancestry & data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            n=nrow(d_sub);
            counts=table(d_sub$repeat_size);
            # counts=counts[repeat_sizes];
            # names(counts)=repeat_sizes;
            # print(counts)
            # counts_df=rbind(counts_df,counts);
            # colnames(counts_df)=repeat_sizes;
            probs=counts/n;
            print(probs)
            gd=(1-sum(probs**2))*(n/(n-1));
            # print(gd)
            # gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),as.numeric(gsub('%','',names(quantiles)[i+1])),
            #     paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),n,gd,ancestry));
            gd_df=rbind(gd_df,cbind(structure,n,gd,ancestry));
        }
    }}else{
        for(structure in unique(data_new$structure)){
            # percentile=quantiles[i];
            # next_percentile=quantiles[i+1];
            d_sub=data_new[data_new$structure==structure,];
            # d_sub=data_new[(data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            n=nrow(d_sub);
            counts=table(d_sub$repeat_size);
            # counts=counts[repeat_sizes];
            # names(counts)=repeat_sizes;
            # counts_df=rbind(counts_df,counts);
            # colnames(counts_df)=repeat_sizes;
            probs=counts/n;
            gd=(1-sum(probs**2))*(n/(n-1));
            # gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),as.numeric(gsub('%','',names(quantiles)[i+1])),
            #     paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),n,gd,'All'));
            gd_df=rbind(gd_df,cbind(structure,n,gd,'All'));
        }
    }
    colnames(gd_df)=c('structure','number','gene_diversity','ancestry');
    # gd_df=gd_df[order(as.numeric(gd_df[,1]),decreasing=F),];
    # gd_df[,1]=as.numeric(gd_df[,1]);
    # gd_df$Q=with(gd_df,cut(repeat_size,breaks=quantile(gd$repeat_size,probs=c(0,0.25,0.5,0.75,0.9,0.95,1))));
    # gd_df=data.frame(cbind(gd_df,counts_df),stringsAsFactors=F);
    return(gd_df);
}



gene_diversity=function(data_new,facet_ancestries=F){
    # library(gtools);
    data_new=data_new[order(as.numeric(data_new$repeat_size),decreasing=F),];
    data_new$repeat_size=as.numeric(data_new$repeat_size)
    structures=unique(na.omit(data_new$structure))
    # data_new$Q=as.integer(quantcut(as.integer(data_new$repeat_size),q=200))/2;
    gd_df=data.frame();
    counts_df=data.frame();
    percentiles=c(0,0.25,0.50,0.75,0.90,1);
    to=length(percentiles)-1;
    quantiles=quantile(as.numeric(data_new$repeat_size),probs=percentiles);
    if(facet_ancestries){
    for(ancestry in unique(data_new$ancestry)){
        # for(uniq_length in unique(data_new$repeat_size)){
            # percentile=quantiles[i];
            # next_percentile=quantiles[i+1];
            d_sub=data_new[(data_new$ancestry==ancestry),];
            # d_sub=data_new[(data_new$ancestry==ancestry & data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            n=nrow(d_sub);
            counts=table(d_sub[,c('repeat_size','structure')]);
            # counts=counts[structures];
            # names(counts)=structures;
            # counts_df=rbind(counts_df,counts);
            # colnames(counts_df)=structures;
            probs=counts/n;
            gd=(1-sum(probs**2))*(n/(n-1));
            # gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),as.numeric(gsub('%','',names(quantiles)[i+1])),
            #     paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),n,gd,ancestry));
            gd_df=rbind(gd_df,cbind(uniq_length,n,probs,gd,ancestry));
        # }
    }}else{
        # for(uniq_length in unique(data_new$repeat_size)){
            # percentile=quantiles[i];
            # next_percentile=quantiles[i+1];
            # d_sub=data_new[data_new$repeat_size==uniq_length,];
            # d_sub=data_new[(data_new$repeat_size>=percentile & data_new$repeat_size<next_percentile),];
            n=nrow(data_new);
            counts=table(data_new[,c('repeat_size','structure')]);
            # counts=counts[structures];
            # names(counts)=structures;
            # counts_df=rbind(counts_df,counts);
            # colnames(counts_df)=structures;
            probs=counts/n;
            gd=(1-sum(probs**2))*(n/(n-1));
            # gd_df=rbind(gd_df,cbind(as.numeric(gsub('%','',names(quantiles)[i])),as.numeric(gsub('%','',names(quantiles)[i+1])),
            #     paste0(gsub('%','',names(quantiles)[i]),'-',names(quantiles)[i+1]),n,gd,'All'));
            gd_df=rbind(gd_df,cbind(uniq_length,n,probs,gd,'All'));
        # }
    }
    colnames(gd_df)=c('repeat_size','number','freq','gene_diversity','ancestry');
    # gd_df=gd_df[order(as.numeric(gd_df[,1]),decreasing=F),]
    # gd_df[,1]=as.numeric(gd_df[,1])
    # gd_df$Q=with(gd_df,cut(repeat_size,breaks=quantile(gd$repeat_size,probs=c(0,0.25,0.5,0.75,0.9,0.95,1))))
    gd_df=data.frame(cbind(gd_df,counts_df),stringsAsFactors=F);
    return(gd_df);
}


gene_diversity_counts_repeat_size=function(data_new,percentiles=NULL,categories=NULL,allowed_gap=20){
    if(!is.null(categories)){
        data_new$category=NA
        data_new$structure=as.character(data_new$structure)
        for(name in names(categories)){
            data_new$category[data_new$structure%in%categories[[name]]]=name
        }
        data_new$category[is.na(data_new$category)]='rare'
        levels=names(categories)
        if('rare'%in%unique(data_new$structure)){
            levels=c(levels,'rare')
        }
        data_new$category=factor(data_new$category,levels=levels)
        data_new$structure=data_new$category
    }
    if(is.null(percentiles)){
        percentiles=seq(0,1,0.01);
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    data_new=data_new[order(data_new$repeat_size,decreasing=F),];
    percentages=data_new%>%group_by(repeat_size,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired")[c(2,1,4,3,6,5,8,7,10,9,12,11)],"#8DD3C7","#9E0142",'#D37731');
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'black')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    g_freq_stacked=ggplot(percentages,aes(x=as.numeric(repeat_size),y=freq,fill=structure))+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
        xlab('Repeat Size Percentile')+
        ylab('Percentages')+fillScale+colScale+theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        theme(title=element_text(size=40),axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position='top',legend.title=element_blank(),
            legend.key.height=unit(1,'cm'),legend.key.width=unit(1,'cm'),legend.text=element_text(face='bold',size=30),strip.text=element_text(size=30))+
            guides(fill=guide_legend(ncol=2))
    return(g_freq_stacked)
}


gene_diversity_counts_handle_percentile_gaps=function(data_new,percentiles=NULL,categories=NULL,allowed_gap=20){
    if(!is.null(categories)){
        data_new$category=NA
        data_new$structure=as.character(data_new$structure)
        for(name in names(categories)){
            data_new$category[data_new$structure%in%categories[[name]]]=name
        }
        data_new$category[is.na(data_new$category)]='rare'
        levels=names(categories)
        if('rare'%in%unique(data_new$structure)){
            levels=c(levels,'rare')
        }
        data_new$category=factor(data_new$category,levels=levels)
        data_new$structure=data_new$category
    }
    if(is.null(percentiles)){
        percentiles=seq(0,1,0.01);
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    data_new=data_new[order(data_new$repeat_size,decreasing=F),];
    quantiles=quantile(data_new$repeat_size,include.lowest=T,right=T,probs=percentiles);
    uniq_names=names(quantiles)[!duplicated(quantiles)];
    quantiles=quantiles[uniq_names];
    # print(quantiles);
    data_new$Q=NA;
    data_new$difference=NA;
    for(r in unique(data_new$repeat_size)){
        quantiles_ordered=quantiles[order(abs(r-quantiles),decreasing=F)];
        closest=as.numeric(gsub('%','',names(quantiles_ordered)[1]))+1;
        second_closest=as.numeric(gsub('%','',names(quantiles_ordered)[2]));

        difference=abs(closest-second_closest);
        if(difference>=allowed_gap){
            data_new$Q[data_new$repeat_size==r]=paste0(min(c(closest,second_closest)),'-',max(c(closest,second_closest)));
        }else{
            data_new$Q[data_new$repeat_size==r]=paste0(closest);
        }
        data_new$difference[data_new$repeat_size==r]='';
    }
    percentages=data_new%>%group_by(Q,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    print(unique(data_new$Q))
    percentages$Q=factor(percentages$Q,levels=unique(data_new$Q))
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731');
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    x_labels=unique(as.character(percentages$Q))
    print(x_labels)
    ranges=x_labels[grep(pattern='-',x_labels)]
    print(ranges)
    x_labels[(x_labels%in%ranges)]=''
    print(x_labels)
    g_freq_stacked=ggplot(percentages,aes(x=as.numeric(Q),y=freq,fill=structure))+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
        xlab('Repeat Size Percentile')+
        ylab('Percentages')+fillScale+colScale+theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        scale_x_discrete(labels=x_labels)+theme(title=element_text(size=40),axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position='top',legend.title=element_blank(),
            legend.key.height=unit(1,'cm'),legend.key.width=unit(1,'cm'),legend.text=element_text(face='bold',size=30),strip.text=element_text(size=30))+
            guides(fill=guide_legend(ncol=2))
    return(g_freq_stacked)
}




gene_diversity_counts=function(data_new,percentiles=NULL,categories=NULL,allowed_gap=20){
    if(!is.null(categories)){
        data_new$category=NA
        data_new$structure=as.character(data_new$structure)
        for(name in names(categories)){
            data_new$category[data_new$structure%in%categories[[name]]]=name
        }
        data_new$category[is.na(data_new$category)]='rare'
        levels=names(categories)
        if('rare'%in%unique(data_new$structure)){
            levels=c(levels,'rare')
        }
        data_new$category=factor(data_new$category,levels=levels)
        data_new$structure=data_new$category
    }
    print(levels(data_new$structure))
    if(is.null(percentiles)){
        percentiles=seq(0,1,0.01);
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    data_new=data_new[order(data_new$repeat_size,decreasing=F),];
    quantiles=quantile(unique(data_new$repeat_size),include.lowest=T,right=T,probs=percentiles);
    uniq_names=names(quantiles)[!duplicated(quantiles)];
    quantiles=quantiles[uniq_names];
    # print(quantiles);
    data_new$Q=NA;
    data_new$difference=NA;
    for(r in unique(data_new$repeat_size)){
        quantiles_ordered=quantiles[order(abs(r-quantiles),decreasing=F)];
        closest=as.numeric(gsub('%','',names(quantiles_ordered)[1]));
        second_closest=as.numeric(gsub('%','',names(quantiles_ordered)[2]));
        difference=abs(closest-second_closest);
        if(difference>=allowed_gap){
            data_new$Q[data_new$repeat_size==r]=paste0(min(c(closest,second_closest)),'-',max(c(closest,second_closest)));
        }else{
            data_new$Q[data_new$repeat_size==r]=paste0(closest);
        }
        data_new$difference[data_new$repeat_size==r]='';
    }
    percentages=data_new%>%group_by(Q,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    percentages$Q=factor(percentages$Q,levels=unique(data_new$Q))
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731');
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'darkgrey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    # x_labels=unique(as.character(percentages$Q))
    # print(x_labels)
    # ranges=x_labels[grep(pattern='-',x_labels)]
    # print(ranges)
    # x_labels[!(x_labels%in%ranges)]=''
    # print(x_labels)
    g_freq_stacked=ggplot(percentages,aes(x=as.numeric(as.character(Q)),y=freq,fill=structure))+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
        xlab('Repeat Size Percentile')+
        ylab('Frequency (%)')+fillScale+colScale+theme_cowplot()+
        theme(title=element_text(size=40),axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position='top',legend.title=element_blank(),
              legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),
              legend.key.height=unit(1,'cm'),legend.key.width=unit(1,'cm'),
              legend.text=element_text(face='bold',size=40),strip.text=element_text(size=30))+
            guides(fill=guide_legend(ncol=2))
    return(g_freq_stacked)
}


gene_diversity_counts_multiple=function(data_new,percentiles=NULL){
    if(is.null(percentiles)){
        percentiles=seq(0,1,0.01);
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    data_new=data_new[order(data_new$repeat_size,decreasing=F),];
    quantiles=quantile(data_new$repeat_size,include.lowest=T,right=T,probs=percentiles);
    uniq_names=names(quantiles)[!duplicated(quantiles)];
    quantiles=quantiles[uniq_names];
    print(quantiles);
    data_new$Q=NA;
    data_new$difference=NA
    for(r in unique(data_new$repeat_size)){
        # print(r);
        quantiles_ordered=quantiles[order(abs(r-quantiles),decreasing=F)];
        closest=as.numeric(gsub('%','',names(quantiles_ordered)[1]));
        second_closest=as.numeric(gsub('%','',names(quantiles_ordered)[2]));
        difference=abs(closest-second_closest);

        if(difference>=1){
            data_new$Q[data_new$repeat_size==r]=closest
            tmp=data_new[data_new$repeat_size==r,]
            print('gap')
            print(closest)
            print(second_closest)
            gap_filler=seq(min(c(closest,second_closest)),max(c(closest,second_closest)),1)
            gap_filler=gap_filler[2:length(gap_filler)]
            for(d in gap_filler){
                print(d)
                tmp$Q=d
                data_new=rbind(data_new,tmp)
            }
        }else{
            data_new$Q[data_new$repeat_size==r]=closest
        }
        data_new$difference[data_new$repeat_size==r]=difference
    }
    # data_new$Q=cut(data_new$repeat_size,breaks=quantiles,include.lowest=T,right=T,labels=uniq_names[2:length(uniq_names)]);
    percentages=data_new%>%group_by(Q,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    percentages$Q=factor(percentages$Q,levels=unique(data_new$Q))
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731');
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    g_freq_stacked=ggplot(percentages,aes(x=as.numeric(Q),y=freq,fill=structure))+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        xlab('Repeat Quantile')+
        ylab('Percentages')+fillScale+colScale+theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position='none')
    return(g_freq_stacked)
}




gene_diversity_counts_new=function(data_new,percentiles=NULL){
    if(is.null(percentiles)){
        percentiles=seq(0,1,0.01);
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    data_new=data_new[order(data_new$repeat_size,decreasing=F),];
   
    # data_new$Q=cut(data_new$repeat_size,breaks=quantiles,include.lowest=T,right=T,labels=uniq_names[2:length(uniq_names)]);
    percentages=data_new%>%group_by(repeat_size,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    # percentages$Q=factor(percentages$Q,levels=unique(data_new$Q))
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731');
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    g_freq_stacked=ggplot(percentages,aes(x=repeat_size,y=freq,fill=structure))+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        xlab('Repeat Quantile')+
        ylab('Percentages')+fillScale+colScale+theme_cowplot()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position='none')
    return(g_freq_stacked)
}




gene_diversity_counts2=function(data_new,percentiles=NULL){
    if(is.null(percentiles)){
        # percentiles=seq(0,1,0.01)
        percentiles=c(0,.75,seq(0.80,1,0.05))
        continuous=F
    }else if(percentiles=='0-100'){
        percentiles=seq(0,1,0.01)
        continuous=T
    }else{
        continuous=F
    }
    data_new$repeat_size=as.numeric(data_new$repeat_size);
    quantiles=quantile(data_new$repeat_size,probs=percentiles);
    uniq_names=names(quantiles)[!duplicated(quantiles)];
    quantiles=quantiles[uniq_names];
    labels=c()
    to=length(labels)-1
    for(i in 1:to){
        labels=c(labels,paste0(gsub('%','',uniq_names[i]),'-',uniq_names[i+1]))
    }
    data_new$Q=cut(data_new$repeat_size,breaks=quantiles,include.lowest=T,right=T,labels=labels);
    percentages=data_new%>%group_by(Q,ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100)
    print(percentages)
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731')
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    allele_order=levels(data_new$structure)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    # as.numeric(gsub('%','',as.character(Q))
    g_freq_stacked=ggplot(percentages,aes(x=Q),y=freq,fill=structure)+
        geom_bar(stat='identity')+facet_wrap(~ancestry,nrow=1)+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        xlab('Repeat Quantile')+
        ylab('Percentages')+fillScale+colScale+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20))
    return(g_freq_stacked)
}

melt_columns_and_plot_structure_across_ancestries_specialised=function(data,struct_cols,rc_cols,
                                        eh_cols,get_rid=NULL,allele_n_threshold=0,
                                        name,do_comparisons=NULL,file_suffix='_realigned.bam',
                                        percentage_threshold=0,allele_order=NULL,
                                        pathogenic_repeats=NULL,facet_ancestries=FALSE,title=NULL,get_rid_of_numbers=TRUE){
    out=data.frame()
    data=bind_ancestries(data,file_suffix)
    print(head(data))
    for(i in 1:length(struct_cols)){
        if(!is.null(rc_cols)){
            d=cbind(data[,struct_cols[i]],data[,rc_cols[i]],data[,eh_cols[i]],data$superpopu)
        }else{
            d=cbind(data[,struct_cols[i]],data[,eh_cols[i]],data$superpopu)
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    # print(head(out))
    if(!is.null(rc_cols)){
        colnames(out)=c('structure','rc_rep','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
    }else{
        print(head(out))
        colnames(out)=c('structure','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    out=na.omit(out)
    percentage=data.frame(out,stringsAsFactors=F) %>%
        group_by(ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100,
            upperE=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100,
            lowerE=(-1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100)
    percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    if(get_rid_of_numbers){percentage$structure=gsub("(?<!x)\\d","",percentage$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        percentage$structure[((percentage$freq)<percentage_threshold) & !(percentage$structure%in%pathogenic_repeats)]='rare'
        percentage$structure[(percentage$structure%in%pathogenic_repeats)]=paste('Pathogenic:',percentage$structure[(percentage$structure%in%pathogenic_repeats)])
    }else{
        a=data.frame(percentage[(percentage$freq)>=percentage_threshold,])
        print(a[!(a$structure%in%allele_order),])
        percentage$structure[(percentage$freq)<percentage_threshold]='rare'
    }
    # percentage=na.omit(percentage)
    out$structure_annotated=NA
    out$ancestry_annotated=NA
    out_new=data.frame()
    for(t in unique(out$structure)){
        out_sub=out[out$structure==t,]
        # if(length(unique(out_sub$ancestry))>1){
            for(ancestry in unique(out_sub$ancestry)){
                out_sub2=out_sub[out_sub$ancestry==ancestry,]
                # if(ancestry!='rare'){
                    n=nrow(out_sub2)
                # }else{
                #     n=nrow(out[out$ancestry==ancestry,])
                # }
                if(n>0){
                    out_sub2=out_sub2[out_sub2$structure==t,]
                    # if(n<allele_n_threshold){
                    #     out_sub2=out_sub2[-which(out_sub2$structure==t),]
                    # }else{
                        out_sub2$structure_annotated[out_sub2$structure==t]=paste0(t,' N=',n)
                        out_sub2$ancestry_annotated=paste0(ancestry,' N=',n)
                    # }
                    out_new=rbind(out_new,out_sub2)
                }
            }
        # }
    }
    out=out_new
    if(get_rid_of_numbers){out$structure=gsub("(?<!x)\\d","",out$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        pathogenic_sub=out[out$structure%in%pathogenic_repeats,]
        pathogenic_sub$structure=paste('Pathogenic:',pathogenic_sub$structure)
        pathogenic_sub$structure_annotated=paste('Pathogenic:',pathogenic_sub$structure_annotated)
        out=out[!(out$structure%in%pathogenic_repeats),]
        # print(nrow(out))
    }
    if(!is.null(percentage_threshold)){
        percentage2=percentage[percentage$freq>=percentage_threshold,]
        out_new=data.frame()
        for(ancestry in unique(percentage2$ancestry)){
            percentage_sub=percentage2[percentage2$ancestry==ancestry,]
            out_sub=out[out$ancestry==ancestry,]
            out_sub$structure[!(out_sub$structure%in%unique(percentage_sub$structure))]='rare'
            out_sub$ancestry_annotated[out_sub$structure=='rare']=paste0(ancestry,' N=',sum(out_sub$structure=='rare'))
            out_new=rbind(out_new,out_sub)
        }
        out=out_new
    }
    if(!is.null(pathogenic_repeats)){
        out=rbind(out,pathogenic_sub)
        print(nrow(out))
    }
    if(!is.null(allele_order)){
        if(get_rid_of_numbers){allele_order=gsub("(?<!x)\\d","",allele_order,perl=T)}
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }else{
        allele_order=unique(out$structure)
        # print('s')
        print(allele_order)
        if(('rare'%in%allele_order)){
            allele_order=allele_order[allele_order!='rare']
            allele_order=c(allele_order,'rare')
        }
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731')
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    if(facet_ancestries){
        g=ggplot(na.omit(out),aes(y=ancestry_annotated,x=as.numeric(as.character(repeat_size)),fill=structure))+
            geom_violin()+
            geom_boxplot()+
            theme_cowplot()+
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))+
            facet_wrap(~structure,ncol=1,scales='free_y')
    }else{
        g=ggplot(na.omit(out),aes(y=factor(structure,levels=rev(allele_order)),x=as.numeric(as.character(repeat_size)),fill=structure))+
            geom_violin()+
            geom_boxplot()+
            theme_cowplot()+
            theme(legend.key.size = unit(2, "cm"),axis.title.y=element_text(size=20),axis.title.x=element_text(size=25),axis.text.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))
        for(a in allele_order){
            print(a)
            nalleles=nrow(out[as.character(out$structure)==a,])
            print(nalleles)
            if(nalleles>0){
                g=g+annotate(geom='text',y=a, x=round(max(as.numeric(as.character(out$repeat_size)))*0.9,0),
                            size=8, label= paste0('N=',nalleles),vjust=1.6,hjust='inwards')
            }
            # print(nalleles)
        }
    }
    g=g+fillScale
    if(!is.null(do_comparisons)){
        comparisons=na.omit(unique(out$structure))
        comparisons=t(combn(comparisons,2))
        combinations=list()
        for(i in 1:nrow(comparisons)){
            combinations[[i]]=c(comparisons[i,])
        }
        # print(combinations)
        g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1)
    }
    # g=g+theme(legend.key.height=unit(0.2,'cm'),legend.key.width=unit(0.2,'cm'),legend.height=element_text(size=5),legend.text=element_text(size=5))
    # print(percentage[is.na(percentage$ancestry),])
    # print(percentage[is.na(percentage$structure),])
    # percentage$structure[is.na(percentage$structure)]='rare'
    percentage$label=paste0(round(percentage$freq,1),'%')
    percentage$label2=paste0(round(percentage$freq,2),'%')
    percentage$label[percentage$freq<4]=NA
    percentage$label2[!is.na(percentage$label)]=NA
    g_freq_stacked=ggplot(percentage,aes(y=ancestry,x=freq,col=structure,fill=structure,label=label))+
        geom_bar(stat='identity')+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black") +
        ylab('Population')+
        xlab('Percentages')+theme_cowplot()+geom_text(color='white', position = position_stack(vjust = 0.5))+
            ggrepel::geom_text_repel(aes(label=label2),position = position_stack(vjust = 0.5))+
            theme(legend.position='none',axis.title.y=element_text(size=10),axis.title.x=element_text(size=25),
                axis.text.y = element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20))
    g_freq_stacked=g_freq_stacked+fillScale+colScale
    # percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    # g_freq=ggplot(na.omit(percentage),aes(y=ancestry_annotated,x=freq))+
    #     geom_bar(stat='identity',col='black',fill='white')+
    #     geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
    #     ylab('Repeat structure')+
    #     xlab('Percentages')+
    #     facet_wrap(~factor(structure),ncol=1,scales='free')+theme_cowplot()
    
    g=ggarrange(g_freq_stacked+ggtitle(name)+theme(plot.title=element_text(size=20)),g,nrow=2,ncol=1,heights=c(1,3))
    # if(!is.null(title)){
    #     g=g+ggtitle(title)
    # }
    
    #f=paste0('tmp2/structure_vs_eh_repsize_ancestries_',name,'.pdf')
    #file.remove(f)
    #pdf(f,height=10,width=15)
    #print(g)
    #dev.off()
    #print(f)
    return(list(sizes=out,proportions=percentage,plot=g,plot_proportions=g_freq_stacked))
}



chisq_ancestries=function(data,out_name,thresh=0,categories=NULL){
    population_meta=fread('/re_gecip/neurology/Kristina/POPULATION/list_34191_unrel_notNeuro_withSuperpopu.tsv',sep='\t',header=T)
    # print(head(data))
    data$bam_file=basename(gsub(paste0('EH_'),"",gsub(paste0('_',out_name),"",gsub(paste0('_realigned.bam'),"",data$bam_file))))
    # print(head(data))
    data=merge(data,population_meta,by.x=c('bam_file'),by.y=c('platekey'))
    # print(head(data))
    counts_pop=data.frame()
    for(p in unique(data$superpopu)){
        data_sub=data[data$superpopu==p,]
        data_sub$first_second_gt=paste(data_sub$GT1,data_sub$GT2,sep='/')
        tab=table(as.character(data_sub$first_second_gt))
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
        total=nrow(data_sub)*2
        counts$freq=counts$count/total
        counts=counts[counts$freq>thresh,]
        counts=counts[order(counts$count,decreasing=T),]
        counts$population=p
        counts_pop=rbind(counts_pop,counts)
    }
    
    counts_pop$count=as.numeric(counts_pop$count)
    data_new=counts_pop
    # print(head(data_new))
    if(!is.null(categories)){
        data_new$category=NA
        data_new$GT=as.character(data_new$GT)
        print(unique(data_new$GT))
        for(name in names(categories)){
            print(name)
            data_new$category[data_new$GT%in%categories[[name]]]=name
        }
        data_new=data_new[!is.na(data_new$category),]
        # levels=names(categories)
        # if('rare'%in%unique(data_new$GT)){
        #     levels=c(levels,'rare')
        # }
        # data_new$category=factor(data_new$category,levels=levels)
        data_new$GT=data_new$category
        data_new$category=as.character(data_new$category)
        # data_new=na.omit(data_new)
    }
    cols=unique(data_new$population)
    rows=na.omit(unique(data_new$GT))
    rows=rows[rows!='NA']
    tab=matrix(ncol=length(cols),nrow=length(rows))
    colnames(tab)=cols
    rownames(tab)=rows
    tab_compare=matrix(ncol=length(cols),nrow=length(rows))
    colnames(tab_compare)=cols
    rownames(tab_compare)=rows
    # print(tab_compare)
    for(p in cols){
        data_sub=data_new[data_new$population==p,]
        # print(data_sub)
        n=data_sub$count
        ci=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)
        data_sub$freq=round((data_sub$freq)*100,2)
        ci=round((ci)*100,2)
        data_sub$freq=data_sub$count
        # data_sub$freq=paste0(data_sub$freq,'%\nN=',data_sub$count)
        for(cat in rows){
            if(nrow(data_sub[data_sub$GT==cat,])>0){
                tab[cat,p]=data_sub$freq[data_sub$GT==cat]
                tab_compare[cat,p]=data_sub$count[data_sub$GT==cat]
            }
        }
    }
    tab_compare[is.na(tab_compare)]=0
    tab_compare=tab_compare[order(rownames(tab_compare)),]
    tab_compare=tab_compare[!is.na(rownames(tab_compare)),]
    print(tab_compare)
    test=chisq.test(tab_compare)
    pval=test$p.value
    if(grepl(pattern='e',as.character(pval))){
        split=unlist(strsplit(as.character(pval),'e'))
        pval=paste0(round(as.numeric(split[1]),1),split[2])
    }else{
        pval=as.character(round(as.numeric(pval),1))
    }
    residuals=test$residuals    
    plot_data=data.frame()
    for(i in colnames(residuals)){
        for(j in rownames(residuals))
            plot_data=rbind(plot_data,cbind(j,i,residuals[j,i]))
    }
    library(viridis)
    colnames(plot_data)=c('structure','ancestry','residual')
    plot_data$residual=as.numeric(as.character(plot_data$residual))
    g=ggplot(plot_data,aes(x=ancestry,y=structure,size=residual,col=residual))+geom_point()+scale_color_viridis_c()+
        ggtitle(paste0(out_name,' X^2=',round(as.numeric(test$statistic),0),' P-value=',pval))+
        theme_minimal()+ggplot2::scale_size(range = c(4, 15)) +guides(size=F)+theme(title=element_text(size=40),
              legend.position='top',axis.title=element_text(size=20),axis.text=element_text(size=20),
              legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,0,0),
              legend.key.height=unit(1,'cm'),legend.key.width=unit(1,'cm'),
              legend.text=element_text(face='bold',size=20))+xlab('')+ylab('')
    pdf(paste0('chisq/',out_name,'_chisqed.pdf'),height=10,width=26)
    # corrplot(test$residuals,is.corr=F,col=colorRampPalette(c("blue","white","red"))(9))
    print(g)
    dev.off()
    return(g)
}




melt_columns_and_plot_structure_across_ancestries_and_undone=function(data,struct_cols,rc_cols,
                                        eh_cols,get_rid=NULL,allele_n_threshold=0,
                                        name,do_comparisons=NULL,file_suffix='_realigned.bam',
                                        percentage_threshold=0,allele_order=NULL,
                                        pathogenic_repeats=NULL,facet_ancestries=FALSE,
                                        title=NULL,get_rid_of_numbers=TRUE,phase_to=NULL,
                                        full_list_of_platekeys=NULL){
    if(!is.null(full_list_of_platekeys)){
        data$bam_file=basename(data$bam_file)
        data$bam_file=gsub('EH_','',data$bam_file)
        data$bam_file=gsub(file_suffix,'',data$bam_file)
        platekeys_not_annotated=full_list_of_platekeys[!(full_list_of_platekeys%in%data$bam_file)]
        print('not annotated')
        for(p in platekeys_not_annotated){
            data=rbind(data,c(p,rep(NA,ncol(data)-1)))
        }
    }
    print(tail(data))
    if(!is.null(phase_to)){
        data=switcharoo(data,phase_to)
    }
    out=data.frame()
    data=bind_ancestries(data,file_suffix)
    print(tail(data[data$bam_file%in%platekeys_not_annotated,]))
    for(i in 1:length(struct_cols)){
        if(!is.null(rc_cols)){
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,rc_cols[i]],data[,eh_cols[i]],data$superpopu)
        }else{
            d=cbind(data[,'bam_file'],data[,struct_cols[i]],data[,eh_cols[i]],data$superpopu)
        }
        # d[d[,3]==d[,4],4]=NA
        out=rbind(out,d)
    }
    if(!is.null(rc_cols)){
        colnames(out)=c('bam_file','structure','rc_rep','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
        # print('hereeee')
        # print(dim(out))
        # print(head(out[as.integer(out$rc_rep)!=as.integer(out$repeat_size),]))
        out=out[as.integer(out$rc_rep)==as.integer(out$repeat_size),]
        print(dim(out))
    }else{
        # print(nrow(out))
        colnames(out)=c('bam_file','structure','repeat_size','ancestry')
        out=data.frame(out,stringsAsFactors=F)
    }
    if(!is.null(get_rid)){
        for(s in get_rid){
            print(s)
            out$structure=gsub(paste0(s,'\\|'),'',out$structure)
        }
    }
    out=unique(out)
    out=out[!is.na(out$bam_file),]
    out_extra_exploration=out
    print(head(out_extra_exploration[is.na(out_extra_exploration[,2]),]))
    percentage=data.frame(out,stringsAsFactors=F) %>%
        group_by(ancestry,structure) %>%
        summarise(n = n()) %>%
        mutate(total=sum(n),freq = (n / sum(n))*100,
            upperE=(1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100,
            lowerE=(-1.96*sqrt(n/sum(n))*(1-(n/sum(n)))/n)*100)
    percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    if(get_rid_of_numbers){percentage$structure=gsub("(?<!x)\\d","",percentage$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        percentage$structure[((percentage$freq)<percentage_threshold) & !(percentage$structure%in%pathogenic_repeats)]='rare'
        percentage$structure[(percentage$structure%in%pathogenic_repeats)]=paste('Pathogenic:',percentage$structure[(percentage$structure%in%pathogenic_repeats)])
    }else{
        a=data.frame(percentage[(percentage$freq)>=percentage_threshold,])
        print(a[!(a$structure%in%allele_order),])
        percentage$structure[(percentage$freq)<percentage_threshold]='rare'
    }
    # percentage=na.omit(percentage)
    out=na.omit(out)
    out$structure_annotated=NA
    out$ancestry_annotated=NA
    out_new=data.frame()
    for(t in unique(out$structure)){
        out_sub=out[out$structure==t,]
        # if(length(unique(out_sub$ancestry))>1){
            for(ancestry in unique(out_sub$ancestry)){
                out_sub2=out_sub[out_sub$ancestry==ancestry,]
                # if(ancestry!='rare'){
                    n=nrow(out_sub2)
                # }else{
                #     n=nrow(out[out$ancestry==ancestry,])
                # }
                if(n>0){
                    out_sub2=out_sub2[out_sub2$structure==t,]
                    if(n<allele_n_threshold){
                        out_sub2=out_sub2[-which(out_sub2$structure==t),]
                    }else{
                        out_sub2=unique(out_sub2)
                        out_sub2$structure_annotated[out_sub2$structure==t]=paste0(t,' N=',n)
                        out_sub2$ancestry_annotated=paste0(ancestry,' N=',n)
                    }
                    out_new=rbind(out_new,out_sub2)
                }
            }
        # }
    }
    out=out_new
    out=unique(out)

    if(get_rid_of_numbers){out$structure=gsub("(?<!x)\\d","",out$structure,perl=T)}
    if(!is.null(pathogenic_repeats)){
        pathogenic_sub=out[out$structure%in%pathogenic_repeats,]
        pathogenic_sub$structure=paste('Pathogenic:',pathogenic_sub$structure)
        pathogenic_sub$structure_annotated=paste('Pathogenic:',pathogenic_sub$structure_annotated)
        out=out[!(out$structure%in%pathogenic_repeats),]
        # print(nrow(out))
    }
    if(!is.null(percentage_threshold)){
        percentage2=percentage[percentage$freq>=percentage_threshold,]
        out_new=data.frame()
        for(ancestry in unique(percentage2$ancestry)){
            percentage_sub=percentage2[percentage2$ancestry==ancestry,]
            out_sub=out[out$ancestry==ancestry,]
            out_sub$structure[!(out_sub$structure%in%unique(percentage_sub$structure))]='rare'
            out_sub=unique(out_sub)
            out_sub$ancestry_annotated[out_sub$structure=='rare']=paste0(ancestry,' N=',sum(out_sub$structure=='rare'))
            out_new=rbind(out_new,out_sub)
        }
        out=out_new
    }
    out=unique(out)
    if(!is.null(pathogenic_repeats)){
        out=rbind(out,pathogenic_sub)
        print(nrow(out))
    }
    if(!is.null(allele_order)){
        if(get_rid_of_numbers){allele_order=gsub("(?<!x)\\d","",allele_order,perl=T)}
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }else{
        allele_order=unique(out$structure)
        # print('s')
        print(allele_order)
        if(('rare'%in%allele_order)){
            allele_order=allele_order[allele_order!='rare']
            allele_order=c(allele_order,'rare')
        }
        out$structure=factor(out$structure,levels=allele_order)
        # print(percentage)
        percentage$structure=factor(percentage$structure,levels=rev(allele_order))
        # print(percentage)
    }
    myColors <- c(RColorBrewer::brewer.pal(12, "Paired"),"#8DD3C7","#9E0142",'#D37731')
    # colScale <- scale_colour_manual(name = "ancestry",values = myColors)
    if('rare' %in% allele_order){
        to=length(allele_order)-1
        myColors=c(myColors[1:to],'grey')
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = myColors)
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = myColors)
    }else{
        myColors=c(myColors[1:length(allele_order)])
        fillScale <- scale_fill_manual(name = "structure",breaks=allele_order,values = c(myColors))
        colScale <- scale_color_manual(name = "structure",breaks=allele_order,values = c(myColors))
    }
    if(facet_ancestries){
        g=ggplot(na.omit(out),aes(y=ancestry_annotated,x=as.numeric(as.character(repeat_size)),fill=structure,group=ancestry_annotated))+
            geom_violin()+
            geom_boxplot(width=0.1)+
            theme_cowplot()+
            theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))+
            facet_wrap(~structure,ncol=1,scales='free_y')
    }else{
        g=ggplot(na.omit(out),aes(y=factor(structure,levels=rev(allele_order)),x=as.numeric(as.character(repeat_size)),fill=structure))+
            geom_violin()+
            geom_boxplot(width=0.1)+
            theme_cowplot()+
            theme(legend.key.size=unit(2, "cm"),axis.title.y=element_text(size=20),axis.title.x=element_text(size=25),axis.text.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
            ylab('')+
            xlab(paste('Expansion Hunter Repeat size'))
        for(a in allele_order){
            print(a)
            nalleles=nrow(out[as.character(out$structure)==a,])
            print(nalleles)
            if(nalleles>0){
                g=g+annotate(geom='text',y=a, x=round(max(as.numeric(as.character(out$repeat_size)))*0.9,0),
                            size=8, label= paste0('N=',nalleles),vjust=1.6,hjust='inwards')
            }
            # print(nalleles)
        }
    }
    g=g+fillScale
    if(!is.null(do_comparisons)){
        print('compare')
        if(!facet_ancestries){
            comparisons=na.omit(unique(as.character(out$structure)))
            comparisons=t(combn(comparisons,2))
            combinations=list()
            for(i in 1:nrow(comparisons)){
                combinations[[i]]=c(comparisons[i,])
            }
        }else{
            comparisons=na.omit(unique(as.character(out$ancestry_annotated)))
            comparisons=t(combn(comparisons,2))
            combinations=list()
            for(i in 1:nrow(comparisons)){
                combinations[[i]]=c(comparisons[i,])
            }
        }
        print(combinations)

        g=g+geom_signif(comparisons=combinations,map_signif_level=TRUE,step_increase=0.1,orientation='y')
        # g=g+stat_compare_means(method = "wilcox.test", ref.group='.all.')
    }
    g_freq_stacked=ggplot(percentage,aes(y=ancestry,x=freq,fill=structure))+
        geom_bar(stat='identity')+
        # geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
        ylab('Population')+
        xlab('Percentages')+theme_cowplot()+
            theme(legend.position='none',axis.title.y=element_text(size=10),axis.title.x=element_text(size=25),axis.text.y = element_text(size=20),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=20))
    g_freq_stacked=g_freq_stacked+fillScale
    # percentage$ancestry_annotated=paste0(percentage$ancestry,' N=',percentage$n)
    # g_freq=ggplot(na.omit(percentage),aes(y=ancestry_annotated,x=freq))+
    #     geom_bar(stat='identity',col='black',fill='white')+
    #     geom_errorbar(aes(xmin = freq+lowerE, xmax = freq+upperE), width = .2, col = "black")+
    #     ylab('Repeat structure')+
    #     xlab('Percentages')+
    #     facet_wrap(~factor(structure),ncol=1,scales='free')+theme_cowplot()
    
    g1=ggarrange(g_freq_stacked+ggtitle(name)+theme(plot.title=element_text(size=20)),g,nrow=2,ncol=1,heights=c(1,3))
    # if(!is.null(title)){
    #     g=g+ggtitle(title)
    # }
    
    f=paste0('structure_vs_eh_repsize_ancestries_',name,'.pdf')
    file.remove(f)
    pdf(f,height=20,width=15)
    print(g1)
    dev.off()
    print(f)
    return(list(sizes=out,extra_exploration=out_extra_exploration,proportions=percentage,plot=g,plot_proportions=g_freq_stacked,platekeys_not_annotated=platekeys_not_annotated))
}



