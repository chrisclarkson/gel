library(data.table)
library(ggplot2)
library(ggpubr)


sharp_plot=function(repeat_size_data,platekey_col,gene,gene_col,
                    LongAllele_col,breaks,case_solved_family_data,vcf_name,title,
                    control_data=NULL,control_platekey_col=NULL,eh_classifier_col=NULL){
    #read and organise repeat size table
    repeat_size_data=data.frame(fread(repeat_size_data),stringsAsFactors=F)
    repeat_size_data=repeat_size_data[repeat_size_data[,gene_col]==gene,]
    repeat_size_data=repeat_size_data[!is.na(repeat_size_data[,LongAllele_col]),]
    if(!is.null(eh_classifier_col)){
        repeat_size_data=repeat_size_data[repeat_size_data[,eh_classifier_col]=='T',]
    }
    repeat_size_data[,platekey_col]=gsub(paste0(vcf_name,'.vcf'),'',repeat_size_data[,platekey_col])
    case_solved_family_data=fread(case_solved_family_data,sep='\t',header=T)
    case_solved_family_data$case_solved_family[case_solved_family_data$case_solved_family!='yes']='no'
    # cases solved column is supplied via the `case_solved_family_data` table
    repeat_size_data=merge(repeat_size_data,case_solved_family_data,by.x=c(platekey_col),by.y=c('plate_key'),all=T)
    unsolved_cases=na.omit(repeat_size_data[repeat_size_data$case_solved_family=='no',])
    unsolved_cases=unsolved_cases[!is.na(unsolved_cases$case_solved_family),]
    solved_cases=na.omit(repeat_size_data[repeat_size_data$case_solved_family=='yes',])
    solved_cases=solved_cases[!is.na(solved_cases$case_solved_family),]
    #if control list is supplied- use it to separate controls in repeat data
    if(is.null(control_data)){
        controls=repeat_size_data[is.na(repeat_size_data$case_solved_family),]
    }else{
        if(is.null(control_platekey_col)){control_platekey_col=platekey_col}
        controls=data.frame(fread(control_data),stringsAsFactors=F)[,control_platekey_col]
        controls=repeat_size_data[repeat_size_data[,platekey_col]%in%controls,]
    }
    print(paste(nrow(repeat_size_data[!is.na(repeat_size_data$case_solved_family),]),'cases'))
    print(paste(nrow(unsolved_cases),'unsolved cases'))
    print(paste(nrow(solved_cases),'solved cases'))
    print(paste(nrow(controls),'controls'))
    #create histogram dataframe and bins in which odds ratio will be calculated
    h=hist(repeat_size_data[,LongAllele_col],breaks=breaks,plot=F)
    breaks=h$breaks
    breaks=breaks[-length(breaks)]
    counts=h$counts
    hist_data=data.frame(breaks,counts,stringsAsFactors=F)
    as_far_as=length(breaks)
    # construct odds ratio table
    odds_ratio_data=data.frame()
    for(i in 1:as_far_as){
        unsolved_cases_gt=na.omit(unsolved_cases[(unsolved_cases[,LongAllele_col]>=breaks[i]),])
        unsolved_cases_lt=na.omit(unsolved_cases[(unsolved_cases[,LongAllele_col]<breaks[i]),])
        solved_cases_gt=na.omit(solved_cases[(solved_cases[,LongAllele_col]>=breaks[i]),])
        solved_cases_lt=na.omit(solved_cases[(solved_cases[,LongAllele_col]<breaks[i]),])
        controls_gt=controls[controls[,LongAllele_col]>=breaks[i],]
        controls_lt=controls[controls[,LongAllele_col]<breaks[i],]
        # odds_ratio=(sum(d_sub$case_solved_family=='no')/sum(d_sub$case_solved_family=='yes'))
        ratio1=(nrow(unsolved_cases_gt)+1)/(nrow(unsolved_cases_lt)+1)
        ratio2=(nrow(controls_gt)+1)/(nrow(controls_lt)+1)
        odds_ratio=ratio1/ratio2
        test1=fisher.test(matrix(
            c(  nrow(unsolved_cases_gt)+1,nrow(unsolved_cases_lt)+1,
                nrow(controls_gt)+1,nrow(controls_lt)+1),
                byrow=T,nrow=2))
        test2=fisher.test(matrix(
            c(  nrow(solved_cases_gt)+1,nrow(solved_cases_lt)+1,
                nrow(controls_gt)+1,nrow(controls_lt)+1),
                byrow=T,nrow=2))
        pval=test1$p.value
        df=rbind(cbind(breaks[i],test1$estimate,'Odds Ratio Unsolved Probands'),
            cbind(breaks[i],test2$estimate,'Odds Ratio Solved Probands'),
            cbind(breaks[i],pval,'-log10 P-value'))
        odds_ratio_data=rbind(odds_ratio_data,df)
    }
    colnames(odds_ratio_data)=c('breaks','number','metric');
    odds_ratio_data=data.frame(odds_ratio_data,stringsAsFactors=F);
    write.table(odds_ratio_data,paste0(gene,'_OR_table.tsv'),sep='\t',col.names=F,row.names=F,quote=F)

    #plot data
    g_hist=ggplot(hist_data,aes(x=breaks,y=counts))+geom_col()+scale_y_log10()+ggtitle(title)+xlab('')+ylab('Count')
    # g_hist=ggplot(d,aes(x=A2,col=case_solved_family))+
    #     geom_histogram()+scale_y_log10()
    g1=ggplot(odds_ratio_data,aes(x=as.numeric(breaks),y=as.numeric(number),col=metric))+
        geom_point()+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position=c(.25,.75))+
        xlab('Repeat sizes')+
        ylab('Odds Ratio/ -log10 P-value')+scale_y_log10()
    g=ggpubr::ggarrange(g_hist,g1,ncol = 1, nrow = 2)
    return(g)
}



sharp_plot_delia=function(repeat_size_data,platekey_col,gene,gene_col,
                    LongAllele_col,breaks,case_solved_family_data,vcf_name,title,
                    control_data=NULL,control_platekey_col=NULL,
                    legend_position=c(.75,.75),include_ci=TRUE,handle_infinite_values='set to zero',start_from=0,export_used_data=FALSE,eh_classifier_col=NULL){
    # read and organise repeat size table
    repeat_size_data=data.frame(fread(repeat_size_data),stringsAsFactors=F)
    repeat_size_data=repeat_size_data[repeat_size_data[,gene_col]==gene,]
    repeat_size_data=repeat_size_data[!is.na(repeat_size_data[,LongAllele_col]),]
    print(sum(is.na(repeat_size_data[,LongAllele_col])))
    repeat_size_data=repeat_size_data[as.numeric(repeat_size_data[,LongAllele_col])>=start_from,]
    if(!is.null(eh_classifier_col)){
        repeat_size_data=repeat_size_data[repeat_size_data[,eh_classifier_col]=='T',]
    }
    repeat_size_data[,platekey_col]=basename(gsub(paste0(vcf_name,'.vcf'),'',repeat_size_data[,platekey_col]))
    case_solved_family_data=fread(case_solved_family_data,sep='\t',header=T,fill=T)
    case_solved_family_data$case_solved_family[is.na(case_solved_family_data$case_solved_family)]='control'
    case_solved_family_data$case_solved_family[which(case_solved_family_data$case_solved_family!='yes')]='no'
    # cases solved column is supplied via the `case_solved_family_data` table
    repeat_size_data=merge(repeat_size_data,case_solved_family_data,by.x=c(platekey_col),by.y=c('plate_key'),all=T)
    
    if('case_or_control'%in%colnames(repeat_size_data)){
        print('hisdf')
        print(head(repeat_size_data))
        print(table(repeat_size_data$case_solved_family))
            unsolved_cases=repeat_size_data[which(repeat_size_data$case_solved_family=='no' & repeat_size_data$case_or_control=='case') ,]
            unsolved_cases=unsolved_cases[!is.na(unsolved_cases$case_solved_family),]
            solved_cases=repeat_size_data[which(repeat_size_data$case_solved_family=='yes' & repeat_size_data$case_or_control=='case'),]
            solved_cases=solved_cases[!is.na(solved_cases$case_solved_family),]
            
        }else{
            print('in here')
            unsolved_cases=repeat_size_data[which(repeat_size_data$case_solved_family=='no'),]
            unsolved_cases=unsolved_cases[!is.na(unsolved_cases$case_solved_family),]
            solved_cases=repeat_size_data[which(repeat_size_data$case_solved_family=='yes'),]
            solved_cases=solved_cases[!is.na(solved_cases$case_solved_family),]
        }
    # if control list is supplied- use it to separate controls in repeat data
    if(is.null(control_data)){
        if('case_or_control'%in%colnames(repeat_size_data)){
            controls=repeat_size_data[which(repeat_size_data$case_or_control=='control'),]
        }else{
            controls=repeat_size_data[is.na(repeat_size_data$case_solved_family),]
            print('dddd')
            controls=controls[!(controls[,platekey_col]%in%c(unsolved_cases[,platekey_col],solved_cases[,platekey_col])),]
        }
    }else{
        if(is.null(control_platekey_col)){control_platekey_col=platekey_col}
        controls=data.frame(fread(control_data),stringsAsFactors=F)[,control_platekey_col]
        controls=repeat_size_data[repeat_size_data[,platekey_col]%in%controls,]
    }
    print(paste(nrow(repeat_size_data[!is.na(repeat_size_data$case_solved_family),]),'cases'))
    print(paste(nrow(unsolved_cases),'unsolved cases'))
    print(paste(nrow(solved_cases),'solved cases'))
    print(paste(nrow(controls),'controls'))
    if(export_used_data){
        unsolved_cases$category='unsolved cases'
        solved_cases$category='solved cases'
        controls$category='controls'
        export=rbind(unsolved_cases,solved_cases,controls)
        write.table(export[order(export[,LongAllele_col],decreasing=T),],paste0(gene,'_used_data.tsv'),sep='\t',row.names=F,col.names=F,quote=F)
    }
    # create histogram dataframe and bins in which odds ratio will be calculated
    all_breaks=unique(na.omit(repeat_size_data[,LongAllele_col]))
    all_breaks=all_breaks[order(all_breaks)]
    # print(all_breaks)
    make_histogram_df=function(repeat_size_data,LongAllele_col,breaks,name){
        h=hist(repeat_size_data[,LongAllele_col],breaks=breaks,plot=F)
        breaks=h$breaks
        breaks=breaks[-1]
        counts=h$counts
        hist_data=data.frame(breaks,counts,stringsAsFactors=F)
        hist_data$lt=cumsum(hist_data$counts)
        hist_data$gt=rev(cumsum(rev(hist_data$counts)))
        hist_data$category=name
        return(hist_data)
    }
    hist_data_controls=make_histogram_df(controls,LongAllele_col,all_breaks,'Controls')
    hist_data_unsolved=make_histogram_df(unsolved_cases,LongAllele_col,all_breaks,'Unsolved Probands')
    hist_data_solved=make_histogram_df(solved_cases,LongAllele_col,all_breaks,'Solved Probands')
    as_far_as=length(all_breaks)-1
    # construct odds ratio table
    odds_ratio_data1=data.frame()
    odds_ratio_data2=data.frame()
    for(i in 1:as_far_as){
        test1=fisher.test(matrix(
            c(  hist_data_unsolved$gt[i],hist_data_controls$gt[i],
                hist_data_unsolved$lt[i],hist_data_controls$lt[i]),
                byrow=T,nrow=2))
        test2=fisher.test(matrix(
            c(  hist_data_unsolved$gt[i],hist_data_solved$gt[i],
                hist_data_unsolved$lt[i],hist_data_solved$lt[i]),
                byrow=T,nrow=2))
        pval1=-log10(test1$p.value)
        pval2=-log10(test2$p.value)
        ci_lower1=test1$conf.int[1]
        ci_lower2=test2$conf.int[1]
        ci_upper1=test1$conf.int[2]
        ci_upper2=test2$conf.int[2]
        df1=rbind(
            cbind(all_breaks[i],test1$estimate,ci_lower1,ci_upper1,'Odds Ratio Unsolved Probands vs controls'),
            cbind(all_breaks[i],pval1,NA,NA,'-log10 P-value')
        )
        df2=rbind(
            cbind(all_breaks[i],test2$estimate,ci_lower1,ci_upper1,'Odds Ratio Unsolved vs Solved Probands'),
            cbind(all_breaks[i],pval2,NA,NA,'-log10 P-value')
        )
        odds_ratio_data1=rbind(odds_ratio_data1,df1)
        odds_ratio_data2=rbind(odds_ratio_data2,df2)
    }
    colnames(odds_ratio_data1)=c('breaks','number','CI_lower','CI_upper','metric');
    odds_ratio_data1=data.frame(odds_ratio_data1,stringsAsFactors=F);
    # write.table(odds_ratio_data1,'odds_ratio_data.txt',sep='\t',col.names=F,row.names=F,quote=F)
    colnames(odds_ratio_data2)=c('breaks','number','CI_lower','CI_upper','metric');
    odds_ratio_data2=data.frame(odds_ratio_data2,stringsAsFactors=F);
    if(!is.null(handle_infinite_values)){
        odds_ratio_data1$number=as.numeric(odds_ratio_data1$number)
        if(any(is.infinite(odds_ratio_data1$number))){
            maximum=max(odds_ratio_data1$number[!is.infinite(odds_ratio_data1$number)])
            if(handle_infinite_values=='set to maximum'){
                odds_ratio_data1$number[is.infinite(odds_ratio_data1$number)]=maximum
            }else{
                odds_ratio_data1$number[is.infinite(odds_ratio_data1$number)]=0
            }
        }
        odds_ratio_data2$number=as.numeric(odds_ratio_data2$number)
        if(any(is.infinite(odds_ratio_data2$number))){
            maximum=max(odds_ratio_data2$number[!is.infinite(odds_ratio_data2$number)])
            if(handle_infinite_values=='set to maximum'){
                odds_ratio_data2$number[is.infinite(odds_ratio_data2$number)]=maximum
            }else{
                odds_ratio_data2$number[is.infinite(odds_ratio_data2$number)]=0
            }
        }
        # print(odds_ratio_data1)
        if(any(is.infinite(odds_ratio_data1$CI_upper))){
            print(odds_ratio_data1)
            maximum=max(odds_ratio_data1$CI_upper[!is.infinite(odds_ratio_data1$CI_upper)])
            odds_ratio_data1$CI_upper[is.infinite(odds_ratio_data1$CI_upper)]=0
        }
        odds_ratio_data2$CI_upper=as.numeric(odds_ratio_data2$CI_upper)
        if(any(is.infinite(odds_ratio_data2$CI_upper))){
            maximum=max(odds_ratio_data2$CI_upper[!is.infinite(odds_ratio_data2$CI_upper)])

            odds_ratio_data2$CI_upper[is.infinite(odds_ratio_data2$CI_upper)]=0
        }
    }
    print(odds_ratio_data1)
    write.table(cbind(odds_ratio_data1,odds_ratio_data2),paste0(gene,'_OR_table.tsv'),sep='\t',col.names=T,row.names=F,quote=F)
    #plot data
    hist_data=rbind(hist_data_controls,hist_data_unsolved,hist_data_solved)
    g_hist=ggplot(hist_data,aes(x=breaks,y=counts,fill=category,col=category))+geom_col()+
        scale_y_log10()+ggtitle(title)+xlab('')+ylab('Count')+theme(legend.position=legend_position)
    # g_hist=ggplot(d,aes(x=A2,col=case_solved_family))+
    #     geom_histogram()+scale_y_log10()
    g1=ggplot(odds_ratio_data1,aes(x=as.numeric(breaks),y=as.numeric(number)))+
        geom_point()+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Repeat sizes')+
        ylab(' ')+ggtitle('Unsolved cases vs Controls')+xlim(0,max(hist_data$breaks,na.rm=T))+
        facet_wrap(~factor(metric,levels=c('Odds Ratio Unsolved Probands vs controls','-log10 P-value')),ncol=1,scales='free_y',strip.position='left')
    g2=ggplot(odds_ratio_data2,aes(x=as.numeric(breaks),y=as.numeric(number)))+
        geom_point()+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Repeat sizes')+
        ylab(' ')+ggtitle('Unsolved cases vs Solved')+xlim(0,max(hist_data$breaks,na.rm=T))+
        facet_wrap(~factor(metric,levels=c('Odds Ratio Unsolved vs Solved Probands','-log10 P-value')),ncol=1,scales='free_y',strip.position='left')
    if(include_ci){
        g1=g1+geom_errorbar(aes(ymin=as.numeric(CI_lower), ymax=as.numeric(CI_upper)))
        g2=g2+geom_errorbar(aes(ymin=as.numeric(CI_lower), ymax=as.numeric(CI_upper)))
    }
    g=ggpubr::ggarrange(g_hist,g1,g2,ncol = 1, nrow = 3)
    return(g)
}





sharp_plot_delia_latest=function(repeat_size_data,
                    gene='PURA',gene_col='V2',
                    LongAllele_col='A2',breaks=50,case_solved_col='solved',
                    control_data,
                    include_ci=TRUE,handle_infinite_values='set to zero',start_from=0,log10=FALSE,export_used_data=FALSE){
    # read and organise repeat size table
    print(paste(basename(repeat_size_data),'vs',basename(control_data)))
    repeat_size_data=data.frame(fread(repeat_size_data),stringsAsFactors=F)
    repeat_size_data=repeat_size_data[repeat_size_data[,gene_col]==gene,]
    repeat_size_data=repeat_size_data[!is.na(repeat_size_data[,LongAllele_col]),]
    repeat_size_data[,LongAllele_col]=as.numeric(repeat_size_data[,LongAllele_col])
    repeat_size_data=repeat_size_data[as.numeric(repeat_size_data[,LongAllele_col])>=start_from,]
    repeat_size_data[is.na(repeat_size_data[,case_solved_col]),case_solved_col]='no'
    # print(repeat_size_data[,case_solved_col])
    repeat_size_data[repeat_size_data[,case_solved_col]!='yes',case_solved_col]='no'
    # cases solved column is supplied via the `case_solved_family_data` table
    unsolved_cases=repeat_size_data[repeat_size_data[,case_solved_col]=='no',]
    unsolved_cases=unsolved_cases[!is.na(unsolved_cases[,case_solved_col]),]
    solved_cases=na.omit(repeat_size_data[repeat_size_data[,case_solved_col]=='yes',])
    solved_cases=solved_cases[!is.na(solved_cases[,case_solved_col]),]
    # if control list is supplied- use it to separate controls in repeat data
    controls=data.frame(fread(control_data),stringsAsFactors=F)
    controls=controls[!is.na(controls[,LongAllele_col]),]

    # print(head(unsolved_cases))
    # print(head(solved_cases))
    # print(head(controls))
    # controls=repeat_size_data[repeat_size_data[,platekey_col]%in%controls,]
    controls[,case_solved_col]=NA
    repeat_size_data=repeat_size_data[,colnames(repeat_size_data)%in%colnames(controls)]
    controls=controls[,(colnames(controls)%in%colnames(repeat_size_data))]
    repeat_size_data=rbind(repeat_size_data,controls)
    # write.table(repeat_size_data,'tmp.tsv',sep='\t',row.names=F,quote=F)
    # print(tail())
    repeat_size_data=repeat_size_data[!is.na(repeat_size_data[,LongAllele_col]),]
    # print(repeat_size_data)
    print(paste(nrow(repeat_size_data[!is.na(repeat_size_data[,case_solved_col]),]),'cases'))
    print(paste(nrow(unsolved_cases),'unsolved cases'))
    print(paste(nrow(solved_cases),'solved cases'))
    print(paste(nrow(controls),'controls'))
    # print(unique(repeat_size_data[,LongAllele_col]))
    if(export_used_data){
        # unsolved_cases$category='unsolved cases'
        # solved_cases$category='solved cases'
        # controls$category='controls'
        export=repeat_size_data
        # print(head(repeat_size_data))
        write.table(export[order(as.numeric(export[,LongAllele_col]),decreasing=T),],paste0(gene,'_used_data.tsv'),sep='\t',row.names=F,col.names=F,quote=F)
    }


    breaks=min(as.numeric(repeat_size_data[,LongAllele_col]),na.rm=T):max(as.numeric(repeat_size_data[,LongAllele_col]),na.rm = T)
    print(breaks)
    output_table=data.frame()
    for(repeat_size in breaks){
        count_unsolved=nrow(unsolved_cases[as.numeric(unsolved_cases[,LongAllele_col])==repeat_size,])
        # print(head(unsolved_cases[as.numeric(unsolved_cases[,LongAllele_col])==repeat_size,]))
        lt_unsolved=nrow(unsolved_cases[as.numeric(unsolved_cases[,LongAllele_col])<repeat_size,])
        gt_unsolved=nrow(unsolved_cases[as.numeric(unsolved_cases[,LongAllele_col])>=repeat_size,])
        # counts_data=rbind(counts_data,cbind('unsolved',count,lt,gt,repeat_size))
        #solved
        count_solved=nrow(solved_cases[as.numeric(solved_cases[,LongAllele_col])==repeat_size,])
        lt_solved=nrow(solved_cases[as.numeric(solved_cases[,LongAllele_col])<repeat_size,])
        gt_solved=nrow(solved_cases[as.numeric(solved_cases[,LongAllele_col])>=repeat_size,])
        # counts_data=rbind(counts_data,cbind('solved',count,lt,gt,repeat_size))
        #controls
        count_controls=nrow(controls[as.numeric(controls[,LongAllele_col])==repeat_size,])
        lt_controls=nrow(controls[as.numeric(controls[,LongAllele_col])<repeat_size,])
        gt_controls=nrow(controls[as.numeric(controls[,LongAllele_col])>=repeat_size,])
        # counts_data=rbind(counts_data,cbind('controls',count,lt,gt,repeat_size))
        # print(gt_controls)
        # print(gt_unsolved)
        test1=fisher.test(matrix(
            c(  gt_unsolved,gt_controls,
                lt_unsolved,lt_controls),
                byrow=T,nrow=2))
        test2=fisher.test(matrix(
            c(  gt_unsolved,gt_solved,
                lt_unsolved,lt_solved),
                byrow=T,nrow=2))
        pval1=test1$p.value
        pval2=test2$p.value
        if(pval1<0.05){
            print(paste('Repeat size',repeat_size,'has a significant OR of',test1$estimate,'; P-value of',pval1));
        }
        if(pval2<0.05){
            print(paste('Repeat size',repeat_size,'has a significant OR of',test2$estimate,'; P-value of',pval2));
        }
        ci_lower1=test1$conf.int[1]
        ci_lower2=test2$conf.int[1]
        ci_upper1=test1$conf.int[2]
        ci_upper2=test2$conf.int[2]
        output_table=rbind(output_table,cbind(repeat_size,count_controls,lt_controls,gt_controls,
                                               count_unsolved,lt_unsolved,gt_unsolved,test1$estimate,pval1,ci_lower1,ci_upper1,
                                               count_solved,lt_solved,gt_solved,test2$estimate,pval2,ci_lower2,ci_upper2))
    }
    colnames(output_table)=c('repeat_size','controls_count','cumsum_controls','bigger_than_controls',
                            'count_unsolved','cumsum_unsolved','bigger_than_unsolved','OR_unsolved_vs_controls','Pval_unsolved_vs_controls','CI_lower_unsolved_vs_controls','CI_upper_unsolved_vs_controls',
                            'count_solved','cumsum_solved','bigger_than_solved','OR_unsolved_vs_solved','Pval_unsolved_vs_solved','CI_lower_unsolved_vs_solved','CI_upper_unsolved_vs_solved')
    tmp1=cbind(output_table[,c('repeat_size','controls_count')],'controls')
    colnames(tmp1)=c('repeat_size','counts','category')
    tmp2=cbind(output_table[,c('repeat_size','count_unsolved')],'unsolved')
    colnames(tmp2)=c('repeat_size','counts','category')
    tmp3=cbind(output_table[,c('repeat_size','count_solved')],'solved')
    colnames(tmp3)=c('repeat_size','counts','category')
    hist_data=rbind(tmp1,tmp2,tmp3)
    colnames(hist_data)=c('repeat_size','counts','category')
    g_hist=ggplot(hist_data,aes(x=repeat_size,y=counts,fill=category,col=category))+geom_col()+
        xlab('')+ylab('Count')+theme(legend.position=c(.75,.75))
    if(log10){g_hist=g_hist+scale_y_log10()}
    tmp1=cbind(output_table[,c('repeat_size','OR_unsolved_vs_controls','CI_lower_unsolved_vs_controls','CI_upper_unsolved_vs_controls')],'Odds Ratio Unsolved Probands vs controls')
    colnames(tmp1)=c('repeat_size','number','CI_lower','CI_upper','metric')
    tmp2=cbind(output_table[,c('repeat_size','Pval_unsolved_vs_controls')],NA,NA,'-log10 P-value')
    tmp2[,'Pval_unsolved_vs_controls']=-log10(tmp2[,'Pval_unsolved_vs_controls'])
    colnames(tmp2)=c('repeat_size','number','CI_lower','CI_upper','metric')
    odds_ratio_data1=rbind(tmp1, tmp2)
    colnames(odds_ratio_data1)=c('repeat_size','number','CI_lower','CI_upper','metric')
    tmp1=cbind(output_table[,c('repeat_size','OR_unsolved_vs_solved','CI_lower_unsolved_vs_solved','CI_upper_unsolved_vs_solved')],'Odds Ratio Unsolved vs Solved Probands')
    colnames(tmp1)=c('repeat_size','number','CI_lower','CI_upper','metric')
    tmp2=cbind(output_table[,c('repeat_size','Pval_unsolved_vs_solved')],NA,NA,'-log10 P-value')
    tmp2[,'Pval_unsolved_vs_solved']=-log10(tmp2[,'Pval_unsolved_vs_solved'])
    colnames(tmp2)=c('repeat_size','number','CI_lower','CI_upper','metric')
    odds_ratio_data2=rbind(tmp1,tmp2)
    colnames(odds_ratio_data2)=c('repeat_size','number','CI_lower','CI_upper','metric')
    odds_ratio_data1$number=as.numeric(odds_ratio_data1$number)
    # if(any(is.infinite(odds_ratio_data1$number))){
    #     maximum=max(odds_ratio_data1$number[!is.infinite(odds_ratio_data1$number)])
    #     print(maximum)
    #     odds_ratio_data1$number[is.infinite(odds_ratio_data1$number)]=0
    # }
    odds_ratio_data2$number=as.numeric(odds_ratio_data2$number)
    if(!is.null(handle_infinite_values)){
        odds_ratio_data1$number=as.numeric(odds_ratio_data1$number)
        if(any(is.infinite(odds_ratio_data1$number))){
            maximum=max(odds_ratio_data1$number[!is.infinite(odds_ratio_data1$number)])
            if(handle_infinite_values=='set to maximum'){
                odds_ratio_data1$number[is.infinite(odds_ratio_data1$number)]=maximum
            }else{
                odds_ratio_data1$number[is.infinite(odds_ratio_data1$number)]=0
            }
        }
        odds_ratio_data2$number=as.numeric(odds_ratio_data2$number)
        if(any(is.infinite(odds_ratio_data2$number))){
            maximum=max(odds_ratio_data2$number[!is.infinite(odds_ratio_data2$number)])
            if(handle_infinite_values=='set to maximum'){
                odds_ratio_data2$number[is.infinite(odds_ratio_data2$number)]=maximum
            }else{
                odds_ratio_data2$number[is.infinite(odds_ratio_data2$number)]=0
            }
        }
        # print(odds_ratio_data1)
        if(any(is.infinite(odds_ratio_data1$CI_upper))){
            # print(odds_ratio_data1)
            maximum=max(odds_ratio_data1$CI_upper[!is.infinite(odds_ratio_data1$CI_upper)])
            odds_ratio_data1$CI_upper[is.infinite(odds_ratio_data1$CI_upper)]=0
        }
        odds_ratio_data2$CI_upper=as.numeric(odds_ratio_data2$CI_upper)
        if(any(is.infinite(odds_ratio_data2$CI_upper))){
            maximum=max(odds_ratio_data2$CI_upper[!is.infinite(odds_ratio_data2$CI_upper)])

            odds_ratio_data2$CI_upper[is.infinite(odds_ratio_data2$CI_upper)]=0
        }
    }
    g1=ggplot(odds_ratio_data1,aes(x=as.numeric(repeat_size),y=as.numeric(number)))+
        geom_point()+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Repeat sizes')+
        ylab(' ')+ggtitle('Unsolved cases vs Controls')+
        facet_wrap(~factor(metric,levels=c('Odds Ratio Unsolved Probands vs controls','-log10 P-value')),ncol=1,scales='free_y',strip.position='left')
    
    g2=ggplot(odds_ratio_data2,aes(x=as.numeric(repeat_size),y=as.numeric(number)))+
        geom_point()+geom_line()+
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
        xlab('Repeat sizes')+
        ylab(' ')+ggtitle('Unsolved cases vs Solved')+
        facet_wrap(~factor(metric,levels=c('Odds Ratio Unsolved vs Solved Probands','-log10 P-value')),ncol=1,scales='free_y',strip.position='left')
    if(include_ci){
        g1=g1+geom_errorbar(aes(ymin=as.numeric(CI_lower), ymax=as.numeric(CI_upper)))
        g2=g2+geom_errorbar(aes(ymin=as.numeric(CI_lower), ymax=as.numeric(CI_upper)))
    }
    return(list(hist_plot=g_hist,unsolved_vs_controls_plot=g1,unsolved_vs_solved_plot=g2,table=output_table))
}



