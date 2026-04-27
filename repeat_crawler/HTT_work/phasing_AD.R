# set working directory 

setwd("/nas/weka.gel.zone/re_gecip/shared_allGeCIPs/AD_CC/repeat_crawler")

# load datasets 

data_v2=read.table('v2_post_processing/fmtd_v2.tsv',sep='\t',
                   stringsAsFactors = F,header=T)

data_v0=read.table('v0_post_processing/fmtd_v0.tsv',sep='\t',
                   stringsAsFactors = F,header=T)


#RCv2 data:
eha1_rca1=data_v2[data_v2$EH_A1 == data_v2$RC_A1,]

eha2_rca2=data_v2[data_v2$EH_A2==data_v2$RC_A2,]

eha1_rca2=data_v2[data_v2$EH_A1 == data_v2$RC_A2,]
  
eha2_rca1=data_v2[data_v2$EH_A2 == data_v2$RC_A1,]
  

eha1_rca2=eha1_rca2[!(eha1_rca2$bam_file%in%eha1_rca1$bam_file),]

eha1_rca2=eha1_rca2[!is.na(eha1_rca2$bam_file),]

eha2_rca2=eha2_rca2[!(eha2_rca2$bam_file%in%eha1_rca1$bam_file) &
                      !(eha2_rca2$bam_file%in%eha1_rca2$bam_file),]


eha2_rca2=eha2_rca2[!is.na(eha2_rca2$bam_file),]

eha2_rca1=eha2_rca1[!(eha2_rca1$bam_file%in%eha1_rca1$bam_file) &
                      !(eha2_rca1$bam_file%in%eha1_rca2$bam_file) &
                      !(eha2_rca1$bam_file %in% eha2_rca2$bam_file),]

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

data_v2_sorted<-rbind(eha1_rca1, eha1_rca2, eha2_rca1, 
                      eha2_rca2)

# put the allele types into categories 
#a1

data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CAACAG|CCGCCA|CCG|CCT)"]<-"canonical"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|(CAACAG)x2|CCGCCA|CCG|CCT)"]<-"Q2_dup"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CAACAG|CCG|CCT)"]<-"P1_loss"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CAACAG|CCGCCA|CCT)"]<-"P2_loss"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CAACAG|CCT)"]<-"P1_P2_dup"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CCG|CCT)"]<-"Q1_P1_loss"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CCGCCA|CCG|CCT)"]<-"Q2_loss"
data_v2_sorted$GT1_category[data_v2_sorted$GT1=="(CAG|CAA|CCGCCA|CCG|CCT)"]<-"partial_Q2_loss"

# add v0 genotypes 

colnames(data_v0)<-c("bam_file", "GT1_v0", "GT2_v0")

data_v2_v0<-merge(data_v2_sorted, data_v0,
                  by = "bam_file")

# Checking which v0 genotypes match the phased genotype

matching_g1=data_v2_v0[data_v2_v0$GT1==data_v2_v0$GT1_v0,]
matching_g2=data_v2_v0[data_v2_v0$GT1==data_v2_v0$GT2_v0,]
matching_g2=matching_g2[!(matching_g2$bam_file%in%matching_g1$bam_file),]

# Assigning the second unphased genotype
matching_g1$GT2=matching_g1$GT2_v0
matching_g2$GT2=matching_g2$GT1_v0

out=rbind(matching_g1,matching_g2)

# 53 observations did not match 
unmatched<-data_v2_v0[! data_v2_v0$bam_file %in% 
                        out$bam_file,]

# Giving categories to GT2
out$GT2_category[out$GT2=="(CAG|(CAACAG)x2|CCGCCA|CCG|CCT)" |
                   out$GT2=="(CAG|(CAACAG)x2|CCGCCA|CCG)"]<-"Q2_dup"
out$GT2_category[out$GT2=="(CAG|CAACAG|CCG|CCT)" |
                   out$GT2=="(CAG|CAACAG|CCG)"]<-"P1_loss"
out$GT2_category[out$GT2=="(CAG|CAACAG|CCGCCA|CCG|CCT)" |
                   out$GT2=="(CAG|CAACAG|CCGCCA|CCG)"]<-"canonical"

out$GT2_category[out$GT2=="(CAG|CAACAG|CCGCCA)"]<-"P2_loss"

out$GT2_category[out$GT2=="(CAG|CCG)" |
                   out$GT2=="(CAG|CCG|CCT)"]<-"Q2_P1_loss"

out$GT2_category[out$GT2=="(CAG|CCGCCA|CCG|CCT)"]<-"Q2_loss"

out$GT2_category[out$GT2=="(CAG|CAA|CCG|CCT)"]<-"partial_Q2_loss"

write.table(unmatched, file = "unmatched_AD.txt",
            sep = "\t", quote = F, row.names = F)

write.table(out, file = "final_phased_genotypes_AD.txt",
            sep = "\t", quote = F, row.names = F)

# For the "unmatched" individuals with homozyous calls from v0
unmatched_homo.gt.v0<-unmatched[unmatched$GT1_v0 == unmatched$GT2_v0,]

unmatched_homo.gt.v0$GT1 = unmatched_homo.gt.v0$GT1_v0
unmatched_homo.gt.v0$GT2 = unmatched_homo.gt.v0$GT2_v0


# Giving categories to GTs
unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|(CAACAG)x2|CCGCCA|CCG|CCT)" |
                                    unmatched_homo.gt.v0$GT2=="(CAG|(CAACAG)x2|CCGCCA|CCG)"]<-"Q2_dup"
unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CAACAG|CCG|CCT)" |
                                    unmatched_homo.gt.v0$GT2=="(CAG|CAACAG|CCG)"]<-"P1_loss"
unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CAACAG|CCGCCA|CCG|CCT)" |
                                    unmatched_homo.gt.v0$GT2=="(CAG|CAACAG|CCGCCA|CCG)"]<-"canonical"

unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CAACAG|CCGCCA)"]<-"P2_loss"

unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CCG)" |
                                    unmatched_homo.gt.v0$GT2=="(CAG|CCG|CCT)"]<-"Q2_P1_loss"

unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CCGCCA|CCG|CCT)"]<-"Q2_loss"

unmatched_homo.gt.v0$GT2_category[unmatched_homo.gt.v0$GT2=="(CAG|CAA|CCG|CCT)"]<-"partial_Q2_loss"

unmatched_homo.gt.v0$GT1_category<-unmatched_homo.gt.v0$GT2_category

out_2<-rbind(out, unmatched_homo.gt.v0)


# manually checking the single unmatched heterozygous call by v0
unmatched_het.gt.v0<-unmatched[unmatched$GT1_v0 != unmatched$GT2_v0,]

unmatched_het.gt.v0$GT1<-unmatched_het.gt.v0$GT2_v0
unmatched_het.gt.v0$GT1_category<-"canonical"

unmatched_het.gt.v0$GT2<-unmatched_het.gt.v0$GT1_v0
unmatched_het.gt.v0$GT2_category<-"Q2_P1_loss"

out_3<-rbind(out_2, unmatched_het.gt.v0)

#sort the final dataset to ensure all Q1_1 <= Q1_2
out_final_1<-out_3[out_3$Q1_1<=out_3$Q1_2,]

out_final_2<-out_3[out_3$Q1_1>out_3$Q1_2,]

out_final_2[, c("Q1_1", "Q1_2")] <- out_final_2[, c("Q1_2", "Q1_1")]
out_final_2[, c("GT1", "GT2")] <- out_final_2[, c("GT2", "GT1")]
out_final_2[, c("GT1_category", "GT2_category")] <- 
  out_final_2[, c("GT2_category", "GT1_category")]

out_final<-rbind(out_final_1, out_final_2)

