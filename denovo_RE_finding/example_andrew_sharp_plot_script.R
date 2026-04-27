setwd('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots')
source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots_no_ggpubr.R')
output=sharp_plot_delia_no_ggpubr()
print(output$table)
print(output$unsolved_vs_controls_plot)
print(output$unsolved_vs_solved_plot)
print(output$hist_plot)


g=sharp_plot_delia_no_ggpubr(
        repeat_size_data='/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/PURA_neuro_100k_all.tsv',
        gene='PURA',gene_col='V2',
        LongAllele_col='A2',breaks=50,
        control_data='/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')

print(g$unsolved_vs_solved_plot)



g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'ataxia_patients_case_solved.txt','_PURA_EP400.vcf','Ataxia','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_ataxia.png')
print(g)
dev.off()
g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'ataxia_patients_case_solved.txt','_PURA_EP400.vcf','Ataxia','/re_gecip/neurology/DeliaG/PURA_EP400/100k/EP400/control_EP400_100k.tsv')
png('EP400_ataxia.png')
print(g)
dev.off()


g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'HSP_case_solved.txt','_PURA_EP400.vcf','HSP','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_HSP.png')
print(g)
dev.off()
g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'HSP_case_solved.txt','_PURA_EP400.vcf','HSP','/re_gecip/neurology/DeliaG/PURA_EP400/100k/EP400/control_EP400_100k.tsv')
png('EP400_HSP.png')
print(g)
dev.off()

g=sharp_plot_delia_no_ggpubr('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'ataxia_case_solved_gmcq.txt','_PURA_EP400.vcf','Ataxia','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_ataxia_gmcq.png')
print(g)
dev.off()
g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'ataxia_case_solved_gmcq.txt','_PURA_EP400.vcf','Ataxia','/re_gecip/neurology/DeliaG/PURA_EP400/100k/EP400/control_EP400_100k.tsv')
png('EP400_ataxia_gmcq.png')
print(g)
dev.off()

g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'HSP_case_solved_gmcq.txt','_PURA_EP400.vcf','HSP','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_HSP_gmcq.png')
print(g)
dev.off()
g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'HSP_case_solved_gmcq.txt','_PURA_EP400.vcf','HSP','/re_gecip/neurology/DeliaG/PURA_EP400/100k/EP400/control_EP400_100k.tsv')
png('EP400_HSP_gmcq.png')
print(g)
dev.off()




g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'intellectual_disability_case_solved.txt','_PURA_EP400.vcf','Intellectual disability','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'intellectual_disability_case_solved.txt','_PURA_EP400.vcf','Intellectual disability','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_ID.png')
print(g)
dev.off()

g=sharp_plot_delia('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'neurological_disorders_noID_case_solved_gmcq.txt','_PURA_EP400.vcf','PURA: Neurological disorders','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_neuro.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'neurological_disorders_noID_case_solved.txt','_PURA_EP400.vcf','PURA: Neurological disorders (minus Intellectual disability)','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_neuro_noID.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'CMT_disease_case_solved.txt','_PURA_EP400.vcf','PURA: Charcot-Marie tooth disease','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_CMT_disease.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'CMT_HN_PMN_ALS_panel_case_solved.txt','_PURA_EP400.vcf','PURA: CMT+HN+PMN+ALS','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_neuro_CMT_HN_PMN_ALS_panel.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'CMT_HN_PMN_panel_case_solved.txt','_PURA_EP400.vcf','PURA: CMT+HN+PM','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_neuro_CMT_HN_PMN_panel.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','PURA','V2','V4',50,'CMT_PMN_case_solved.txt','_PURA_EP400.vcf','PURA: CMT+PMN','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('PURA_neuro_CMT_PMN_panel.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'neurological_disorders_case_solved.txt','_PURA_EP400.vcf','EP400 Neurological disorders','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('EP400_neuro.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'neurological_disorders_noID_case_solved.txt','_PURA_EP400.vcf','EP400 Neurological disorders (minus Intellectual disability)','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('EP400_neuro_noID.png')
p1       rint(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'neurological_disorders_noID_panel_case_solved.txt','_PURA_EP400.vcf','EP400 Neurological disorders (minus Intellectual disability)','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('EP400_neuro_noID_panel.png')
print(g)
dev.off()

g=sharp_plot('../100k_PURA_EP400_reps.tsv','V1','EP400','V2','V4',50,'intellectual_disability_case_solved.txt','_PURA_EP400.vcf','Intellectual disability','/re_gecip/neurology/DeliaG/PURA_EP400/100k/PURA/control_PURA_100k.tsv')
png('EP400_ID.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/CC_AT/EH/RNF213/expansion_sizes_omit_NAs__same_line.txt','V1','KCTD1','V2','V4',50,'ataxia_patients_case_solved.txt','_RNF213','')
png('RNF213.png')
print(g)
dev.off()
source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/CC_AT/EH/RBFOX1/expansion_sizes_omit_NAs__same_line.txt','V1','RBFOX1','V2','V4',50,'ataxia_patients_case_solved.txt','_RBFOX1','')
png('RBFOX1_2.png')
print(g)
dev.off()
source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/CC_AT/EH/RBFOX1_smaller/expansion_sizes_omit_NAs__same_line.txt','V1','RBFOX1','V2','V4',50,'ataxia_patients_case_solved.txt','_RBFOX1_smaller','')
png('RBFOX1_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('expansion_sizes_omit_NAs_FRMPD4_larger_same_line.txt','V1','FRMPD4','V2','V4',50,'intellectual_disability_case_solved.txt','_FRMPD4_larger','',legend_position=c(.2,.75),include_ci=FALSE,export_used_data=T)
png('FRMPD4_larger_aaarg_larger.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('expansion_sizes_omit_NAs_FRMPD4_larger_aaag_same_line.txt',
        'V1','FRMPD4','V2','V3',50,'intellectual_disability_case_solved.txt','_FRMPD4_larger_aaag','')
png('FRMPD4_largeraaag_smaller.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/EH/FRMPD4_larger/expansion_sizes_omit_NAs_FRMPD4_sex_same_line.txt',
        'V1','FRMPD4','V2','V3',50,'intellectual_disability_case_solved.txt','_FRMPD4_larger','',include_ci=FALSE,export_used_data=T)
png('FRMPD4_larger_sex_used_a1_irrs.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/EH/FRMPD4/expansion_sizes_omit_NAs__same_line.txt',
        'V1','FRMPD4','V2','V14',50,'intellectual_disability_case_solved.txt','_FRMPD4','',include_ci=FALSE,export_used_data=T)
png('FRMPD4_sex_used_a1_irrs.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('FRMPD4_larger_sex_used_males_only.tsv',
        'V1','FRMPD4','V2','V3',50,'intellectual_disability_case_solved.txt','_FRMPD4_larger','',include_ci=FALSE,export_used_data=T)
png('FRMPD4_larger_sex_males_only_used_a1.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/THAP11/expansion_sizes_omit_NAs_THAP11_same_line.txt','V1',
      'THAP11_CAG','V2','V4',50,
      'neurological_disorders_case_solved.txt','_THAP11','',
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('THAP11_ND.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/ABCC1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'ABCC1','V2','V4',50,
      'intellectual_disability_case_solved.txt','_ABCC1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ABCC1_ID.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/ABCC1_longer/expansion_sizes_omit_NAs__same_line.txt','V1',
      'ABCC1_longer','V2','V4',50,
      'intellectual_disability_case_solved.txt','_ABCC1_longer','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ABCC1_longer_ID.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/ABCC1_os/expansion_sizes_omit_NAs__same_line.txt','V1',
      'ABCC1_os','V2','V3',50,
      'intellectual_disability_case_solved.txt','_ABCC1_os','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ABCC1_os_ID_shorter_allele.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/SHANK2/expansion_sizes_omit_NAs__same_line.txt','V1',
      'SHANK2','V2','V4',50,
      'ataxia_patients_case_solved.txt','_SHANK2','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('SHANK2_ataxia.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/SHANK2/expansion_sizes_omit_NAs_shank2_gms_included_same_line.txt','V1',
      'SHANK2','V2','V3',50,
      'ataxia_patients_case_solved.txt','_SHANK2','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('SHANK2_ataxia_gms_included.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/SHANK2_aaatg/expansion_sizes_omit_NAs_SHANK2_aaatg_same_line.txt','V1',
      'SHANK2_aaatg','V2','V4',50,
      'ataxia_patients_case_solved.txt','_SHANK2_aaatg','',export_used_data=T,
      handle_infinite_values='set to maximum',include_ci=FALSE)

png('SHANK2_aaatg_ataxia_gms_cancer_included2.png')
print(g)
dev.off()




source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('expansion_sizes_omit_NAs_SHANK2_same_line_aaatg.txt','V1',
      'SHANK2','V2','V4',50,
      'ataxia_patients_case_solved.txt','_SHANK2','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('SHANK2_ataxia_gms_included_aaatg.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/ATG10/expansion_sizes_omit_NAs__same_line.txt','V1',
      'ATG10','V2','V4',50,
      'ataxia_patients_case_solved.txt','_ATG10','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ATG10_ataxia_gms_included_tttkc.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('atg10_tttgc_present.txt','V1',
      'ATG10','V2','V4',50,
      'ataxia_patients_case_solved.txt','_ATG10','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ATG10_ataxia_gms_included_tttgc.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/ATG10/expansion_sizes_omit_NAs__same_line.txt','V1',
      'ATG10','V2','V4',50,
      'neurological_disorders_case_solved.txt','_ATG10','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('ATG10_ND_gms_included_tttkc.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/BIN1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'BIN1','V2','V4',50,
      'neurological_disorders_case_solved.txt','_BIN1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('BIN1_ND_gms_and_cancer_included.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/BIN1/expansion_sizes_omit_NAs__eh_classifier_procesed_data.txt','SampleId',
      'BIN1','VARID','LongAllele',50,
      'ataxia_patients_case_solved.txt','_BIN1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE,eh_classifier_col='PredictedLabels')

png('BIN1_ataxia_gms_included_cancer_eh_classified.png')
print(g)
dev.off()





source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/COLEC11/expansion_sizes_omit_NAs__same_line.txt','V1',
      'COLEC11','V2','V4',50,
      'ataxia_patients_case_solved.txt','_COLEC11','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('COLEC11_ataxia_gms_included_cancer.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/TTC1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'TTC1','V2','V4',50,
      'ataxia_patients_case_solved.txt','_TTC1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('TTC1_ataxia_gms_included_cancer.png')
print(g)
dev.off()




source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/RC/output/COLEC11_african.txt',
      'V1','COLEC11','V2','V4',50,'ataxia_patients_case_solved.txt','_COLEC11','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('COLEC11_ataxia_african_gms_included_cancer.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MBOAT7/expansion_sizes_omit_NAs__same_line.txt','V1',
      'MBOAT7','V2','V4',50,
      'ataxia_patients_case_solved.txt','_MBOAT7','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('MBOAT7_ataxia_gms_included_cancer.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MBOAT7/expansion_sizes_omit_NAs__eh_classifier_procesed_data.txt','SampleId',
      'MBOAT7','VARID','LongAllele',50,
      'ataxia_patients_case_solved.txt','_MBOAT7','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE,eh_classifier_col='PredictedLabels')

png('MBOAT7_ataxia_gms_included_cancer_eh_classified.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/SMARCA1/expansion_sizes_omit_NAs__eh_classifier_procesed_data.txt','SampleId',
      'SMARCA1','VARID','LongAllele',50,
      'ataxia_patients_case_solved.txt','_SMARCA1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE,eh_classifier_col='PredictedLabels')

png('SMARCA1_ataxia_gms_included_cancer_eh_classified.png')
print(g)
dev.off()





source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MBOAT7/expansion_sizes_omit_NAs__same_line.txt','V1',
      'MBOAT7','V2','V4',50,
      'ataxia_patients_case_solved.txt','_MBOAT7','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('MBOAT7_ND_gms_included_cancer.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/VPS13B/expansion_sizes_omit_NAs__same_line.txt','V1',
      'VPS13B','V2','V4',50,
      'ataxia_patients_case_solved.txt','_VPS13B','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('VPS13B_ataxia_gms_included_cancer.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/OCM/expansion_sizes_omit_NAs__same_line.txt','V1',
      'OCM','V2','V3',50,
      'ataxia_patients_case_solved.txt','_OCM','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('OCM_smaller_ataxia_gms_included_cancer.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/SMARCA1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'SMARCA1','V2','V3',50,
      'ataxia_patients_case_solved.txt','_SMARCA1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('SMARCA1_smaller_ataxia_gms_included.png')
print(g)
dev.off()

data=fread('../EH/SMARCA1/expansion_sizes_omit_NAs__same_line.txt')
data$V1=basename(gsub('_SMARCA1.vcf','',data$V1))
meta=fread('../combined_jobs.tsv',fill=T)
data=merge(data,meta,by='V1')
write.table(data[tolower(data$V5.y)=='female',],'smarca1_females.tsv',sep='\t',row.names=F,col.names=F,quote=F)
write.table(data[tolower(data$V5.y)=='male',],'smarca1_males.tsv',sep='\t',row.names=F,col.names=F,quote=F)


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('smarca1_females.tsv','V1',
      'SMARCA1','V2','V4',50,
      'ataxia_patients_case_solved.txt','_SMARCA1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('SMARCA1_ataxia_gms_included_females.png')
print(g)
dev.off()
file.rename('SMARCA1_used_data.tsv','SMARCA1_used_data_females.tsv')

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('smarca1_males.tsv','V1',
      'SMARCA1','V2','V4',50,
      'ataxia_patients_case_solved.txt','_SMARCA1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('SMARCA1_ataxia_gms_included_males.png')
print(g)
dev.off()
file.rename('SMARCA1_used_data.tsv','SMARCA1_used_data_males.tsv')


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/COL1A1_ot/expansion_sizes_omit_NAs__same_line.txt','V1',
      'COL1A1_ot','V2','V4',50,
      'ataxia_patients_case_solved.txt','_COL1A1_ot','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('COL1A1_ataxia_gms_included.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('COL1A1_ot_all.txt','V1',
      'COL1A1_ot_all','V2','V4',50,
      'ataxia_patients_case_solved.txt','_COL1A1_ot_all','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('COL1A1_ot_all_ataxia_gms_included.png')
print(g)
dev.off()





source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/THAP11/expansion_sizes_omit_NAs_THAP11_same_line.txt',
        'V1','THAP11_CAG','V2','V4',50,
        'neurological_disorders_case_solved.txt','_THAP11','',
        legend_position=c(.2,.75),include_ci=FALSE)
png('THAP11_ND_new.png')
print(g)
dev.off()






source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/EH/DISC1FP1/expansion_sizes_omit_NAs__same_line.txt',
      'V1','DISC1FP1','V2',
      'V4',50,'intellectual_disability_case_solved.txt','_DISC1FP1','',
      include_ci=F,handle_infinite_values='set to maximum',export_used_data=T)
png('DISC1FP1.png')
print(g)
dev.off()




source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/EH/SHROOM2/expansion_sizes_omit_NAs_SHROOM2_same_line.txt','V1','SHROOM2','V2','V3',50,'intellectual_disability_case_solved.txt','_SHROOM2','',include_ci=F,handle_infinite_values='set to maximum',start_from=30)
png('SHROOM2.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/EH/ADARB1/expansion_sizes_omit_NAs_ADARB1_same_line.txt','V1','ADARB1','V2','V4',50,'intellectual_disability_case_solved.txt','_ADARB1','',include_ci=F,handle_infinite_values='set to maximum',start_from=30)
png('ADARB1.png')
print(g)
dev.off()


source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/EH/USH2A/expansion_sizes_omit_NAs_USH2A_same_line.txt','V1','USH2A',
'V2','V4',50,'intellectual_disability_case_solved.txt','_USH2A','',include_ci=F,handle_infinite_values='set to maximum',start_from=30)
png('USH2A.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/EH/NPAS1NM/expansion_sizes_omit_NAs__same_line.txt',
        'V1','NPAS1NM','V2','V4',50,'intellectual_disability_case_solved.txt',
        '_NPAS1NM','',include_ci=F,handle_infinite_values='set to maximum',start_from=0,export_used_data=T)
png('NPAS1NM.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/EH/NPAS1NM/expansion_sizes_omit_NAs_NPAS1NM_same_line.txt',
        'V1','NPAS1NM','V2','V3',50,'intellectual_disability_case_solved.txt',
        '_NPAS1NM','',include_ci=F,handle_infinite_values='set to maximum',start_from=0,export_used_data=T)
png('NPAS1NM_smaller.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/STX7_SAS_only.tsv',
        'V1','STX7','V2','V4',50,'intellectual_disability_case_solved.txt',
        '_STX7','',include_ci=F,handle_infinite_values='set to maximum',start_from=0,export_used_data=T)
png('STX7_SAS_only.png')
print(g)
dev.off()






source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MYOCD/expansion_sizes_omit_NAs_MYOCD_same_line.txt','V1','MYOCD','V2','V4',50,'neurological_disorders_case_solved.txt','_MYOCD','',include_ci=F,handle_infinite_values='set to maximum',start_from=0,export_used_data=T)
png('MYOCD_ND.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MYOCD/expansion_sizes_omit_NAs_MYOCD_same_line.txt','V1','MYOCD','V2','V3',50,'neurological_disorders_case_solved.txt','_MYOCD','',include_ci=F,handle_infinite_values='set to maximum',start_from=0)
png('MYOCD_ND_small_allele.png')
print(g)
dev.off()



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/EH/TMEM52_wtf/expansion_sizes_omit_NAs_TMEM52_wtf_same_line.txt',
        'V1','TMEM52','V2','V4',50,'intellectual_disability_case_solved.txt',
        '_TMEM52_wtf','',include_ci=F,handle_infinite_values='set to maximum',
        start_from=0,export_used_data=T)

png('TMEM52.png')
print(g)
dev.off()

setwd('/re_gecip/neurology/Kristina/research/MIAMI/chris_tables/')
source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
for(gene in c('BAIAP','ADARB1','BCLAF3','E2F4','KCNN3','MEF2A','ZNF384')){
        g=sharp_plot_delia_latest(
                repeat_size_data=paste0('/re_gecip/neurology/Kristina/research/MIAMI/chris_tables/',gene,'_100K_and_GMS_cases.tsv'),
                control_data=paste0('/re_gecip/neurology/Kristina/research/MIAMI/chris_tables/',gene,'_100K_and_GMS_controls.tsv'),
                gene,'gene',LongAllele_col='a2',case_solved_col='Case.Solved.Family',
                include_ci=F,handle_infinite_values='set to maximum',log10=T,
                start_from=0,export_used_data=TRUE)
        png(paste0(gene,'_log.png'))
        print(ggarrange(g$hist_plot,g$unsolved_vs_solved_plot,g$unsolved_vs_controls_plot,ncol=1))
        dev.off()
}



source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/CLEC2B/expansion_sizes_omit_NAs__same_line.txt','V1',
      'CLEC2B','V2','V3',50,
      'ataxia_patients_case_solved.txt','_CLEC2B','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('CLEC2B_ataxia_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/LITAF/expansion_sizes_omit_NAs__same_line.txt','V1',
      'LITAF','V2','V4',50,
      'ataxia_patients_case_solved.txt','_LITAF','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('LITAF_ataxia_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/LITAF_aaact/expansion_sizes_omit_NAs__same_line.txt','V1',
      'LITAF','V2','V4',50,
      'ataxia_patients_case_solved.txt','_LITAF_aaact','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('LITAF_aaact_ataxia.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH_asian_population_specific/LITAF_aaact/expansion_sizes_omit_NAs__same_line.txt',
      'V1','LITAF','V2','V4',50,
      'ataxia_patients_case_solved.txt','_LITAF_aaact','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('LITAF_aaact_asian_ataxia.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/DAB1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'DAB1','V2','V3',50,
      'ataxia_patients_case_solved.txt','_DAB1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('DAB1_ataxia_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/DACH1/expansion_sizes_omit_NAs__same_line.txt','V1',
      'DACH1','V2','V3',50,
      'ataxia_patients_case_solved.txt','_DACH1','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('DACH1_ataxia_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/C1orf168/expansion_sizes_omit_NAs__same_line.txt','V1',
      'C1orf168','V2','V3',50,
      'ataxia_patients_case_solved.txt','_C1orf168','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)
png('C1orf168_ataxia_smaller.png')
print(g)
dev.off()

source('/re_gecip/shared_allGeCIPs/AT_CC/GMS_and_100K_combined/andrew_sharp_plots/andrew_sharp_plots.R')
g=sharp_plot_delia('../EH/MYH6/expansion_sizes_omit_NAs__same_line.txt','V1',
      'MYH6','V2','V4',50,
      'ataxia_patients_case_solved.txt','_MYH6','',export_used_data=T,
      handle_infinite_values='set to maximum',legend_position=c(.2,.75),include_ci=FALSE)

png('MYH6_ataxia_gms_included_cancer.png')
print(g)
dev.off()





