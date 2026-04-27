#!/bin/bash

# module load bio/illumina/ExpansionHunterDenovo/0.9.0

/tools/aws-workspace-apps/re_admin/source_code/ExpansionHunterDenovo/0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64/scripts/casecontrol.py locus \
    --manifest manifest_gms \
    --multisample-profile gms_and_100k_no_relations.multisample_profile.json \
    --output gms_and_100k_no_relations_casecontrol_locus_min5.tsv \
    --min-inrepeat-reads 2

/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/analysis_scripts/annotate_ehdn.sh \
    --ehdn-results gms_and_100k_no_relations_casecontrol_locus_min5.tsv \
    --ehdn-annotated-results gms_and_100k_no_relations_casecontrol_locus_min5_annotated.tsv \
    --annovar-annotate-variation ./annotate_variation.pl \
    --annovar-humandb /public_data_resources/annovar_databases/ \
    --annovar-buildver hg38



/tools/aws-workspace-apps/re_admin/source_code/ExpansionHunterDenovo/0.9.0/ExpansionHunterDenovo-v0.9.0-linux_x86_64/scripts/casecontrol.py locus \
    --manifest manifest_merge_gms_no_relations_europeans \
    --multisample-profile gms_and_100k_no_relations_europeans.multisample_profile.json \
    --output gms_and_100k_no_relations_europeans_casecontrol_locus_min5.tsv \
    --min-inrepeat-reads 2

/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/analysis_scripts/annotate_ehdn.sh \
    --ehdn-results gms_and_100k_no_relations_europeans_casecontrol_locus_min5.tsv \
    --ehdn-annotated-results gms_and_100k_no_relations_europeans_casecontrol_locus_min5_annotated.tsv \
    --annovar-annotate-variation ./annotate_variation.pl \
    --annovar-humandb /public_data_resources/annovar_databases/ \
    --annovar-buildver hg38





/re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/scripts/outlier.py locus \
  --manifest manifest_merge_gms \
  --multisample-profile gms_and_100k.multisample_profile.json \
  --output gms_and_100k_outlier_locus_min5.tsv

/re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/analysis_scripts/annotate_ehdn.sh \
    --ehdn-results gms_and_100k_outlier_locus_min5.tsv \
    --ehdn-annotated-results gms_and_100k_outlier_locus_min5_annotated.tsv \
    --annovar-annotate-variation ./annotate_variation.pl \
    --annovar-humandb /public_data_resources/annovar_databases/ \
    --annovar-buildver hg38


