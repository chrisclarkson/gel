ExpansionHunterDenovo merge \
    --reference genome.fa \
    --manifest manifest_merge_gms_no_relations_europeans \
    --output-prefix /re_gecip/shared_allGeCIPs/FA_AT/EHdenovo/gms_and_100k_no_relations_europeans \
    --min-unit-len 3 --max-unit-len 12
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

# ls /re_gecip/shared_allGeCIPs/FA_AT/Chris_Clarkson/intellectual_disability_analysis/GMS/non_neuros/str_profiles/*json |awk -F'str_profiles/' '{
#   split($1, split_field, /.str_profile.json/)
#   print split_field[0]"\t"control"\t"$1".str_profile.json"
# }'
