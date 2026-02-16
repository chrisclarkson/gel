# code adapted from https://github.com/PacificBiosciences/apps-scripts/tree/master/RepeatAnalysisTools
python pacbio_repeat_analysis_tools/extractRegion.py /gel_data_resources/LRS_cohort_genomes/path_to_lrs.fa \
    /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa \
    'chrX:19990538-19991358' \
    python pacbio_repeat_analysis_tools/waterfall.py -m CCG,CCA -o waterfallBCLAF3 --format pdf

python pacbio_repeat_analysis_tools/extractRegion.py /gel_data_resources/LRS_cohort_genomes/path_to_lrs2.fa \
    /public_data_resources/reference/GRCh38/GRCh38Decoy_no_alt.fa \
    'chrX:19990538-19991358' \
    python pacbio_repeat_analysis_tools/waterfall.py -m CCG,CCA -o waterfallBCLAF3_control --format pdf

