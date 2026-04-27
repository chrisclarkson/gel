`find_regex_genotypes.py` is version 0 of repeat crawler. It accepts reads that span from 1-6 and only doesn't count the CAG stretch it was used to generate `most_likely_gts_parents.tsv` 

It can be run with the latest version of RC as follows:

```
python find_regex_genotypes_parallel_vcf.py --bam_files test_files_backups --span 1-6 --count CAA CAACAG "(CAACAG)x2" CCGCCA CCG CCT --output out_1-6_CAG-CCT_CAG_not_counted.tsv
```
