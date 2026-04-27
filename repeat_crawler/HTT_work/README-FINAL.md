
Importance: Short repeated sequences on the right flank of the CAG in the HTT locus are thought modify disease pathogenesis thorough molecular instability. 

------------------------
The structure of the HTT locus: 
 LeftFlank -- (CAG)n (CAACAG)0-2 (CCGCCA)0-1 (CCG)n (CCT)1-3 -- RightFlank

Also referred to as:
 LeftFlank -- Q1 Q2 P1 P3 P3 -- RightFlank
------------------------

Question: Can we genotype those sequences using short-read WGS (150bp), and can they be phased to the CAG tract?


While ExpansionHunter can accurately size the CAG tract up to 35 repeats, it struggles to detect polymorphisms in each component (CAACAG deletions, for example), and is not designed to phase alleles.

We here develop a workflow, Repeat Crawler (RC), to be run alongside ExpansionHunter (EH) to genotype and phase the short repeat sequences flanking the CAG (Q1) in the HTT locus. RC takes BAM files (output by EH) and summarises the content of complex repeats in the reads. 
It does so by taking a list of repeat component structures (for example., Q2, P1, P2 etc.) as input to return the length of each component (which can be 0 in cases of an absent repeat component). 

For a more in-depth description of how RC works- see README_RC_indepth.md.

Note that RC requires EH to be run with a custom configuration (input) file.
Here, the two json files used as input to EH have been designed for 2 genome builds: GRCh37 and GRCh38. 
See below the coordinates used for both files:

```
	GRCh37: 
	{
        "LocusId": "HTT",
        "LocusStructure": "(CAG)*(CAACAG)*(CCGCCA)*(CCG)*(CCT)*",
        "ReferenceRegion": [
            "4:3076603-3076660",
	          "4:3076661-3076666",
            "4:3076667-3076672",
            "4:3076673-3076693",
	          "4:3076694-3076699"
        ],
        "VariantId": [
            "HTT",
	          "HTT_CAACAG",
	          "HTT_CCGCCA",
            "HTT_CCG",
	          "HTT_CCT"
        ],
        "VariantType": [
            "Repeat",
            "Repeat",
            "Repeat",
	          "Repeat",
            "Repeat"
        ]
    }
    GRCh38:
    {
        "LocusId": "HTT",
        "LocusStructure": "(CAG)*(CAACAG)*(CCGCCA)*(CCG)*(CCT)*",
        "ReferenceRegion": [
            "chr4:3074876-3074933",
	          "chr4:3074934-3074939",
            "chr4:3074940-3074945",
            "chr4:3074946-3074966",
	          "chr4:3074967-3074972"
        ],
        "VariantId": [
            "HTT",
	          "HTT_CAACAG",
	          "HTT_CCGCCA",
            "HTT_CCG",
	          "HTT_CCT"
        ],
        "VariantType": [
            "Repeat",
            "Repeat",
            "Repeat",
	          "Repeat",
            "Repeat"
        ]
    }
```

The RC workflow is currently divided into 5 steps:
1. RCv0 (find_regex_genotypes_parallel_vcf.py). This is a python script that aims to accurately determine the two most likely genotypes for the HTT repeat structure (Q2-P3) in each genome. First, it scans a BAM file produced by EH to detect the reads that span Q2- P1-P2- P3 (i.e. indices 1-6, according to EH). Then, it uses regular expressions to summarise the HTT structure in each read. It then writes the two most likely HTT structures (by taking the ones supported by the highest number of reads) and the total number of reads supporting each.  When RCv0 finds one repeat structure only, it will output that one structure and the total number of reads that support it.   
```
python find_regex_genotypes_parallel_vcf.py --bam_files BAM_files_list --span 1-6 --count CAA CAACAG "(CAACAG)x2" CCGCCA CCG CCT --output out_1-6_CAG_not_counted.tsv
```
2. V0_postprocessing (clean_data_count_alleles.R). This is a quality control step written in R. It takes the output produced by RCv0 and removes all genotypes supported by 1 read only, or where there the ratio between the total number of reads that support each genotype is > 1:4. 
```
Rscript clean_data_count_alleles.R out.tsv out_cleaned.tsv # cleans data based on read counts
```
3. RCv2 (find_regex_genotypes_parallel_vcf.py). 

### homework - finish this paragraph. "~95% of the genomes have at least one allele with up to 25 CAG repeats. Given that whole genome sequencing (WGS) data uses ~150bp reads..."

This is a python script that aims to determine ONE accurate HTT repeat structure, phased with Q1 from each genome. Like RCv0,  it scans the BAM file produced by EH to detect the reads that span Q1-Q2- P1-P2- P3 (i.e. indices 0-6, according to EH). Then, it uses regular expressions to summarise the HTT structure in each read.   It then writes the two most likely HTT structures (by taking the ones supported by the highest number of reads) and the total number of reads supporting each.  When RCv2 finds one repeat structure only, it will output that one structure and the total number of reads that support it.   
```
python find_regex_genotypes_parallel_vcf.py --bam_files test_files_backups --output out_0-6_CAG-CCT.tsv
```
4. V2_postprocessing (format_RC_data.R). This is an R script that indexes the CAG sizes predicted by EH and those predicted by RCv2, putting the smaller Q1 allele as A1 and the larger allele as A2. It outputs a data-frame with Q1 sizes predicted  by EH, RCv2 as well as the repeat structures predicted by RCv2.
```
Rscript format_RC_data.R out.tsv out_fmtd.tsv # formats CAG alleles in interpretable way
```

5. Phasing (phasing_AD.R). This is an R script that executes the following steps:
	- Matching the EH CAG allele size with the RCv2 CAG allele size giving priority to the smaller allele when both match.
	- Assigning the phased sequence structure to the CAG size using results from RCv2
	- Extracting the second CAG allele size i.e., the unmatched EH CAG genotype
	- Extracting the second unphased sequence structure from RCv0 and phasing it to the unmatched second CAG allele size 



#In summary, we produced a workflow that combines the CAG size predictions by EH together with a custom-made set of scripts to accurately genotype and phase the HTT repeat stuctures found in a normal population (i.e. <35 CAG repeat)
#HTT  = CAG1-SS1 / CAG2-SS2
#Where:
#CAG1 == EH1 == RCv2 - SS1 == RCv2 
#CAG2 == EH2 - SS2 == RCv0 



