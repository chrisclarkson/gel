# Code that I have used on the Genomics england platform when assessing repeat expansion loci (specifically for short read data- at this time).
## Repeat crawler is the tool that I developed to characterise STR motif constitution- look in `repeat_crawler/` for details.

## A separate project that I work on is to identify novel pathogenic Repeat expansion loci -using tools documented in `denovo_RE_finding/`

## running expansion hunter on GEL:
A generic run of ExpansionHunter on GEL will be done using this bash script (for a given `$gene`, on a set of genomes contatined in a job array file):
```
for gene in gene1 gene2 gene3
do
bsub -q short -P re_gecip_neurology -J "Array[1-20262:100]" -R rusage[mem=48000] -M 48000 -n 3 \
-o "${PWD}/jobs/myJob.out%I" -e "${PWD}/jobs/myJob.err%I" \
"${PWD}/runEH.sh $gene ${PWD}/100k_paths_annotated ${PWD}/EH"
done
```
In the above example- there will be a `gene1.json`,`gene2.json` as per ExpansionHunter format to specify STR parameters etc.
