import vcf
import pandas as pd
import sys
import argparse
import os

parser = argparse.ArgumentParser(description='A script to count motif occurences in BAM file')
parser.add_argument('--input', type=str,default=None, help='Output data')
parser.add_argument('--output', type=str,default=None, help='File containing the list of BAM files to be analysed')
parser.add_argument('--gene', type=str,default=None, help='BAM file to be analysed')
parser.add_argument('--motifs',nargs='*', type=str,default=None, help='BAM file to be analysed')
parser.add_argument('--input_dir', type=str,default=None, help='BAM file to be analysed')
parser.add_argument('--output_dir', type=str,default=None, help='BAM file to be analysed')

args = parser.parse_args()

file=args.input
df=pd.read_csv(file,sep='\t')

files=df['bam_file'].str.replace('_realigned.bam','.vcf')

print(args.motifs)
gts=[]
if args.input_dir is not None:
    input_dir=args.input_dir
else:
    input_dir=''

motif_dict={}
for m in args.motifs:
    motif_dict['EH_'+m]=[]

for n in range(len(files)):
    v=files[n]
    print(v)
    v=os.path.join(input_dir,v)
    if n%1000==0:
        print(n/len(files))
    try:
        vcf_reader = vcf.Reader(open(v, 'r'))
    except:
        for m in args.motifs:
            motif_dict['EH_'+m].append('NA')
    try:
        for i in vcf_reader:
            selected,selected_info=i,i.INFO
            if selected_info['RU'] in args.motifs and args.gene in selected_info['VARID']:
                hash_fields=dict(i.INFO)
                hash_fields.update(dict(zip(i.samples[0].data._fields,i.samples[0].data)))
                sizes=hash_fields['REPCN']
                sizes=sizes.replace(', ','/')
                motif=selected_info['RU']
                motif_dict['EH_'+motif].append(sizes)
    except:
        for m in args.motifs:
            motif_dict['EH_'+m].append('NA')

print(motif_dict)
eh_df=pd.DataFrame(motif_dict)
print(eh_df)
print(df)
df=pd.concat([df,eh_df],axis=1)
print(df)
if args.output_dir is not None:
    output=os.path.join(args.output_dir,args.output)
print(output)
df.to_csv(output,sep='\t',index=False)

