import pysam
import re
import sys
import os
import pandas as pd
import numpy as np

# create dictionary with a regex for all possible genotypes of each of the HTT locus components
components={
    'Q1':["(CAG){2,}"],
    'Q2':["(?<!CAACAG)CAACAG(?!CAACAG)","(CAACAG){2,}","(CAA)"],
    'P1':["(CCGCCA){1,}"],
    'P2':["(CCG){1,}"],
    'P3':["(CCT){1,}"]
}

bam_file=sys.argv[1] # input bam file

# This function searches for a component pattern in a spanning read- then when its found it- 
# reports the number repeats and cuts the read string such that the pattern 
# that comes next will be at the start of the read string- which allows one to build the regex 
# that most suits the read
# e.g. 
# read_string=CAGCAGCAG CAACAG CCGCCA...
# regex='' -> search CAG -> regex='CAG|' -> cut
# read_string=CAACAG CCGCCA...
# regex='' -> search CAACAG -> regex='CAG|CAACAG|' -> cut
# read_string=CCGCCA...
# regex='' -> search CCGCCA -> regex='CAG|CAACAG|CCGCCA' -> cut...

def query_count_and_cut(regex,read_seq,component_stretch_counts,find_start=False):
    if find_start: # find initial CAG stretch or insist that string must start with the next component
        query=re.search(regex,read_seq)
    else:
        query=re.match(regex,read_seq)
    if query:
        span=query.span()
        length=span[1]-span[0]
        if regex=='(CAACAG){2,}':
            regex_stripped_down='(CAACAG)x2'
            motif_length=len('CAACAG')
            rep_no=length/motif_length
        elif regex=="(?<!CAACAG)CAACAG(?!CAACAG)":
            regex_stripped_down='CAACAG'
            motif_length=len('CAACAG')
            rep_no=length/motif_length
        else:
            regex_stripped_down=re.sub(r'[^A-Z]', '', regex)
            motif_length=len(regex_stripped_down)
            rep_no=length/motif_length
        if regex_stripped_down not in component_stretch_counts:
            component_stretch_counts[regex_stripped_down]=[]
        component_stretch_counts[regex_stripped_down].append(rep_no)
        read_seq=read_seq[span[1]:]
    return query,read_seq,component_stretch_counts


# This function takes the read string to be searched and spits 
# out a dictionary counting the stretch of each component and the most likely gt 
def construct_regex_and_count_stretches(read_seq,components):
    output_regex=[]
    component_stretch_counts={}
    query,read_seq,component_stretch_counts=query_count_and_cut("(CAG){2,}",read_seq,component_stretch_counts,find_start=True)
    if query:
        for i in range(1,len(components.keys())):
            component=list(components.keys())[i]
            for regex in components[component]:
                query,read_seq,component_stretch_counts=query_count_and_cut(regex,read_seq,component_stretch_counts)
        output_regex="|".join(component_stretch_counts.keys())
    else:
        output_regex=None
    component_stretch_counts['genotype']=output_regex
    # if output_regex=='CAG|CCGCCA|CCG|CCT':
    #     print(read.get_forward_sequence())
    #     print(revc(read.get_forward_sequence()))
    #     print(component_stretch_counts)
    #     print(read.to_dict()['qual'])
    return component_stretch_counts,output_regex

def modify_string(read_seq,qual):
    letters=[]
    for i in range(len(read_seq)):
        if qual[i]=='!':
            letters.append(read_seq[i].lower())
        else:
            letters.append(read_seq[i])
    return ''.join(letters)

def revc(seq):
    seq = seq.replace("A", "t").replace("C", "g").replace("T", "a").replace("G", "c")
    seq = seq.upper()
    seq=seq[::-1]
    return seq

genotypes=[]
genotype_count={}
component_stretch_counts_list=[]
input_file = pysam.AlignmentFile(bam_file, "rb")
# loop through reads and for each one- create a dictionary containing 
# the stretch count of each component of the HTT locus and then propose the most likely gt 
# each dictionary will be appended to a list
for i,read in enumerate(input_file):
    read_seq=read.get_forward_sequence()
    if read.is_reverse:
        read_seq=revc(read_seq)
    if re.search("1\[.*\]6",read.to_dict()['tags'][0].split(',')[2]): # ensure Q1-P3
        # read_seq=modify_string(read_seq,read.to_dict()['qual'])
        component_stretch_counts,genotype=construct_regex_and_count_stretches(read_seq,components)
        # if component_stretch_counts['genotype']=='CAG|CAA|CCG|CCT':
        #     print(read.get_forward_sequence())
        #     print(revc(read.get_forward_sequence()))
        #     print(component_stretch_counts)
        #     print(read.to_dict()['qual'])
        # if "CAACAGCAACAG" in read_seq:
        #     print(read_seq)
        #     print(component_stretch_counts)
        if genotype is not None:
            component_stretch_counts_list.append(component_stretch_counts)
            if genotype not in genotype_count:
                genotype_count[genotype]=1
            else:
                genotype_count[genotype]+=1

input_file.close()



required_keys=['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT']
# ensure all dictionaries have all header info
for dictionary in component_stretch_counts_list:
    for key in required_keys:
        if key not in dictionary:
            dictionary[key] = [0]

# put all dictionaries containing stretch counts in a dataframe each row represents a read
# count the number of reads supporting a given GT with a certain composition
df=pd.concat([pd.DataFrame(dictionary) for dictionary in component_stretch_counts_list], ignore_index=True)
df=df[['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT','genotype']]
df=df[df['genotype'].str.startswith('CAG')]
df=df[df['genotype'].str.endswith('CCT')]
df=df.groupby(['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT','genotype']).size().reset_index(name='Count').sort_values(by='Count',ascending=False).reset_index(drop=True)
df=df[df['CAG']>=6].reset_index(drop=True)
if df.shape[0]>=10:
    df=df[df["Count"]>1]

# print(df)
if df.shape[0]>1:
    if len(df['genotype'].drop_duplicates())>1:
        df_out=pd.DataFrame()
        for df_sub in df.groupby('genotype'):
            # df_sub=df_sub[1].sort_values(by=['CAG','CCT'],ascending=False).reset_index(drop=True).loc[0,:].T
            add=df_sub[1].iloc[0,:]
            # if add['Count']>1:
            df_out=df_out.append(add,ignore_index=True)
        df=df_out.reset_index(drop=True).sort_values(by='Count',ascending=False).reset_index(drop=True)

df['bam_file']=bam_file
# print(df)
# write 2 most common allele types and stretch counts to file 
# for stretch counts there is a separate column for CAACAG and (CAACAG)x2 

if not os.path.exists(sys.argv[2]):
    with open(sys.argv[2],'w') as out:
        out.write('bam_file\tfirst_second_gt\tread_counts\tCAG\tCAA\tCAACAG\t(CAACAG)x2\tCCGCCA\tCCG\tCCT\n')

if df.shape[0]>1:
    with open(sys.argv[2],'a') as out:
        out.write('{0}\t({1})/({2})\t{3}/{4}\t{5}/{6}\t{7}/{8}\t{9}/{10}\t{11}/{12}\t{13}/{14}\t{15}/{16}\t{17}/{18}\n'.format(df.loc[0,'bam_file'],df.loc[0,'genotype'],df.loc[1,'genotype'],df.loc[0,'Count'],df.loc[1,'Count'],df.loc[0,'CAG'],df.loc[1,'CAG'],df.loc[0,'CAA'],df.loc[1,'CAA'],df.loc[0,'CAACAG'],df.loc[1,'CAACAG'],df.loc[0,'(CAACAG)x2'],df.loc[1,'(CAACAG)x2'],df.loc[0,'CCGCCA'],df.loc[1,'CCGCCA'],df.loc[0,'CCG'],df.loc[1,'CCG'],df.loc[0,'CCT'],df.loc[1,'CCT']))
else:
    with open(sys.argv[2],'a') as out:
        out.write('{0}\t({1})\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n'.format(df.loc[0,'bam_file'],df.loc[0,'genotype'],df.loc[0,'Count'],df.loc[0,'CAG'],df.loc[0,'CAA'],df.loc[0,'CAACAG'],df.loc[0,'(CAACAG)x2'],df.loc[0,'CCGCCA'],df.loc[0,'CCG'],df.loc[0,'CCT']))


