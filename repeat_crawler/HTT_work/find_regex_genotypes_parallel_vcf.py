import pysam
import re
import vcf
import sys
import os
import pandas as pd
import numpy as np
import argparse
# create dictionary with a regex for all possible genotypes of each of the HTT locus components
components={
    'Q1':["(CAG){2,}"],
    'Q2':["(?<!CAACAG)CAACAG(?!CAACAG)","(CAACAG){2,}","(CAA)"],
    'P1':["(CCGCCA){1,}"],
    'P2':["(CCG){1,}"],
    'P3':["(CCT){1,}"]
}

parser = argparse.ArgumentParser(description='A script to count motif occurences in BAM file')
parser.add_argument('--output', type=str,default=None, help='Output data')
parser.add_argument('--bam_files', type=str,default=None, help='File containing the list of BAM files to be analysed')
parser.add_argument('--bam_file', type=str,default=None, help='BAM file to be analysed')
parser.add_argument('--vcf_informed', type=str,default=None, help='Record vcf genotypes')
parser.add_argument('--n_mismatches', type=int,default=0, help='Number of mismatches allowed in a read')
parser.add_argument('--span', type=str,default='0-6', help='filter reads that span components- default is 0-6')
parser.add_argument('--must_contain',nargs='*',default=['CAG','CCT'], help='Repeat structure must contain')
parser.add_argument('--count',nargs='*',default=['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT'], help='Components to consider when counting read type frequency')

args = parser.parse_args()



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

def allow_mismatches(current_motif,next_motif,read_seq,component_stretch_counts):
    from_here,to_here=len(next_motif),len(next_motif)+len(current_motif)
    if read_seq[:from_here]!=next_motif and read_seq[from_here:to_here]==current_motif:
        query,read_seq_new,component_stretch_counts=query_count_and_cut("(CAG){2,}",read_seq[from_here:],{},find_start=True)
        if query and len(read_seq_new)<len(read_seq):
            difference=len(read_seq)-len(read_seq_new)
            extra_reps=difference/len(current_motif)
        else:
            query,extra_reps,difference=None,None,None
    else:
        query,extra_reps,difference=None,None,None
    return query,extra_reps,read_seq

# This function takes the read string to be searched and spits 
# out a dictionary counting the stretch of each component and the most likely gt 
def construct_regex_and_count_stretches(read_seq,components,n_mismatches_allowed=0):
    output_regex=[]
    component_stretch_counts={}
    query,read_seq,component_stretch_counts=query_count_and_cut("(CAG){2,}",read_seq,component_stretch_counts,find_start=True)
    for mis in range(n_mismatches_allowed):
        query_new,extra_reps,read_seq_new=allow_mismatches('CAG','CAACAG',read_seq,component_stretch_counts)
        if query_new and extra_reps:
            # print(read_seq)
            new_cut=len('CAG')*int(extra_reps)
            read_seq=read_seq_new[new_cut:]
            # print(read_seq)
            # print(component_stretch_counts['CAG'])
            # print(extra_reps)
            component_stretch_counts['CAG'][0]=int(component_stretch_counts['CAG'][0])+extra_reps
    
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

def extract_vcf_info(v):
    vcf_reader = vcf.Reader(open(v, 'r'))
    for i in vcf_reader:
        selected,selected_info=i,i.INFO
        if i.INFO['RU']=='CAG':
            hash_fields=dict(i.INFO)
            hash_fields.update(dict(zip(i.samples[0].data._fields,i.samples[0].data)))
            gt=hash_fields['REPCN'].replace('/',',')
    return gt


def para(bam_files,vcf_informed,from_component,to_component,n_mismatches_allowed=0):
    counts,EH_rep_nos=[],[]
    for bam_file in bam_files:
        result=process_bam_file(bam_file,n_mismatches_allowed,from_component,to_component)
        if vcf_informed is not None:
            eh_reps=extract_vcf_info(bam_file.replace('_realigned.bam','.vcf',))
            EH_rep_nos.append(eh_reps)
        counts.append(result)
        # pool.close()
        # pool.join()
    # return [stretch_count.get() for stretch_count in counts]
    return counts,EH_rep_nos

# loop through reads and for each one- create a dictionary containing 
# the stretch count of each component of the HTT locus and then propose the most likely gt 
# each dictionary will be appended to a list

def process_bam_file(input_file,n_mismatches_allowed,from_component,to_component):
    genotypes=[]
    genotype_count={}
    component_stretch_counts_list=[]
    input_file = pysam.AlignmentFile(input_file, "rb")
    for i,read in enumerate(input_file):
        read_seq=read.get_forward_sequence()
        if read.is_reverse:
            read_seq=revc(read_seq)
        if from_component!='None' and re.search(str(from_component)+"\[.*\]"+str(to_component),read.to_dict()['tags'][0].split(',')[2]): # ensure Q1-P3
            # read_seq=modify_string(read_seq,read.to_dict()['qual'])
            component_stretch_counts,genotype=construct_regex_and_count_stretches(read_seq,components,n_mismatches_allowed)
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
        else:
            component_stretch_counts,genotype=construct_regex_and_count_stretches(read_seq,components,n_mismatches_allowed)
            if genotype is not None:
                component_stretch_counts_list.append(component_stretch_counts)
                if genotype not in genotype_count:
                    genotype_count[genotype]=1
                else:
                    genotype_count[genotype]+=1
    input_file.close()
    return component_stretch_counts_list


def format_header(component_stretch_counts_list,required_keys=['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT','genotype']):
    # ensure all dictionaries have all header info
    for dictionary in component_stretch_counts_list:
        for key in required_keys:
            if key not in dictionary:
                dictionary[key] = [0]
    df=pd.concat([pd.DataFrame(dictionary) for dictionary in component_stretch_counts_list], ignore_index=True)
    df=df[required_keys]
    return df

def count_and_write(bam_file,component_stretch_counts_list,count_list,print_df=False):
    # put all dictionaries containing stretch counts in a dataframe each row represents a read
    # count the number of reads supporting a given GT with a certain composition
    if args.vcf_informed is not None:
        eh_reps=EH_rep_nos[i].split(',')
        df=format_header(component_stretch_counts_list,required_keys=['CAG','CAA','CAACAG','(CAACAG)x2','CCGCCA','CCG','CCT','genotype'])
        # df=df[df['genotype'].str.startswith('CAG')]
        for must in args.must_contain:
            df=df[df['genotype'].str.contains(must)]
        df=df.loc[df['CAG'].astype(int).isin([int(j) for j in eh_reps])]
        df=df.groupby(count_list+['genotype']).size().reset_index(name='Count').sort_values(by='Count',ascending=False).reset_index(drop=True)
    else:
        df=format_header(component_stretch_counts_list)
        # df=df[df['genotype'].str.startswith('CAG')]
        for must in args.must_contain:
            df=df[df['genotype'].str.contains(must)]
        df=df[df['CAG']>=6].reset_index(drop=True)
        df=df.groupby(count_list+['genotype']).size().reset_index(name='Count').sort_values(by='Count',ascending=False).reset_index(drop=True)
        # if len(df['Count'].drop_duplicates())==1:
        #     df=df.sort_values(by='CAG',ascending=False).reset_index(drop=True)
    if print_df:
        print(df)
    df['bam_file']=bam_file
    df=df[['bam_file','genotype','Count']+count_list]
    df=df.rename(columns={'genotype':'first_second_gt','Count':'read_counts'})
    output_header_str="\t".join(df.columns.to_list())+"\n"
    data_output_str=[]
    data_output_str.append(df.loc[0,'bam_file'])
    for h in ['first_second_gt','read_counts']:
        if df.shape[0]>1:
            add=str(df.loc[0,h])+'/'+str(df.loc[1,h])
        else:
            add=str(df.loc[0,h])
        if h=='first_second_gt':
            add='('+add+')'
            add=add.replace('/',')/(')
        data_output_str.append(add)
    for h in count_list:
        if df.shape[0]>1:
            add=str(df.loc[0,h])+'/'+str(df.loc[1,h])
        else:
            add=str(df.loc[0,h])
        data_output_str.append(add)
    data_output_str="\t".join(data_output_str)+'\n'
    # write 2 most common allele types and stretch counts to file 
    # for stretch counts there is a separate column for CAACAG and (CAACAG)x2 
    if not os.path.exists(args.output):
        with open(args.output,'w') as out:
            out.write(output_header_str)
            out.write(data_output_str)
    else:
        with open(args.output,'a') as out:
            out.write(data_output_str)

if args.span=='None':
    from_component,to_component='None','None'
else:
    from_component,to_component=args.span.split('-')

count_list=args.count

if args.bam_files is not None:
    bam_files=pd.read_csv(args.bam_files,header=None)[0].to_list() # input bam file
    lists,EH_rep_nos=para(bam_files,args.vcf_informed,from_component,to_component,args.n_mismatches)
    for i in range(len(lists)):
        bam_file=bam_files[i]
        component_stretch_counts_list=lists[i]
        count_and_write(bam_file,component_stretch_counts_list,count_list)
else:
    bam_file=args.bam_file
    component_stretch_counts_list=process_bam_file(bam_file,args.n_mismatches,from_component,to_component)
    count_and_write(bam_file,component_stretch_counts_list,count_list,True)


