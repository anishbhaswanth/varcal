#!/usr/bin/env python
# coding: utf-8

"""
Author: Anish Chakka
Email: anishbhaswanth@gmail.com
Last updated: 07/25/21


This code takes bam and fasta as input.
Variant calls are written in tabular txt format.
"""

# import modules
import os
import pysam
import pandas as pd
import re
from itertools import dropwhile, chain
import operator
import numpy

# function to convert bam to sam file
def convert_to_sam(file):
   
    # try reading the input bam
    try:
        bamfile = pysam.AlignmentFile(file, 'rb')

        # get seq names from bam header
        bam_refs = bamfile.header['SQ']
        bam_ref_chr = []
        for i in range(len(bam_refs)):
            bam_ref_chr.append(bam_refs[i]['SN'])
    
    except FileNotFoundError as fnf_error:
        print(fnf_error)


    # create sam file
    sam_outfile = os.path.splitext(file)[0]+'.sam' # replace .bam with .sam
    samfile = pysam.AlignmentFile(sam_outfile, "w", template=bamfile)

    # convert it to sam and write to a file
    for i in bamfile:
        samfile.write(i)

    return(sam_outfile, bam_ref_chr)


# function to read and store fasta file
def read_fasta(ref_file, refs_chr):
    ref_fasta = {}
    chr_name = None

    try:
        with open(ref_file) as f:
            for line in f:
                if line.startswith('>'):
                    chr_name = line[1:].rstrip()
                    ref_fasta[chr_name] = []
                else:
                    ref_fasta[chr_name].append(line.rstrip())

            for chr_name, chr_seq in ref_fasta.items():
                ref_fasta[chr_name] = ''.join(chr_seq)

    except FileNotFoundError as fnf_error:
        print(fnf_error)

    # check if BAM ref matches with fasta ref
    if(set(refs_chr).issubset(set(ref_fasta.keys()))):
        pass

    else:
        print("Reference sequences DOES NOT match with BAM reference!\n\n")
        print("Reference fasta has: ", list(ref_fasta.keys()), "\n\n")
        print("While, BAM file has: ", refs_chr)
        quit()

    return(ref_fasta)


# check for SAM header
def is_header(s):
    """ function to check for 
        SAM file's header
    """
    # return true if a line starts with @
    return s.startswith('@')


# process CIGAR strings in SAM file
def process_cigar(cigar_string, start_pos):
    
    # store Del, Ins, M/Mis and Soft clippings
    dims = []
    for i, cigar in re.findall('(\d+)([MDIS])', cigar_string):
        if cigar == 'I':
            dims.append(list(int(i) * 'I'))
        if cigar == 'M':
            dims.append(list(range(start_pos, start_pos+int(i))))
            start_pos = start_pos + int(i)
        elif cigar in 'DS':
            start_pos = start_pos + int(i)

    return(list(chain.from_iterable(dims)))

# read sam file to initialize variants
def parse_sam(sam_file):
    # list for storing variants
    var_list = []
    try:
        # read sam file and filter out multiple mapped reads
        with open(sam_file, 'r') as file:
            # loop through each line
             for line in dropwhile(is_header, file):
                    # parse only uniquely mapped regions
                    if re.search("XA:Z:", line) is None:
                        # split columns by \t
                        cols = line.split('\t')
                
                        # parse cigar string to get the right seq length
                        cigar_pos = process_cigar(cols[5], int(cols[3]))
                        
                        # loop through each nucleotide of the seq
                        i = 0
                        for pos in cigar_pos:
                            # list to store variant info
                            var = []
                            # skip inserstions in seq
                            if pos == "I":
                                i = i + 1
                                continue
                            
                            # add chr 
                            var.append(cols[2])

                            # add ALT pos
                            var.append(pos)

                            # add ALT
                            var.append(cols[9][i])

                            # add BASE QUAL (converting them to phred33)
                            var.append(ord(cols[10][i]) - 33)

                            # add MAP QUAL
                            var.append(int(cols[4]))

                             # strand (+ first, - second)
                            pos_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N':0}
                            neg_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0, 'N':0}
                            
                            if(int(cols[1]) & 16):
                                var.extend(list(pos_counts.values()))
                                neg_counts[cols[9][i]] = 1
                                var.extend(list(neg_counts.values()))
                            else:
                                pos_counts[cols[9][i]] = 1
                                var.extend(list(pos_counts.values()))
                                var.extend(list(neg_counts.values()))

                            # # if reverse strand
                            # if(int(cols[1]) & 16):
                            #     var.extend([0,1])
                            # else:
                            #     var.extend([1,0])

                            # variant position in read
                            var.append(i+1)

                            # coverage
                            var.append(1)
                            
                            i = i+1
                            var_list.append(var)

    except FileNotFoundError as fnf_error:
        print(fnf_error)

    return(var_list)


# get base counts for ALT
def get_base_counts(alt_seq, ref_nucl):
    counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    for nucl in alt_seq:
        if not nucl == 'N':
            counts[nucl] += 1
    
    max_nucl_counts = max(counts.items(), key=operator.itemgetter(1))
    if max_nucl_counts[0] == ref_nucl:
        # delete this base to look for next best one
        del(counts[max_nucl_counts[0]])
        new_max_nucl_counts = max(counts.items(), key=operator.itemgetter(1))
        if new_max_nucl_counts[1] == 0:
            return (ref_nucl, 0)
        else:
            return (new_max_nucl_counts)
    else:
        return max_nucl_counts


# format variants list
def format_list_of_variants(var_list, ref_fasta):
    # rename columns
    var_df = pd.DataFrame(var_list, columns = ['Chr', 
                                           'POS', 
                                           'ALT', 
                                           'BQ', 
                                           'MQ', 
                                           'Pos_AStrand',
                                           'Pos_TStrand',
                                           'Pos_GStrand',
                                           'Pos_CStrand',
                                           'Pos_NStrand',
                                           'Neg_AStrand',
                                           'Neg_TStrand',
                                           'Neg_GStrand',
                                           'Neg_CStrand',
                                           'Neg_NStrand',
                                           'PiR', 
                                           'Coverage'])

    # aggregate
    var_df_mult = var_df.groupby(['Chr', 'POS']).agg({'ALT': ['sum'], 
                                                  'BQ': ['mean'], 
                                                  'MQ': ['mean'],
                                                  'Pos_AStrand': ['sum'],
                                                  'Pos_TStrand': ['sum'],
                                                  'Pos_GStrand': ['sum'],
                                                  'Pos_CStrand': ['sum'],
                                                  'Pos_NStrand': ['sum'],
                                                  'Neg_AStrand': ['sum'],
                                                  'Neg_TStrand': ['sum'],
                                                  'Neg_GStrand': ['sum'],
                                                  'Neg_CStrand': ['sum'],
                                                  'Neg_NStrand': ['sum'],
                                                  'PiR': ['mean'],
                                                  'Coverage': ['sum']})

    # rename columns
    var_df_mult.columns = ['ALT_Bases', 'Avg_BQ', 'Avg_MQ', 
                       'Pos_AStrand', 'Pos_TStrand', 'Pos_GStrand', 'Pos_CStrand', 'Pos_NStrand',
                       'Neg_AStrand', 'Neg_TStrand', 'Neg_GStrand', 'Neg_CStrand', 'Neg_NStrand',
                       'Avg_PiR', 'Coverage']
    var_df_mult = var_df_mult.reset_index()

    # add REF base to the table
    ref_pos = []
    for i, j in zip(var_df_mult['Chr'], var_df_mult['POS']):
        ref_pos.append(ref_fasta[i][j-1].upper())
    
    var_df_mult['REF'] = ref_pos

    var_df_mult = var_df_mult.reset_index(drop=True)
    var_df_mult.head()

    # find best ALT base
    # calculate ALT DP
    # calculate REF DP

    i = 0
    ref_count = []
    alt_nucl = []
    alt_count = []

    while i < len(var_df_mult.index):
        ref_count.append(var_df_mult.loc[i,'ALT_Bases'].count(var_df_mult.loc[i,'REF']))
        
        # find alt base with max count
        max_nucl_counts = get_base_counts(var_df_mult.loc[i,'ALT_Bases'], var_df_mult.loc[i,'REF'])
        
        alt_nucl.append(max_nucl_counts[0])
        alt_count.append(max_nucl_counts[1])
        
        i = i+1

    var_df_mult['ALT'] = alt_nucl
    var_df_mult['ALT_DP'] = alt_count
    var_df_mult['REF_DP'] = ref_count
    var_df_mult['AF'] = var_df_mult['ALT_DP']/var_df_mult['Coverage']

    # add number of total pos and neg strands
    var_df_mult['Total_Pos_Strand'] = var_df_mult.iloc[:, 5:10].sum(axis=1)
    var_df_mult['Total_Neg_Strand'] = var_df_mult.iloc[:, 10:15].sum(axis=1)


    # calculate entropy
    numpy.seterr(divide = 'ignore') 
    var_df_mult['Entropy'] = -(((var_df_mult['Total_Pos_Strand']/var_df_mult['Coverage']) * 
    numpy.log2(var_df_mult['Total_Pos_Strand']/var_df_mult['Coverage'])) + ((var_df_mult['Total_Neg_Strand']/var_df_mult['Coverage']) * 
    numpy.log2(var_df_mult['Total_Neg_Strand']/var_df_mult['Coverage'])))



    return(var_df_mult)

def filter_variants(var_df_mult, min_cov = 5, min_bq = 15, min_adp = 2, min_af = 0.01):
    var_df_mult = var_df_mult[(var_df_mult['Coverage'] >= min_cov) & (var_df_mult['Avg_BQ'] >= min_bq)]
    var_df_mult = var_df_mult[(var_df_mult['ALT_DP'] >= min_adp) & (var_df_mult['AF'] >= min_af)]

    return(var_df_mult)

# my main function
def main(options):
    # convert the file sam
    sam_file, refs_chr = convert_to_sam(options.bam)

     # store fasta file
    my_ref = read_fasta(options.ref, refs_chr)

    # parse the sam file to get variants
    list_of_variants = parse_sam(sam_file)

    # properly format the variant list
    var_df_mult = format_list_of_variants(list_of_variants, my_ref)
    
    # filter variants
    var_df_mult_filt = filter_variants(var_df_mult)

    # output filtered tab delimited variants
    var_df_mult_output = var_df_mult_filt[['Chr', 'POS', 'REF', 'ALT', 'Avg_BQ', 'Avg_MQ', 'Coverage', 
                                    'Avg_PiR', 'ALT_DP', 'AF' , 'Entropy', 'Total_Pos_Strand', 'Total_Neg_Strand',
                                  'Pos_AStrand', 'Pos_TStrand', 'Pos_GStrand', 'Pos_CStrand', 'Pos_NStrand',
                                  'Neg_AStrand', 'Neg_TStrand', 'Neg_GStrand', 'Neg_CStrand', 'Neg_NStrand',]]

                                  
    var_df_mult_output.to_csv(options.out, sep = "\t", index = False, na_rep='NaN')


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Generate variant calls from a BAM file using ref fasta')
    parser.add_argument('-b', '--bam', help='bam file', required=True)
    parser.add_argument('-r', '--ref', help='reference fasta', required=True)
    parser.add_argument('-o', '--out', help='file name to output variant calls', required=True)   
    options = parser.parse_args()
    main(options)
