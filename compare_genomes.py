#!/usr/bin/env python
import sys
import numpy as np
import os
from scipy.stats import mannwhitneyu,ranksums,ttest_ind,ks_2samp

def compare_by_position(bed1,bed2,xmfa):
    pos_dict = {}

    for i,bed in enumerate([bed1,bed2]):
        pos_dict[i] = {}
        with open(bed,'r') as fi:
                for line in fi:
                #2  1892198 1892199 TCMMTMTTMMM 0.5 -   16
                    csome,start,end,motif,perc_meth,strand,num_reads,probabilities = tuple(line.split('\t'))
                    pos_dict[i][(csome,start,end,strand)] = ((perc_meth,num_reads),np.asarray([float(p) for p in probabilities.strip().split(',')]))

    for pos in pos_dict[0]:
        if pos in pos_dict[1]:
            try:
                u,pval = mannwhitneyu(pos_dict[0][pos][1],pos_dict[0][pos][1],alternative='two-sided')
            except ValueError:
                u,pval = 'none','identical'
            u2,pval2 = ranksums(pos_dict[0][pos][1],pos_dict[0][pos][1])
            try:
                t,pval3 = ttest_ind(pos_dict[0][pos][1],pos_dict[0][pos][1])
            except:
                t,pval3 = 'none','missing df'
            d,pval4 = ks_2samp(pos_dict[0][pos][1],pos_dict[0][pos][1])
            if pval4 < 0.9:
                print(pos, pos_dict[0][pos][0], pos_dict[1][pos][0], pval, pval2, pval3, pval4)
                #print pos, pos_dict[0][pos][0], pos_dict[1][pos][0], pval, pval2, pval3, pval4



def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Compare methylation between two genomes by probabilities of methylation for aligned positions')
    parser.add_argument('--bed1',type=str,required=True,help='bed file 1 with verbose output from make_bed.py')
    parser.add_argument('--bed2',type=str,required=True,help='bed file 2 with verbose output from make_bed.py')
    parser.add_argument('-g','--genome_alignment',type=str,required=False,help='an xmfa file from mauve (if absent, alignments assumed to be to the same reference genome)')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print('mCallerNP 0.3')
        sys.exit(0)

    assert os.path.isfile(args.bed1), 'file not found at '+args.bed1
    assert os.path.isfile(args.bed2), 'file not found at '+args.bed2

    compare_by_position(args.bed1,args.bed2,args.genome_alignment)

if __name__ == "__main__":
    main()

