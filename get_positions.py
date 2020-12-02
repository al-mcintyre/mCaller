#!/usr/bin/env python
import sys
import numpy as np
from Bio import SeqIO
from scipy import stats
from scipy.cluster.hierarchy import fcluster,linkage
import scipy.spatial.distance as ssd
from extract_contexts import revcomp
from plotlib import plot_w_labels
import os
import pandas as pd

def make_pos_set(pos_list):
    pos_set = set()
    with open(pos_list,'r') as fi:
        for line in fi:
            if len(line) > 3:
                pos_set.add(tuple(line.strip().split('\t')[:4]))
    return pos_set

def aggregate_by_pos(meth_fi,aggfi,depth_thresh,mod_thresh,pos_list,control,verbose_results,gff,ref,plot,plotdir,plotsummary):
    pos_dict = {}
    if verbose_results:
        pos_dict_verbose = {}

    if pos_list:
        pos_set = make_pos_set(pos_list)

    values_dict = {}
    for line in open(meth_fi,'r'):
        #try:
            #print line
            try:
                csome,read,pos,context,values,strand,label,prob = tuple(line.split('\t'))
            except: #for backwards compatibility; does not work with verbose results
                csome,read,pos,context,values,strand,label = tuple(line.split('\t'))
            nextpos = str(int(pos)+1)
            if pos_list and (csome,pos,nextpos,strand) not in pos_set:
                continue
            if (csome,pos,nextpos,context,strand) not in pos_dict:
                pos_dict[(csome,pos,nextpos,context,strand)] = []
                values_dict[(csome,pos,nextpos,context,strand)] = []
                if verbose_results:
                    pos_dict_verbose[(csome,pos,nextpos,context,strand)] = []
            if (pos_list and (csome,pos,nextpos,strand) in pos_set) or (not pos_list and plot):
                values_dict[(csome,pos,nextpos,context,strand)].append([float(v) for v in values.split(',')][:-1])
            if label[0] == 'm':
                pos_dict[(csome,pos,nextpos,context,strand)].append(1)
            else:
                pos_dict[(csome,pos,nextpos,context,strand)].append(0)

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Produce bed file of methylated positions based on mCaller output')
    parser.add_argument('-d','--min_read_depth',type=int,required=False,help='minimum coverage of position to determine methylation (default = 15)',default=15)
    parser.add_argument('-t','--mod_threshold',type=float,required=False,help='minimum %% of observations at a position to include in report (default = 0.5)',default=0.5)
    parser.add_argument('-f','--mCaller_file',type=str,required=True,help='the output file from mCaller to summarize')
    parser.add_argument('-p','--positions',type=str,required=False,help='~bed file of positions for which to calculate % methylated (chromosome,start,end,strand); ignores other thresholds')
    parser.add_argument('--control',action='store_true',required=False,help='take unmethylated positions as a control for motif detection')
    parser.add_argument('--gff',action='store_true',required=False,help='output PacBio-style gff instead of bed ("identificationQv" score will be average probability of methylation)')
    parser.add_argument('--ref',type=str,required=False,help='use reference fasta to output longer contexts surrounding a base, from -20 to +20')
    parser.add_argument('--plot',action='store_true',required=False,help='plot currents deviations at the positions included (not recommended for many positions)')
    parser.add_argument('--plotsummary',action='store_true',required=False,help='plot currents deviations summarized across the positions included')
    parser.add_argument('--plotdir',type=str,required=False,default='mCaller_position_plots',help='output directory for plots, default=mCaller_position_plots')
    parser.add_argument('--vo',action='store_true',required=False,help='verbose output including probabilities for each position')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.3'
        sys.exit(0)

    assert os.path.isfile(args.mCaller_file), 'file not found at '+args.mCaller_file
    if args.positions:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.positions.summary'
    elif not args.control:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.summary'
    else:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.control.summary'
    if args.gff:
        output_file = output_file+'.gff'
    else:
        output_file = output_file+'.bed'
    if not os.path.isdir(args.plotdir):
        os.mkdir(args.plotdir)

    print args.mCaller_file

    aggregate_by_pos(args.mCaller_file,output_file,args.min_read_depth,args.mod_threshold,args.positions,args.control,args.vo,args.gff,args.ref,args.plot,args.plotdir,args.plotsummary)

if __name__ == "__main__":
    main()
