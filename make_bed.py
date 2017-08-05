#!/usr/bin/env python
import sys
import numpy as np
import os

def make_pos_set(pos_list):
    pos_set = set()
    with open(pos_list,'r') as fi:
        for line in fi:
            if len(line) > 3:
                pos_set.add(tuple(line.split('\n')[0].split('\t')))
    return pos_set

def check_thresh(locus_list,mod_thresh,depth_thresh,control):
    if len(locus_list) >= depth_thresh:
        if not control and np.mean(locus_list) >= mod_thresh:
            return True
        elif control and np.mean(locus_list) < mod_thresh:
            return True
        else:
            return False

def aggregate_by_pos(meth_fi,aggfi,depth_thresh,mod_thresh,pos_list,control):
    pos_dict = {}

    if pos_list:
        pos_set = make_pos_set(pos_list)
    #print pos_set

    for line in open(meth_fi,'r'):
        try:
            try:
                csome,read,pos,context,values,strand,label,prob = tuple(line.split('\t'))
            except:
                csome,read,pos,context,values,strand,label = tuple(line.split('\t'))
            nextpos = str(int(pos)+1)
            if (csome,pos,nextpos,context,strand) not in pos_dict:
                pos_dict[(csome,pos,nextpos,context,strand)] = []
            if label[0] == 'm':
                pos_dict[(csome,pos,nextpos,context,strand)].append(1)
            else:
                pos_dict[(csome,pos,nextpos,context,strand)].append(0)
        except:
            pass

        #print pos_dict
    count = 0
    outfi = open(aggfi,'w')
    for locus in pos_dict:
        if not pos_list:
            if check_thresh(pos_dict[locus],mod_thresh,depth_thresh,control):# np.mean(pos_dict[locus]) > mod_thresh and len(pos_dict[locus]) >= depth_thresh:
                count+=1
                outfi.write('\t'.join(list(locus)[:-1]+[str(np.mean(pos_dict[locus]))]+[locus[-1]]+[str(len(pos_dict[locus]))])+'\n')
        else:
            if (locus[0],locus[1],locus[2],locus[4]) in pos_set and 'A' not in set(locus[4]): #TODO: fix main script for As 
                #print locus[1],locus[4]
                outfi.write('\t'.join(list(locus)[:-1]+[str(np.mean(pos_dict[locus]))]+[locus[-1]]+[str(len(pos_dict[locus]))])+'\n')#+[','.join([str(x) for x in pos_dict[locus]])]+'\n')
    if not pos_list:
        print count, 'methylated loci found with min depth', depth_thresh, 'reads'

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Produce bed file of methylated positions based on mCaller output')
    parser.add_argument('-d','--min_read_depth',type=int,required=False,help='minimum coverage of position to determine methylation (default = 15)',default=15)
    parser.add_argument('-t','--mod_threshold',type=float,required=False,help='minimum %% of observations at a position to include in report (default = 0.5)',default=0.5)
    parser.add_argument('-f','--mCaller_file',type=str,required=True,help='the output file from mCaller to summarize')
    parser.add_argument('-p','--positions',type=str,required=False,help='~bed file of positions for which to calculate % methylated (chromosome,start,end,strand)')
    parser.add_argument('--control',type=str,required=False,help='take unmethylated positions as a control for motif detection',default=False)
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.1'
        sys.exit(0)

    assert os.path.isfile(args.mCaller_file), 'file not found at '+args.mCaller_file
    if args.positions:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.positions.summary.bed'
    elif not args.control:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.summary.bed'
    else:
        output_file = args.mCaller_file.split('.')[0]+'.methylation.control.summary.bed'

    print args.mCaller_file

    aggregate_by_pos(args.mCaller_file,output_file,args.min_read_depth,args.mod_threshold,args.positions,args.control)

if __name__ == "__main__":
    main()

