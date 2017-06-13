#!/usr/bin/env python
import sys
import numpy as np
import os

def aggregate_by_pos(meth_fi,aggfi,depth_thresh,mod_thresh):
    pos_dict = {}

    for line in open(meth_fi,'r').readlines():
        try:
            csome,read,pos,context,values,strand,label = tuple(line.split('\t')) #SHOULD HAVE CSOME
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
        if np.mean(pos_dict[locus]) > mod_thresh and len(pos_dict[locus]) >= depth_thresh:
            count+=1
            outfi.write('\t'.join(list(locus)[:-1]+[str(np.mean(pos_dict[locus]))]+[locus[-1]]+[str(len(pos_dict[locus]))])+'\n')
    print count, 'methylated loci found with min depth', depth_thresh, 'reads'

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Produce bed file of methylated positions based on mCaller output')
    parser.add_argument('-d','--min_read_depth',type=int,required=False,help='minimum coverage of position to determine methylation (default = 15)',default=15)
    parser.add_argument('-t','--mod_threshold',type=float,required=False,help='minimum %% of observations at a position to include in report (default = 0.5)',default=0.5)
    parser.add_argument('-f','--mCaller_file',type=str,required=True,help='the output file from mCaller to summarize')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.1'
        sys.exit(0)

    assert os.path.isfile(args.mCaller_file), 'file not found at '+args.mCaller_file
    output_file = args.mCaller_file.split('.')[0]+'.methylation.summary.bed'
    print args.mCaller_file

    aggregate_by_pos(args.mCaller_file,output_file,args.min_read_depth,args.mod_threshold)

if __name__ == "__main__":
    main()

