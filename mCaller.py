#!/usr/bin/env python
#A program to classify bases as methylated or unmethylated based on long-range signals using the output from nanopolish
#Alexa McIntyre, 2016-2018

import sys
assert sys.version_info >= (3,0), 'please use python3'
import pickle 
from collections import defaultdict
import numpy as np
#import time
import os
import glob
import math
import multiprocessing
import warnings
from Bio import SeqIO

from extract_contexts import *
from train_model import train_classifier,pos2label
from read_qual import extract_read_quality
from load_mCaller_data import tsv2matrix

warnings.filterwarnings("ignore", module="sklearn")

def distribute_threads(positions_list,motif,tsvname,read2qual,refname,num_refs,base,mod,nprocs,nvariables,train,modelfile,skip_thresh,qual_thresh,classifier,training_tsv,plot_training):
    """ distributes list of genomic positions across processes then adds resulting signals to matrix if training"""
    outdir='/'.join(tsvname.split('/')[:-1])
    if len(outdir) > 1:
        outdir = outdir + '/'
    print(outdir)
    if not train:
      tsv_output = '.'.join(tsvname.split('.')[:-1])+'.diffs.'+str(nvariables)
      training_pos_dict = None
    else: 
      tsv_output = '.'.join(tsvname.split('.')[:-1])+'.diffs.'+str(nvariables)+'.train'
      if training_tsv:
        signal_mat, context_array = tsv2matrix(training_tsv,base)
      else:
          training_pos_dict = pos2label(positions_list)

    print('%d contigs' % num_refs)
    print('%d threads' %nprocs)
    #sys.exit(0)

    if not training_tsv:
        bytesize = os.path.getsize(tsvname)

    if nprocs > 1 and not training_tsv:
        nprocs_allocated = 0
        procs = []
        out_q = multiprocessing.Queue()
        def worker(out_q,tsvname,fastaname,read2qual,nvariables,skip_thresh,qual_thresh,modelfile,classifier,startline,endline,train,training_pos_dict,base=base,motif=None,positions_list=None): 
            outtup = extract_features(tsvname,fastaname,read2qual,nvariables,skip_thresh,qual_thresh,modelfile,classifier,startline,endline=endline,train=train,pos_label=training_pos_dict,base=base,motif=motif,positions_list=positions_list)
            out_q.put(outtup)

    if nprocs == 1 and not training_tsv:
        if not train:
            extract_features(tsvname,refname,read2qual,nvariables,skip_thresh,qual_thresh,modelfile,classifier,0,endline=bytesize,train=train,pos_label=training_pos_dict,base=base,motif=motif,positions_list=positions_list) 
        else:
            signal_mat, context_array = extract_features(tsvname,refname,read2qual,nvariables,skip_thresh,qual_thresh,modelfile,classifier,0,endline=bytesize,train=train,pos_label=training_pos_dict,base=base,motif=motif,positions_list=positions_list)

    elif nprocs > 1 and not training_tsv:
        chunksize = int(math.ceil(bytesize/float(nprocs))) 

        for i in range(nprocs):
            p = multiprocessing.Process(
                    target=worker, 
                    args=(out_q,tsvname,refname,read2qual,nvariables,skip_thresh,qual_thresh,modelfile,classifier,chunksize*i,chunksize*(i+1),train,training_pos_dict,base,motif,positions_list))
            procs.append(p)
            p.start()

        if train:
            # Collect all results into a signal matrix and an array of labels
            signal_mat = {}
            context_array = {}
            for i,proc in enumerate(procs):
                tmp_signal_mat,tmp_contexts = out_q.get()
                for twobase in tmp_signal_mat:
                    if twobase not in signal_mat:
                        signal_mat[twobase] = {}
                        context_array[twobase] = {}
                    for label in tmp_signal_mat[twobase]:
                        if label not in signal_mat[twobase]:
                            signal_mat[twobase][label] = []
                            context_array[twobase][label] = []
                        signal_mat[twobase][label].extend(tmp_signal_mat[twobase][label])
                        context_array[twobase][label].extend(tmp_contexts[twobase][label])

        # Wait for all worker processes to finish
        for p in procs:
            p.join ()
       
    print('Finished extracting signals')
    if not training_tsv:
        tsv_out = tsv_output+'.unsorted.tmp'
        tmpfis = glob.glob(".".join(tsvname.split(".")[:-1])+"*.tmp[0-9]*")
        if nprocs > 1:
            print('Merging files...')
            try:
                os.remove(tsv_output)
            except OSError:
                pass
            for tmpfi in tmpfis:
                os.system('cat ' + tmpfi + ' >> '+ tsv_out )
                os.remove(tmpfi)
            os.system('sort -n -k2 ' + tsv_out + '| uniq > ' + tsv_output)
            os.remove(tsv_out)
        else:
            os.rename(tmpfis[0],tsv_output)

    if train: 
       print('Training...')
       #assert len(signal_array) > 5, 'insufficient data aligned to labeled positions for training'
       train_classifier(signal_mat,context_array,modelfile,classifier,plot_training) 
       print('Finished training') 


def main():
    #parse command line options
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Classify bases as methylated or unmethylated',prog='mCaller')
    all_or_some = parser.add_mutually_exclusive_group(required=True)
    all_or_some.add_argument('-p','--positions',type=str,required=False, help='file with a list of positions at which to classify bases (must be formatted as space- or tab-separated file with chromosome, position, strand, and label if training)')
    all_or_some.add_argument('-m','--motif',type=str,required=False, help='classify every base of type --base in the motif specified instead (can be single one-mer)')
    parser.add_argument('-r','--reference',type=str,required=True,help='fasta file with reference aligned to')
    parser.add_argument('-e','--tsv',type=str,required=True,help='tsv file with nanopolish event alignment')
    parser.add_argument('-f','--fastq',type=str,required=True,help='fastq file with nanopore reads')
    parser.add_argument('-t','--threads',type=int,required=False,help='specify number of processes (default = 1)',default=1)
    parser.add_argument('-b','--base',type=str,required=False,help='bases to classify as methylated or unmethylated (A or C, default A)',default='A')
    parser.add_argument('-n','--num_variables',type=int,required=False,help='change the length of the context used to classify (default of 6 variables corresponds to 11-mer context (6*2-1))',default=6)
    parser.add_argument('--train',action='store_true',required=False,help='train a new model (requires labels in positions file)',default=False)
    parser.add_argument('--training_tsv',type=str,required=False,help='mCaller output file for training') 
    parser.add_argument('-d','--modelfile',type=str,required=False,help='model file name')
    parser.add_argument('-s','--skip_thresh',type=int,required=False,help='number of skips to allow within an observation (default 0)',default=0)
    parser.add_argument('-q','--qual_thresh',type=float,required=False,help='quality threshold for reads (default none)',default=0)
    parser.add_argument('-c','--classifier',type=str,required=False,help='use alternative classifier: options = NN (default), RF, LR, or NBC (non-default may significantly increase runtime)',default='NN')
    parser.add_argument('--plot_training',action='store_true',required=False,help='plot probabilities distributions for training positions (requires labels in positions file and --train)',default=False)
    #parser.add_argument('--plot_violin',action='store_true',required=False,help='train a new model (requires labels in positions file)',default=False)
    parser.add_argument('-v','--version',action='version',help='print version',version='%(prog)s v1.0')
    args = parser.parse_args()

    if args.base == 'A':
        mod = 'm6A'
    elif args.base == 'C':
        mod = 'm5C' #TODO: test m4C
    else: 
        print('classification only available for A or C bases so far') 
        sys.exit(0)
 
    if not args.modelfile:
        modelfile = os.path.dirname(os.path.realpath(sys.argv[0])) + '/model_'+args.classifier+'_'+str(args.num_variables)+'_'+mod+'.pkl'
    else:
        modelfile = args.modelfile
    
    if not args.train:
        assert os.path.isfile(modelfile), 'model file not found at '+modelfile

    if args.motif and len(args.motif) == 1:
        base = args.motif
    else:
        base = args.base

    assert (args.skip_thresh < args.num_variables/2), 'too many skips with only '+str(args.num_variables)+' variables - try < half' 

    assert os.path.isfile(args.fastq), 'fastq file not found at '+args.fastq
    read2qual = extract_read_quality(args.fastq)

    if not args.training_tsv:
        training_tsv = None
    else:
        training_tsv = args.training_tsv

    try:
        num_refs = 0
        for ref in SeqIO.parse(args.reference,"fasta"):
            num_refs+=1
    except IOError:
        print('reference file missing')
        sys.exit(0)

    #distribute to multiple threads for main computations
    distribute_threads(args.positions,args.motif,args.tsv,read2qual,args.reference,num_refs,base,mod,args.threads,args.num_variables,
        args.train,modelfile,args.skip_thresh,args.qual_thresh,args.classifier,training_tsv,args.plot_training)

if __name__ == "__main__":
    main()


