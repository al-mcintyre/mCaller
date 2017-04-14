#!/usr/bin/env python
#A program to classify bases as methylated or unmethylated based on long-range signals using the output from nanopolish
#Alexa McIntyre, 2016-2017

from collections import defaultdict
import numpy as np
import cPickle 
import sys
#import time
import os
import multiprocessing
from Bio import SeqIO

from extract_contexts import *
from train_model import train_classifier
from read_qual import extract_read_quality

def distribute_threads(positions_list,motif,tsvname,read2qual,refname,multi_fasta,base,mod,nprocs,nvariables,train,modelfile,skip_thresh,qual_thresh,classifier):
    """ distributes list of genomic positions across processes then adds resulting signals to matrix if training"""

    if not multi_fasta:

       for ref in SeqIO.parse(refname,"fasta"):
          meth_fwd,meth_rev = methylate_references(str(ref.seq).upper(),base,motif=motif,positions=positions_list)

          if nprocs == 1:   
             if not train:
                extract_features(tsvname,read2qual,meth_fwd,meth_rev,nvariables,skip_thresh,qual_thresh,modelfile,classifier,0) #TODO: implement quality thresholding
             else:
                signal_mat, label_array, context_array = extract_features(tsvname,read2qual,meth_fwd,meth_rev,nvariables,modelfile,skip_thresh,qual_thresh,classifier,0,train)

          elif nprocs > 1 and not multi_fasta:
             def worker(out_q,tsvname,read2qual,meth_fwd,meth_rev,nvariables,skip_thresh,qual_thresh,modelfile,classifier,startline,endline,train): 
                outtup = extract_features(tsvname,read2qual,meth_fwd,meth_rev,nvariables,skip_thresh,qual_thresh,modelfile,classifier,startline,endline=endline,train=train)
                out_q.put(outtup)

             def countlines(filename):
                return sum(1 for _ in open(filename))

             nlines = countlines(tsvname)
        
             out_q = multiprocessing.Queue()
             chunksize = int(math.ceil(nlines / float(nprocs)))
             procs = []

             for i in range(nprocs):
                p = multiprocessing.Process(
                   target=worker,
                   args=(out_q,tsvname,read2qual,meth_fwd,meth_rev,nvariables,skip_thresh,qual_thresh,modelfile,classifier,chunksize*i,chunksize*(i+1),train))
                procs.append(p)
                p.start()

             if train:
                # Collect all results into a signal matrix and an array of labels
                signal_mat = []
                label_array = []
                context_array = []
                for i,proc in enumerate(procs):
                   tmp_signal_mat,tmp_label_array,tmp_contexts = out_q.get()
                   print 'updating with results from process',i
                   signal_mat.extend(tmp_signal_mat)
                   label_array.extend(tmp_label_array)
                   context_array.extend(tmp_contexts)

             # Wait for all worker processes to finish
             for p in procs:
                p.join ()
       
    print('Finished extracting signals')

    if train: 
       train_classifier(signal_mat,label_array,context_array,modelfile,classifier) 
       print('Finished training') 


def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Classify bases as methylated or unmethylated')
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
    parser.add_argument('-d','--modelfile',type=str,required=False,help='model file name')
    parser.add_argument('-s','--skip_thresh',type=int,required=False,help='number of skips to allow within an observation (default 0)',default=0)
    parser.add_argument('-q','--qual_thresh',type=float,required=False,help='quality threshold for reads (under development, please sort your own reads for now)',default=0)
    parser.add_argument('-c','--classifier',type=str,required=False,help='use alternative classifier: options = NN (default) RF, LR, or NBC',default='NN')
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.1'
        sys.exit(0)

    if args.base == 'A':
        mod = 'm6A'
    elif args.base == 'C':
        mod = 'm5C' #TODO: test m4C
    else: 
        print('classification only available for A or C bases so far') 
        sys.exit(0)
 
    if not args.modelfile:
        args.modelfile = 'mod.'+args.classifier+'.'+str(args.num_variables)+'.pkl'
    
    if not args.train:
        assert os.path.isfile(args.modelfile), 'model file not found at '+args.modelfile

    assert (args.skip_thresh < args.num_variables/2), 'too many skips with only '+str(args.num_variables)+' variables - try < half' 

    assert os.path.isfile(args.fastq), 'fastq file not found at '+args.fastq
    read2qual = extract_read_quality(args.fastq)

    try:
        num_refs = 0
        for ref in SeqIO.parse(args.reference,"fasta"):
            num_refs+=1
        if num_refs < 2:
            multi_fasta = False
        else:
            multi_fasta = True
            print('multi-fasta base mod calling under development, try again later')
            sys.exit(0)
    except IOError:
        print('reference file missing')
        sys.exit(0)

    #distribute to multiple threads for main computations
    distribute_threads(args.positions,args.motif,args.tsv,read2qual,args.reference,multi_fasta,args.base,mod,args.threads,args.num_variables,
        args.train,args.modelfile,args.skip_thresh,args.qual_thresh,args.classifier)

if __name__ == "__main__":
    main()


