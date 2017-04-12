#!/usr/bin/env python
#A program to classify bases as methylated or unmethylated based on long-range signals using the output from nanopolish
#Alexa McIntyre, 2016-2017

from collections import defaultdict
import numpy as np
import cPickle 
from extract_contexts import *
import sys
from Bio import SeqIO

def distribute_threads(positions_list,motif,tsvname,refname,multi_fasta,base,label,nprocs,nvariables):
    """ distributes list of genomic positions across processes then adds resulting signals to matrix"""
    #def worker(out_q, positions, tsv, fast5dir, refname, multi_fa, base, k=6, label=label):
    #    outtup = extract_features(positions, tsv, fast5dir, refname, multi_fa, base, k, label=label)
    #    out_q.put(outtup)

    if positions_list:
       nproc_min = min(nprocs, len(positions_list))
    else:
       nproc_min = nprocs

    if nproc_min == 1 and not multi_fasta:
        for ref in SeqIO.parse(refname,"fasta"):
           meth_fwd,meth_rev = methylate_references(str(ref.seq).upper(),base,motif=motif,positions=positions_list)
           extract_features(tsvname,meth_fwd,meth_rev,label=label,k=nvariables)
    """
    else:
        print('multiple threads not yet enabled')
        sys.exit(0)
        out_q = multiprocessing.Queue()
        chunksize = int(math.ceil(len(positions_list) / float(nprocs)))
        procs = []

        for i in range(nproc_min):
            p = multiprocessing.Process(
                    target=worker,
                    args=(out_q, positions_list[chunksize * i:chunksize * (i + 1)],
                          bamfi,fast5dir,refname,multi_fa,base,k,lab))
            procs.append(p)
            p.start()

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
            p.join()

    """
    print('Finished extracting signals')

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Classify bases as methylated or unmethylated')
    all_or_some = parser.add_mutually_exclusive_group(required=True)
    all_or_some.add_argument('-p','--positions',type=str,required=False, help='file with a list of positions at which to classify bases (default)')
    all_or_some.add_argument('-m','--motif',type=str,required=False, help='classify every base of type --base in the motif specified instead (can be single one-mer)')
    parser.add_argument('-r','--reference',type=str,required=True,help='fasta file with reference aligned to')
    parser.add_argument('-f','--tsv',type=str,required=True,help='tsv file with nanopolish alignment')
    parser.add_argument('-t','--threads',type=int,required=False,help='specify number of processes (default = 1)',default=1)
    parser.add_argument('-l','--label',type=str,required=False,help='label for bases in positions set (eg. A,C,m6A,m5C)',default='m6A')
    parser.add_argument('-b','--base',type=str,required=False,help='bases to classify as methylated or unmethylated (A or C)',default='A')
    parser.add_argument('-n','--num_variables',type=int,required=False,help='change the length of the context used to classify (default of 6 variables corresponds to 11-mer context (6*2-1))',default=6)
    parser.add_argument('--train',action='store_true',required=False,help='train a new model (requires labels)',default=False)
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    if args.version:
        print 'mCallerNP 0.1'
        sys.exit(0)

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

        #base_list_pos = [(i+1,'+') for i,base in enumerate(list(ref_seq.seq)[:-1]) if base == comp_dict[args.base] and ref_seq.seq[i+1] == 'G']
        #comp_list_pos = [(i+1,'-') for i,base in enumerate(list(ref_seq.seq)[:-1]) if base == args.base and ref_seq.seq[i+1] == 'C']

    #distribute to multiple threads for main computations
    distribute_threads(args.positions,args.motif,args.tsv,args.reference,multi_fasta,args.base,args.label,args.threads,args.num_variables)

if __name__ == "__main__":
    main()


