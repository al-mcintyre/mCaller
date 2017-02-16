#!/usr/bin/python
import sys 
import multiprocessing
import glob
import math
import random

import numpy as np
import pysam
import h5py

from model import *
from alignment import *
from classifier import *
from plots import *
import pyximport; pyximport.install()
from cython_viterbi import cython_viterbiish

def signal_at_pos(pos_list, bamname, fast5_dir, refname, seq, k=6, replace_base='M', label='m6A'):
    """main function: runs viterbi alignment and extracts the signals corresponding to a list of genomic positions"""
    #print 'starting list of',len(pos_list),'positions'
    bam = pysam.AlignmentFile(bamname, "rb")
    t_model = {}
    count = 0
    signal_matrix = []
    sequence_contexts = []
    for pos in pos_list:
        count += 1
        if len(pos_list) > 100 and count%(max(100,int(len(pos_list)/10))) == 0:
            print count    
        #fetch all reads that overlap a potentially methylated position in the reference 
        reads = bam.fetch(refname,pos[0]-1,pos[0])
        for read in reads:
            if read.query_name.split('_')[-1] == 'template' and ((read.is_reverse and pos[1] == '-') or (not read.is_reverse and pos[1] =='+')):
                if len(pos_list) < 4:
                    print read.query_name, pos[1], len(read.query_sequence)
                read_qual = read.query_qualities
                read_seq = read.query_sequence
                read_length = len(read_seq)
                alignment = read.get_aligned_pairs()
                if read.is_reverse:
                    ind_multiple = -1
                    match = 'T'
                else:
                    ind_multiple = 1
                    match = 'A'
                try:
                    alignment_A_pos = [(rr_ind[0],align_ind) for align_ind,rr_ind in enumerate(alignment) if rr_ind[1] == pos[0]-1]
                    read_pos, align_A_pos = alignment_A_pos[0]
                except TypeError: #TODO: error messages
                    break
                except IndexError:
                    break
                preA_rev_alignment = alignment[:align_A_pos][::-1]
                postA_alignment = alignment[align_A_pos:]
                # TODO: **first check the m6A location itself is not in a perfectly aligned kmer** #
                starting_anchor,sa_dist_from_A = find_perfect_match(preA_rev_alignment,k,read_seq,seq.seq)
                ending_anchor,ea_dist_from_A = find_perfect_match(postA_alignment,k,read_seq,seq.seq)
                if starting_anchor is not None and ending_anchor is not None:
                    seq_frag_pos = pos[0]-starting_anchor[-1][1]
                    seq_frag = Sequence(seq.seq[starting_anchor[-1][1]:ending_anchor[-1][1]+1])
                    seq_frag.methylate(seq_frag_pos,replace_base)
                    seq_frag.rev_comp()       
                    context = seq_frag.seq[seq_frag_pos-k:seq_frag_pos+k-1]
                    #print context, pos[1]
                    if len(seq_frag.seq) > 50: #TODO: enable realignment for sequences without close perfect kmer matches
                        break 

                    fast5 = glob.glob(fast5_dir+'/'+'_'.join(read.query_name.split('_')[:-1])+'*')
                    with h5py.File(fast5[0],'r') as hf:        
                        if not t_model:
                            t_model = extract_model(hf)

                        transition_probabilities = extract_transitions(hf)

                        if read.is_reverse:
                            s, e = starting_anchor[0][0]+1, ending_anchor[-1][0]+1
                            pass_seq = seq_frag.revcomp
                            pass_mseq = seq_frag.revmethseq
                        else:
                            s, e = starting_anchor[-1][0], ending_anchor[0][0]
                            pass_seq = seq_frag.seq
                            pass_mseq = seq_frag.methseq

                        signal_frag = extract_signal_fragment(s,e,hf,ind_multiple,read_length)
                        signal_frag.transform()

                        #realigned_signal = cython_viterbiish(signal_frag.means,signal_frag.var,pass_seq,np.array(transition_probabilities),t_model.model,k,pass_mseq)
                        realigned_signal = viterbiish(signal_frag,pass_seq,transition_probabilities,t_model,k,pass_mseq)
                        kmer_array = update_signal_matrix(realigned_signal,t_model.model)
                        signal_matrix.append(kmer_array)
                        sequence_contexts.append(context)
                        #if len(pos_list) < 4:
                            #plot_realignment(realigned_signal,t_model)
    signal_labels = [label]*len(signal_matrix)                        

    return signal_matrix,signal_labels,sequence_contexts


def distribute_genome_positions(positions_list,bamfi,fast5dir,refname,refsequence,lab,nprocs,k=6):
    """ distributes list of genomic positions across processes then adds resulting signals to matrix"""
    def worker(out_q, positions, bamfile, fast5dir, refname, refsequence, k=6,label=lab):
        outtup  = signal_at_pos(positions, bamfile, fast5dir, refname, refsequence, k,label=lab)
        out_q.put(outtup)
    
    nproc_min = min(nprocs, len(positions_list))

    if nproc_min == 1:
        signal_mat,label_array,surrounding_contexts = signal_at_pos(positions_list, bamfi, fast5dir, refname, refsequence,k,label=lab)
        
    else:
        out_q = multiprocessing.Queue()
        chunksize = int(math.ceil(len(positions_list) / float(nprocs)))
        procs = []

        for i in range(nproc_min):
            p = multiprocessing.Process(
                    target=worker,
                    args=(out_q, positions_list[chunksize * i:chunksize * (i + 1)],
                          bamfi,fast5dir,refname,refsequence,k,lab))
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
        
    print('Finished extracting signals')
    return(signal_mat,label_array,context_array)

def main():
    #parse command line options
    from argparse import ArgumentParser
    parser = ArgumentParser(description='Classify adenines as methylated or unmethylated')
    all_or_some = parser.add_mutually_exclusive_group(required=True)
    all_or_some.add_argument('-p','--positions',type=str,required=False, help='file with a list of positions at which to classify adenines (default)')
    all_or_some.add_argument('-w','--whole-genome',action='store_true',required=False, help='classify every adenine in the genome instead')
    parser.add_argument('-r','--reference',type=str,required=True,help='fasta file with reference aligned to')
    parser.add_argument('-b','--bam',type=str,required=True,help='bam file with alignment (must be sorted and indexed)')
    parser.add_argument('-d','--directory',type=str,required=True,help='directory of fast5 files included in the alignment')
    parser.add_argument('-t','--threads',type=int,required=False,help='specify number of processes (default = 1)',default=1)
    parser.add_argument('-l','--label',type=str,required=False,help='label for adenines in positions set (A or m6A)',default='m6A')
    parser.add_argument('--train',action='store_true',required=False,help='train a new model (requires labels)',default=False)
    parser.add_argument('-v','--version',action='store_true',required=False,help='print version')
    args = parser.parse_args()

    try:
        reffi = open(args.reference,'r').read()
        assert len(reffi.split('>')) <= 2
        ref_name = reffi.split('\n')[0][1:]
        ref_seq = Sequence(''.join(reffi.split('\n')[1:])) #TODO: change to accept fastas with multiple sequences
    except IOError:
        print('reference file missing')
        sys.exit(0)
    except AssertionError:
        print('try a fasta with a single sequence, please! we are still working on accepting multifastas')
        sys.exit(0)
    
    if args.positions:
        try:
            pos_list = [(int(x.split()[0]),x.split()[1]) for x in open(args.positions,'r').read().split('\n') if x != '' and x.split()[1] in set(['-','+'])]
            if len(pos_list) == 0:
                raise TypeError
        except TypeError:
            print('positions file does not contain two columns with tab- or space-separated positions and strands')
            sys.exit(0)
        except IOError:
            print('positions file missing')
            sys.exit(0)
    elif args.whole_genome:
        A_list_pos = [(i+1,'+') for i,base in enumerate(list(ref_seq.seq)) if base == 'A']
        T_list_pos = [(i+1,'-') for i,base in enumerate(list(ref_seq.seq)) if base == 'T']
        pos_list = A_list_pos+T_list_pos #random.sample(A_list_pos+T_list_pos,100) #FOR TESTING
    num_threads = args.threads
    signals,labels,contexts = distribute_genome_positions(pos_list,args.bam,args.directory,ref_name,ref_seq,args.label,args.threads)
    meth_calls = random_forest(signals,labels,args.train)   #TODO: save genomic positions of calls 
    egfi = open('lambda_20160907_A.pkl','wb')
    cPickle.dump(zip(signals,labels,contexts),egfi)
    egfi.close()
    print(len([x for x in meth_calls if x[0] > 0.5])*100./len(meth_calls), len(meth_calls))
    print(np.mean([x[0] for x in meth_calls]))

if __name__ == "__main__":
    main()
