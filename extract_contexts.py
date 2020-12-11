from collections import defaultdict
from Bio import SeqIO
import scipy.stats as ss
import numpy as np
#import cPickle
import pickle #_pickle as cPickle
import sys
import re
import os

base_comps = {'A':'T','C':'G','T':'A','G':'C','N':'N','M':'M'}

#@profile
def comp(seq,base_comps=base_comps):
   return ''.join([base_comps[nt] for nt in list(seq)])

#@profile
def revcomp(seq,rev=True):
   if not rev:
      return seq
   else:
      return ''.join(list(comp(seq))[::-1])

#@profile
def strand(rev):
   if rev:
      return '-'
   else:
      return '+'

#find positions of motifs (eg. CG bases) in reference sequence and change to M
#@profile
def methylate_motifs(ref_seq,motif,meth_base,meth_position=None): 
   if meth_position:
      meth_motif = motif[:meth_position]+'M'
      if meth_position < len(motif)-1:
         meth_motif = meth_position+motif[meth_position+1:]
   else:
      meth_motif = 'M'.join(motif.split(meth_base))
   meth_seq = ref_seq.replace(motif,meth_motif)
   return meth_seq

#change specified positions to M in reference sequence
#@profile
def methylate_positions(ref_seq,positions,meth_base):
   meth_seq = ref_seq
   count = 0
   for pos in positions: #changed to 0-based - else have to subtract from pos
      if meth_seq[pos] == meth_base or meth_seq[pos] == 'M':
         meth_seq = meth_seq[:pos]+'M'+meth_seq[pos+1:]
         count+=1
      else:
         print('Base {} does not correspond to methylated base - check reference positions are 0-based - quitting thread now'.format(pos))
         sys.exit(0)
   #print count, 'positions methylated in one strand' 
   return meth_seq

#extract signals around methylated positions from tsv
#@profile
def methylate_references(ref_seq,base,motif=None,positions=None,train=False,contig=None):
   if not positions and motif:
      meth_fwd = methylate_motifs(ref_seq,motif,base)
      meth_rev = methylate_motifs(ref_seq,revcomp(motif),base_comps[base])
      #print len(meth_fwd.split('M')),'Ms in methylated sequence'
   elif positions:
      fwd_pos = [int(pos.split()[1]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[2] == '+' and pos.split()[0] == contig]
      rev_pos = [int(pos.split()[1]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[2] == '-' and pos.split()[0] == contig]
      meth_fwd = methylate_positions(ref_seq,fwd_pos,base)
      meth_rev = methylate_positions(ref_seq,rev_pos,base_comps[base])
   else:  
      print('no motifs or positions specified')
      sys.exit(0)
   return meth_fwd,meth_rev

#@profile
def find_and_methylate(refname,contigname,base,motif,positions_list):
    for ref in SeqIO.parse(refname,"fasta"):
        contigid = ref.id
        if contigid == contigname:
            meth_fwd,meth_rev = methylate_references(str(ref.seq).upper(),base,motif=motif,positions=positions_list,contig=contigname)
            return meth_fwd,meth_rev

def writefi(data,fi):
    with open(fi,'a') as outfi:
        for entry in data:
            outfi.write('\t'.join(entry)+'\n')

def adjust_scores(context_dict,context,diffs,prob,k):
    if context in context_dict['m6A']:
        hmm_score = 1-(1/np.prod([ss.norm(context_dict['m6A'][context]['mean'][i],context_dict['m6A'][context]['sd'][i]).pdf(diffs[i]) for i in range(k)]))
        correlation_score = ss.stats.pearsonr(context_dict['m6A'][context]['mean'],diffs)[0]
        if context in context_dict['A']:
            correlation_diff = correlation_score - ss.stats.pearsonr(context_dict['A'][context]['mean'],diffs)[0]
            frac_meth = context_dict['m6A'][context]['num']*1./context_dict['A'][context]['num']
        else:
            frac_meth = 1
        representation_score = prob + 1 - frac_meth #increases score for contexts not included in methylation training set

def base_models(base,twobase=False):
    if base == 'A' and twobase:
        base_model = {'MG':'MG','MC':'MH','MA':'MH','MT':'MH','MM':'MH','MH':'MH','AT':'MH','AC':'MH','AG':'MG','AT':'MH','AA':'MH','AM':'MH'} #TODO: fix error where sites not methylated
    else:
        base_model = {'M'+nextb:'general' for nextb in ['A','C','G','T','M']}
        base_model.update({'A'+nextb:'general' for nextb in ['A','C','G','T','M']})
        base_model.update({'T'+nextb:'general' for nextb in ['A','C','G','T','M']}) #shouldn't be necessary
    return(base_model)

#determine difference between measurements and model for bases surrounding methylated positions 
#@profile
def extract_features(tsv_input,fasta_input,read2qual,k,skip_thresh,qual_thresh,modelfile,classifier,startline,endline=None,train=False,pos_label=None,base=None,motif=None,positions_list=None):

    #set position variables
    last_read,last_pos,last_pos_in_kmer,last_read_num = '',0,k,0
    last_contig = None
    #set count variables 
    num_observations,w_skips,skipped_skips,pos_set,multi_meth_pos_set,read_set = 0,set(),set(),set(),set(),set()
    #set tracking variables for observation
    mpos = None
    diff_col = [[] for xi in range(k)]

    if not train:
        tsv_output = '.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.tmp'+str(startline)
        modfi = open(modelfile,'rb')
        model = pickle.load(modfi,encoding='latin')
        modfi.close()
        if type(model) != dict:
            model = {'general':model} #for compatibility with previously trained model
            twobase = False
        else:
            twobase = True
        base_model = base_models(base,twobase)
    else:
        base_model = base_models(base,False) #or set to False?
        tsv_output = '.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.train.tmp'+str(startline)
        signals,contexts = {bm:{} for bm in base_model.values()},{bm:{} for bm in base_model.values()}

    towrite = []
    #save only one set of adjoining methylated positions at a time - once the set complete, write the positions to a file 
    #tsv format: ecoli   805 CGCCAT  cc1da58e-3db3-4a4b-93c2-c78e1dbe6aba:1D_000:template    t   1   102.16  0.963   0.00175 CGCCAT  102.23  1.93    -0.03   101.973,100.037,102.403,101.758,104.338,102.618,101.973
    with open(tsv_input,'r') as tsv:
        tsv.seek(max(startline-500,0))
        linepos = max(startline-500,0)
        #startline, endline, and linepos are in characters -- previously used tsv.tell(), but incompatible with python3
        while linepos <= endline-500:  
            #print('current position',linepos)
            lines = tsv.readlines(8000000) #TODO: why 8M? reasonable size for memory consumption, but could change
            for line in lines:
                linepos += len(line)
                try:
                    chrom, read_pos, read_kmer, read_name, x, read_ind, event_current, event_sd, y, ref_kmer, model_current, ref_sd  = line.split()[:12]
                except ValueError:
                    continue
                    
                if chrom != last_contig:
                    try:
                        meth_fwd,meth_rev = find_and_methylate(fasta_input,chrom,base,motif,positions_list)
                        last_contig = chrom
                    except TypeError: #ValueError
                        print('Error: could not find sequence for reference contig',chrom)
                        continue
                if read_name != last_read:
                    first_read_ind = int(read_ind) 
                try: 
                    qual = read2qual[read_name]
                except KeyError:
                    qual = read2qual[read_name.split(':')[0].split('_')[0]]
                if (qual < qual_thresh) or ref_kmer == 'NNNNNN':
                    continue
                if (read_name != last_read and read_kmer == ref_kmer) or (read_name == last_read and int(read_ind) > first_read_ind): #takes into account complementary palindromes - temporarily sets new reads to positive strand
                    rev = False
                    meth_ref = meth_fwd
                else:
                    rev = True
                    meth_ref = meth_rev
                read_pos = int(read_pos)
                reference_kmer = meth_ref[read_pos:read_pos+k]

                #if finished context for previous potentially modified position, save and reset
                if mpos and ((read_pos >= mpos+1 and read_name == last_read) or (read_name != last_read)):
     
                    #write to file
                    num_skips = len([x for x in diff_col if x == []])
                    if num_skips <= skip_thresh: #accept max number of skips within an observation
                        if num_skips > 0:
                            w_skips.add((last_read,mpos))
                        diffs = [np.mean(kmer_pos) if kmer_pos!=[] else 0 for kmer_pos in diff_col] 
                        if not last_rev:
                            diffs = diffs[::-1]
                        try:
                            last_qual = read2qual[last_read]
                        except KeyError:
                            last_qual = read2qual[last_read.split(':')[0].split('_')[0]]
                        diffs = diffs+[last_qual] 
                        context = revcomp(last_ref[mpos-k+1:mpos+k],last_rev)
                        if context[int(len(context)/2)] == 'M':
                            try: 
                                twobase_model = base_model[context[int(len(context)/2):int(len(context)/2)+2]]
                                if not train:
                                    mod_prob = model[twobase_model].predict_proba([diffs]) #TODO: call model only when batch ready to write
                                    if mod_prob[0][1] >= 0.5: 
                                        if base == 'A':
                                            label = 'm6A' #TODO: ensure correct direction + label unmeth/meth as appropriate 
                                        else:
                                            label = 'm'+base
                                    else:
                                        label = base 
                                    label = label+'\t'+str(np.round(mod_prob[0][1],2))
                                else:
                                    mod_prob = ''
                                    label = pos_label[(chrom,mpos,strand(last_rev))] 
                                    if label not in signals[twobase_model]:
                                        signals[twobase_model][label] = []
                                        contexts[twobase_model][label] = []
                                    signals[twobase_model][label].append(diffs)
                                    contexts[twobase_model][label].append(context)
                                towrite.append([chrom,last_read,str(mpos),context,','.join([str(diff) for diff in diffs]),strand(last_rev),label])
                                last_info = last_read+'\t'+str(mpos)+'\t'+context+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev)
                            except (IndexError,KeyError) as e:
                                print(last_read+'\t'+str(mpos)+'\t'+context+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev),'- Index or Key Error')
                                print(model.keys(), base_model.keys(), context[int(len(context)/2):int(len(context)/2)+2])
                                print(e)
                                print(model[twobase_model].predict_proba([diffs]))
                                sys.exit(0)
                        else:
                            print(last_read+'\t'+str(mpos)+'\t'+context+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev))
                            print(read_name, rev, last_read, last_rev, last_first)
                            print(read_kmer,reference_kmer, ref_kmer, last_pos_in_kmer, mspacing, pos_in_kmer)
                            sys.exit(0)
                        num_observations += 1
                        if num_observations%5000 == 0:
                            writefi(towrite,tsv_output)
                            towrite = []
                        pos_set.add(mpos)
                        read_set.add(last_read)
                        if len(read_set)%1000 == 0 and len(read_set) > last_read_num:
                            #print len(read_set), 'reads examined'
                            last_read_num = len(read_set)
                    else:
                        skipped_skips.add((last_read,mpos))
       
                  #reset variables
                    if len(reference_kmer.split('M')) < 2 or read_name != last_read or read_pos > mpos+skip_thresh+1: #allow no more than skip_thresh skips 
                        diff_col = [[] for i in range(k)] 
                        mpos = None
                        last_pos_in_kmer = k 
                    else: 
                        if reference_kmer[0] != 'M':
                            multi_meth_pos_set.add((last_read,mpos))
                        last_mpos = mpos
                        pos_in_kmer = len(reference_kmer.split('M')[0])
                        mpos = read_pos + pos_in_kmer
                        mspacing = min(k,mpos - last_mpos)
                        last_pos_in_kmer = pos_in_kmer
                        last_diff_col = diff_col
                        diffs = [[] for i in range(mspacing)] + diff_col[:-mspacing]
                        diff_col = diffs
                        if len(diff_col) != k:
                            try:
                                print(last_info,'- n diffs off')
                            except:
                                pass
                            #GGCGCM 613883 613878 False 2289b392-746e-4fa0-8226-d3ac661c9620_Basecall_2D_template 2289b392-746e-4fa0-8226-d3ac661c9620_Basecall_2D_template [[], [], [], [], [], [], []] 7

                            print(reference_kmer,last_mpos,mpos,mspacing,read_pos,read_pos-last_mpos,read_name,last_read,diff_col,mspacing, last_diff_col, last_diff_col[:-mspacing])
                            diff_col = [[] for i in range(k)]
                            sys.exit(0)

            #if modified base in reference, save surrounding context to call that position
                if 'M' in set(list(reference_kmer)):
                    pos_in_kmer = [i for i,x in enumerate(list(reference_kmer)) if x == 'M'][0]
                    #if new read, reset differences variable and proceed
                    if mpos:
                        if read_name != last_read:
                            mpos = None
                            diff_col = [[] for i in range(k)]
                        elif rev != last_rev:
                            mpos = None
                    #if new read or new position
                    if not mpos: 
                        mpos = read_pos+pos_in_kmer
                    last_pos_in_kmer = pos_in_kmer
                    last_read = read_name
                    last_rev = rev
                    last_first = first_read_ind
                    last_ref = meth_ref
                    diff_col[pos_in_kmer].append(np.round(float(event_current)-float(model_current),4))
                    last_pos = read_pos 
         
                elif mpos:
                    mpos = None
                    diff_col = [[] for i in range(k)]

    writefi(towrite,tsv_output)

    print('thread finished processing...:')
    print('%d observations' %num_observations)
    num_pos = len(pos_set)
    print('%d positions' %num_pos)
    print('%d regions with multiple methylated bases' %len(multi_meth_pos_set))
    print('%d observations with skips included' %len(w_skips))
    print('%d observations with too many skips' %len(skipped_skips)) 
    if train:
        return signals, contexts

