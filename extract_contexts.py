from collections import defaultdict
from Bio import SeqIO
import scipy.stats as ss
import numpy as np
import cPickle 
import sys
import re

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
   #print motif, meth_base, meth_position
   if meth_position:
      meth_motif = motif[:meth_position]+'M'
      if meth_position < len(motif)-1:
         meth_motif = meth_position+motif[meth_position+1:]
   else:
      meth_motif = 'M'.join(motif.split(meth_base))
   meth_seq = ref_seq.replace(motif,meth_motif)
   #ref_motif_segs = ref_seq.split(motif)
   #meth_seq = meth_motif.join(ref_motif_segs)
   #print len(meth_seq.split(meth_motif))-1, motif+' positions found'
   return meth_seq

#change specified positions to M in reference sequence
#@profile
def methylate_positions(ref_seq,positions,meth_base):
   meth_seq = ref_seq
   #print('next sequence')
   count = 0
   for pos in positions: #changed to 0-based - else have to subtract from pos
      #print meth_seq[pos-6:pos+5]
      if meth_seq[pos] == meth_base:
         meth_seq = meth_seq[:pos]+'M'+meth_seq[pos+1:]
         count+=1
         #print meth_seq[pos-6:pos+5]
      else:
         print count, meth_seq[pos-5:pos+6], pos
         print 'bad methylation - check reference positions are 0-based - quitting thread now'
         sys.exit(0)
   #print count, 'positions methylated in one strand' 
   return meth_seq

#extract signals around methylated positions from tsv
#@profile
def methylate_references(ref_seq,base,motif=None,positions=None,train=False,contig=None):
   #print 'sequence length', len(ref_seq)
   if motif:
      meth_fwd = methylate_motifs(ref_seq,motif,base)
      meth_rev = methylate_motifs(ref_seq,revcomp(motif),base_comps[base])
      #print len(meth_fwd.split('M')),'Ms in methylated sequence'
   elif positions:
      fwd_pos = [int(pos.split()[1]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[2] == '+' and pos.split()[0] == contig]
      rev_pos = [int(pos.split()[1]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[2] == '-' and pos.split()[0] == contig]
      meth_fwd = methylate_positions(ref_seq,fwd_pos,base)
      meth_rev = methylate_positions(ref_seq,rev_pos,base_comps[base])
   else:  
      print 'no motifs or positions specified'
      sys.exit(0)
   return meth_fwd,meth_rev

#@profile
def find_and_methylate(refname,contigname,base,motif,positions_list):
    for ref in SeqIO.parse(refname,"fasta"):
        contigid = ref.id
        if contigid == contigname:
            #print 'contig =',contigid
            #print 'methylating contig', contigname
            meth_fwd,meth_rev = methylate_references(str(ref.seq).upper(),base,motif=motif,positions=positions_list,contig=contigname)
            #print 'methylated reference complete'
            return meth_fwd,meth_rev
        #else:
            #print 'not methylating',contigid

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
        model = cPickle.load(modfi)
        modfi.close()
    else:
        tsv_output = '.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.train.tmp'+str(startline)
        signals,labels,contexts = [],[],[]

    towrite = []
    #save only one set of adjoining methylated positions at a time - once the set complete, write the positions to a file 
    #tsv format: ecoli   805 CGCCAT  cc1da58e-3db3-4a4b-93c2-c78e1dbe6aba:1D_000:template    t   1   102.16  0.963   0.00175 CGCCAT  102.23  1.93    -0.03   101.973,100.037,102.403,101.758,104.338,102.618,101.973
    with open(tsv_input,'r') as tsv:
        tsv.seek(max(startline-500,0))
        tsv.readline() #to start new line
        while tsv.tell() <= endline-500:
            lines = tsv.readlines(8000000)
            #if tsv.tell() > 80374772633:
                #print tsv.tell()
            for line in lines:
                chrom, read_pos, read_kmer, read_name, x, read_ind, event_current, event_sd, y, ref_kmer, model_current, ref_sd, z  = line.split('\t')
                if chrom != last_contig:
                    #print 'loading new contig',chrom
                    try:
                        #print 'new chrom',chrom,'but old was',last_contig,'and byte is',tsv.tell(),'with endbyte',endline,'for thread starting at',startline
                        #print line
                        meth_fwd,meth_rev = find_and_methylate(fasta_input,chrom,base,motif,positions_list)
                        #print 'finished loading reference contig.' #,len(meth_fwd.split('M')),'positions to examine' 
                        last_contig = chrom
                    except ValueError:
                        print 'Error: could not find sequence for reference contig',chrom
                        continue
                if read_name != last_read:
                    first_read_ind = int(read_ind) 
                if (read2qual[read_name] < qual_thresh) or ref_kmer == 'NNNNNN':
                    continue
                if (read_name != last_read and read_kmer == ref_kmer) or (read_name == last_read and int(read_ind) > first_read_ind): #takes into account complementary palindromes
                    rev = False
                    meth_ref = meth_fwd
                else:
                    rev = True
                    meth_ref = meth_rev
                read_pos = int(read_pos)
                reference_kmer = meth_ref[read_pos:read_pos+k]
                #print read_name,rev,ref_kmer,reference_kmer

                #if finished context for previous potentially modified position, save and reset
                if mpos and ((read_pos >= mpos+1 and read_name == last_read) or (read_name != last_read)):
     
                    #write to file
                    num_skips = len([x for x in diff_col if x == []])
                    if num_skips <= skip_thresh: #accept max number of skips within an observation
                        if num_skips> 0:
                            w_skips.add((last_read,mpos))
                        diffs = [np.mean(kmer_pos) if kmer_pos!=[] else 0 for kmer_pos in diff_col] 
                        if not last_rev:
                            diffs = diffs[::-1]
                        diffs = diffs+[read2qual[last_read]] #tried -18 to scale quality.. 
                        context = revcomp(meth_ref[mpos-k+1:mpos+k],last_rev)
                        if not train:
                            mod_prob = model.predict_proba([diffs]) #TODO: call model only when batch ready to write
                            if mod_prob[0][1] >= 0.5: 
                                label = 'm6A' #TODO: ensure correct direction + label unmeth/meth as appropriate 
                            else:
                                label = 'A' 
                            label = label+'\t'+str(np.round(mod_prob[0][1],2))
                        else:
                            mod_prob = ''
                            label = pos_label[(chrom,mpos,strand(last_rev))] 
                            signals.append(diffs)
                            labels.append(label)
                            contexts.append(context)
                        towrite.append([chrom,last_read,str(mpos),context,','.join([str(diff) for diff in diffs]),strand(last_rev),label])
                        #outfi.write(chrom+'\t'+last_read+'\t'+str(mpos)+'\t'+context+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev)+'\t'+label+'\n')
                        last_info = last_read+'\t'+str(mpos)+'\t'+context+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev)+'\t'+label
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
                        mspacing = mpos - last_mpos
                        last_pos_in_kmer = pos_in_kmer
                        diffs = [[] for i in range(mspacing)] + diff_col[:-mspacing]
                        diff_col = diffs
                        if len(diff_col) != k:
                            try:
                                print last_info
                            except:
                                pass
                            print reference_kmer,mpos,read_pos,read_pos>mpos,read_name,last_read,diff_col,mspacing
                            diff_col = [[] for i in range(k)]
                            #break

            #if modified base in reference, save surrounding context to call that position
            #if len([x for x in reference_kmer if x == 'M']) >= 1:
                if 'M' in set(list(reference_kmer)):
                    pos_in_kmer = [i for i,x in enumerate(list(reference_kmer)) if x == 'M'][0]
                    #if new read, reset differences variable and proceed
                    if mpos and read_name != last_read:
                        mpos = None
                        diff_col = [[] for i in range(k)]
                    #if new read or new position
                    if not mpos: #TODO: reject any positions too close to beginning of read automatically
                        mpos = read_pos+pos_in_kmer
                    last_pos_in_kmer = pos_in_kmer
                    last_read = read_name
                    last_rev = rev
                    try:
                        diff_col[pos_in_kmer].append(float(event_current)-float(model_current))
                    except IndexError:
                        print diff_col, mpos, read_pos, reference_kmer, pos_in_kmer
                        diff_col = [[] for i in range(k)]
                        diff_col[pos_in_kmer].append(float(event_current)-float(model_current))
                        #break
                    last_pos = read_pos 
                    #print mpos, reference_kmer, read_pos, diff_col
         
                elif mpos:
                    mpos = None
                    diff_col = [[] for i in range(k)]

    writefi(towrite,tsv_output)
    #with open(tsv_output,'a') as outfi:
    #        for entry in towrite:
    #            outfi.write('\t'.join(entry)+'\n')

    print 'thread finished processing...:'
    print num_observations,'observations'
    num_pos = len(pos_set)
    print num_pos,'positions'
    #cPickle.dump(diffs_by_read,open(tsv_input+'.pkl','wb'))
    print len(multi_meth_pos_set),'regions with multiple methylated bases'
    print len(w_skips), 'observations with skips included'
    print len(skipped_skips), 'observations with too many skips'
    if train:
        return signals, labels, contexts

