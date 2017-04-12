#A program to classify bases as methylated or unmethylated based on long-range signals using the output from nanopolish
#Alexa McIntyre, 2016-2017

from collections import defaultdict
import numpy as np
import cPickle 
import sys

base_comps = {'A':'T','C':'G','T':'A','G':'C','N':'N','M':'M'}

def comp(seq,base_comps=base_comps):
   return ''.join([base_comps[nt] for nt in list(seq)])

def revcomp(seq,rev=True):
   if not rev:
      return seq
   else:
      return ''.join(list(comp(seq))[::-1])

def strand(rev):
   if rev:
      return '-'
   else:
      return '+'

#find positions of motifs (eg. CG bases) in reference sequence and change to M
def methylate_motifs(ref_seq,motif,meth_base,meth_position=None):
   l_motif = len(motif)
   if meth_position:
      meth_motif = motif[:meth_position]+'M'
      if meth_position < len(motif)-1:
         meth_motif = meth_position+motif[meth_position+1:]
   else:
      meth_motif = 'M'.join(motif.split(meth_base))
   ref_motif_segs = ref_seq.split(motif)
   meth_seq = meth_motif.join(ref_motif_segs)
   print len(ref_motif_segs)-1, 'CG positions (not necessarily singletons) in the forward strand'
   return meth_seq

#change specified positions to M in reference sequence
def methylate_positions(ref_seq,positions,meth_base):
   meth_seq = ref_seq
   #print('next sequence')
   count = 0
   for pos in positions:
      #print meth_seq[pos-6:pos+5]
      if meth_seq[pos-1] == meth_base:
         meth_seq = meth_seq[:pos-1]+'M'+meth_seq[pos:]
         count+=1
         #print meth_seq[pos-6:pos+5]
      else:
         print count, meth_seq[pos-6:pos+5]
         print 'bad methylation'
         sys.exit(0)
   return meth_seq

#extract signals around methylated positions from tsv
def methylate_references(ref_seq,base,motif=None,positions=None):
   ref_seq = ref_seq
   if motif:
      meth_fwd = methylate_motifs(ref_seq,motif,base)
      meth_rev = methylate_motifs(ref_seq,comp(motif),base_comps[base])
   elif positions:
      fwd_pos = [int(pos.split()[0]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[1] == '+']
      rev_pos = [int(pos.split()[0]) for pos in open(positions,'r').read().split('\n') if len(pos.split()) > 1 and pos.split()[1] == '-']
      meth_fwd = methylate_positions(ref_seq,fwd_pos,base)
      meth_rev = methylate_positions(ref_seq,rev_pos,base_comps[base])
   else:  
      print 'no motifs or positions specified'
      sys.exit(0)
   return meth_fwd,meth_rev

#determine difference between measurements and model for bases surrounding methylated positions 
def extract_features(tsv_input,meth_fwd,meth_rev,label,k=6):
   #print meth_fwd[14827-k+1:14827+k]
   #print meth_rev[14828-k+1:14828+k]
   last_read = ''
   last_pos = 0
   last_pos_in_kmer = k
   mpos = None
   diff_col = [[] for xi in range(k)]
   num_observations = 0
   w_skips = set()
   pos_set = set()
   read_set = set()
   firstline = True
   last_read_num = 0
   multi_meth_pos_set = set()
   skipped_skips = set()
   tsv_output = '.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.'+label
   print tsv_output

   #only save one set of adjoining methylated positions at a time - once the set complete, write the positions to a file 
   #format: ecoli	805	CGCCAT	cc1da58e-3db3-4a4b-93c2-c78e1dbe6aba:1D_000:template	t	1	102.16	0.963	0.00175	CGCCAT	102.23	1.93	-0.03	101.973,100.037,102.403,101.758,104.338,102.618,101.973
   with open(tsv_input,'r') as tsv:
      for line in tsv:
         linesp = line.split('\t')
         chr, read_pos, read_kmer, read_name, x, read_ind, event_current, event_sd, y, ref_kmer, model_current, ref_sd, z, event_levels = line.split('\t')
         if firstline:
            firstline = False
            continue
         else:
            if read_name != last_read:
               first_read_ind = int(read_ind)
            if ref_kmer == 'NNNNNN':
               continue
            if (read_name != last_read and read_kmer == ref_kmer) or (read_name == last_read and int(read_ind) > first_read_ind): #takes into account complementary palindromes
               rev = False
               meth_ref = meth_fwd
            else:
               rev = True
               meth_ref = meth_rev
            read_pos = int(read_pos)
            reference_kmer = meth_ref[read_pos:read_pos+k]
            #print read_name,rev,ref_kmer, reference_kmer

            #if pass the context for the previous (potentially) modified position, save and reset 
            if mpos and ((read_pos >= mpos+1 and read_name == last_read) or (read_name != last_read)):
        
               #write to file
               num_skips = len([x for x in diff_col if x == []])
               if num_skips <= 2: #accept max of 2 skips for 6mer contexts
                  if num_skips> 0:
                     w_skips.add((last_read,mpos))
                  with open(tsv_output,'a') as outfi:
                     diffs = [np.mean(kmer_pos) if kmer_pos!=[] else 0 for kmer_pos in diff_col]
                     if not last_rev:
                        diffs = diffs[::-1]
                     outfi.write(last_read+'\t'+str(mpos)+'\t'+revcomp(meth_ref[mpos-k+1:mpos+k],last_rev)+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev)+'\t'+label+'\n')
                     last_info = last_read+'\t'+str(mpos)+'\t'+revcomp(meth_ref[mpos-k+1:mpos+k],last_rev)+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(last_rev)+'\t'+label
                  num_observations += 1
                  pos_set.add(last_pos)
                  read_set.add(last_read)
                  if len(read_set)%1000 == 0 and len(read_set) > last_read_num:
                     print len(read_set), 'reads examined'
                     last_read_num = len(read_set)
               else:
                   skipped_skips.add((last_read,mpos))
          
               #reset variables
               if len(reference_kmer.split('M')) < 2 or read_name != last_read or read_pos > mpos+k: #allow no more than k skips
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
                      print last_info
                      print reference_kmer,mpos,read_pos,read_pos>mpos,read_name,last_read,diff_col,mspacing
                      diff_col = [[] for i in range(k)]
                   #break
 
            #if modified base in reference, save surrounding context to call that position
            if len([x for x in reference_kmer if x == 'M']) >= 1:
               pos_in_kmer = [i for i,x in enumerate(list(reference_kmer)) if x == 'M'][0]
               #if new read, reset differences variable and proceed
               if mpos and read_name != last_read:
                  mpos = None
                  diff_col = [[] for i in range(k)]
               #if new read or new position
               if not mpos: #TODO: reject any positions too close to beginning of read
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

   print num_observations,'observations'
   num_pos = len(pos_set)
   print num_pos,'positions'
   #cPickle.dump(diffs_by_read,open(tsv_input+'.pkl','wb'))
   print len(multi_meth_pos_set),'regions with multiple methylated bases'
   print len(w_skips), 'observations with skips'
   print len(skipped_skips), 'observations with too many skips'

