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
def extract_features(tsv_input,meth_fwd,meth_rev,label):
   k = 15 #normally 6
   last_kmer = ''
   last_read = ''
   kmer_level = []
   model_level = 0
   last_pos = 0
   read_col = []
   pos_col = []
   kmer_col = []
   diff_col = []
   last_pos_meth = False
   num_reads = 0
   pos_set = set()
   read_set = set()
   firstline = True
   last_read_num = 0
   multi_meth_pos_set = set()

   #only save one set of adjoining methylated positions at a time - once the set complete, write the positions to a file 
   with open(tsv_input,'r') as tsv:
      for line in tsv:
         linesp = line.split('\t')
         fwd,rev = False,False
         if firstline:
            firstline = False
            continue
         else:
            if linesp[2] == linesp[9]:
               rev = False
               adj = 0
               meth_ref = meth_fwd
            elif linesp[2] == revcomp(linesp[9]):
               rev = True
               adj = 0
               meth_ref = meth_rev
            elif linesp[9] == 'NNNNNN': #nothing changes, this measurement is skipped
               continue
            else:
               print line, linesp[2], revcomp(linesp[9])
            read_pos = int(linesp[1])
            #print linesp[9], meth_ref[read_pos:read_pos+6]
            #print meth_ref[read_pos:read_pos+6], linesp[2], linesp[9], revcomp(meth_ref[read_pos:read_pos+6])
            #break
            if len([x for x in meth_ref[read_pos:read_pos+k] if x == 'M']) > 1 and meth_ref[read_pos] == 'M':
               #print meth_ref[read_pos:read_pos+6], read_pos
               multi_meth_pos_set.add(read_pos)
            if len([x for x in meth_ref[read_pos+adj:read_pos+adj+k] if x == 'M']) == 1: #starting with cases where only a single position is methylated in a kmer
               #print meth_ref[read_pos+adj:read_pos+adj+6],linesp[9],meth_ref == meth_rev
               if linesp[9] == last_kmer:
                  kmer_level.append(float(linesp[6]))
               else:
                  if last_pos_meth:
                     read_col.append(last_read)
                     pos_col.append(last_pos)
                     kmer_col.append(last_kmer)
                     diff_col.append(np.mean(kmer_level)-model_level)
                  #if read_pos == last_pos + 1:
                  last_read = linesp[3]
                  model_level = float(linesp[10])
                  last_pos = read_pos
                  last_kmer = linesp[9]
                  kmer_level = [float(linesp[6])]
                  last_pos_meth = True
            elif last_pos_meth and 'M' not in set(list(meth_ref[read_pos:read_pos+k])):
               last_pos_meth = False
               read_col.append(last_read)
               pos_col.append(last_pos)
               kmer_col.append(last_kmer)
               diff_col.append(np.mean(kmer_level)-model_level)
               if len(set(pos_col)) == k and last_pos == pos_col[0]+k-1: #starting with cases without skips
                  with open('.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.'+label,'a') as outfi:
                     outfi.write(last_read+'\t'+str(last_pos)+'\t'+revcomp(meth_ref[last_pos-5:last_pos+k],rev)+'\t'+','.join([str(diff) for diff in diff_col])+'\t'+strand(rev)+'\t'+label+'\n')
                  num_reads += 1  
                  pos_set.add(last_pos)
                  read_set.add(last_read)
                  if len(read_set)%1000 == 0 and len(read_set) > last_read_num:
                     print len(read_set), 'reads examined'
                     last_read_num = len(read_set)
                  #print rev,revcomp(meth_ref[last_pos-5:last_pos+6],rev), last_pos, zip(pos_col,kmer_col,diff_col)
                  #if num_reads > 10:
                  #    break
               #for pos in pos_col:
               read_col, pos_col, kmer_col, diff_col = [],[],[],[]
            else:
               last_pos_meth = False
               read_col, pos_col, kmer_col, diff_col = [],[],[],[]
               last_kmer = ''
               meth_pos_set = set([])

   print num_reads,'reads represented'
   num_pos = len(pos_set)
   print num_pos,'positions represented'
   #cPickle.dump(diffs_by_read,open(tsv_input+'.pkl','wb'))
   print len(multi_meth_pos_set),'regions with multiple methylated bases'

