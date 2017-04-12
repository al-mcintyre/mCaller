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
   last_kmer = ''
   last_read = ''
   kmer_level = []
   model_level = 0
   last_pos = 0
   last_pos_in_kmer = k
   mpos = None
   diff_col = []
   num_reads = 0
   w_skips = set()
   pos_set = set()
   read_set = set()
   firstline = True
   last_read_num = 0
   multi_meth_pos_set = set()
   tsv_output = 'mCaller_update.'+'.'.join(tsv_input.split('.')[:-1])+'.diffs.'+str(k)+'.'+label
   print tsv_output

   #only save one set of adjoining methylated positions at a time - once the set complete, write the positions to a file 
   #format: ecoli	805	CGCCAT	cc1da58e-3db3-4a4b-93c2-c78e1dbe6aba:1D_000:template	t	1	102.16	0.963	0.00175	CGCCAT	102.23	1.93	-0.03	101.973,100.037,102.403,101.758,104.338,102.618,101.973
   with open(tsv_input,'r') as tsv:
      for line in tsv:
         linesp = line.split('\t')
         chr, read_pos, read_kmer, read_name, base, x, event_current, event_sd, y, ref_kmer, model_current, ref_sd, z, event_levels = line.split('\t')
         fwd,rev = False,False
         if firstline:
            firstline = False
            continue
         else:
            if read_kmer == ref_kmer:
               rev = False
               adj = 0
               meth_ref = meth_fwd
            elif read_kmer == revcomp(ref_kmer):
               rev = True
               adj = 0
               meth_ref = meth_rev
            elif ref_kmer == 'NNNNNN': #nothing changes, this measurement is skipped
               continue
            else:
               print line, read_kmer, revcomp(ref_kmer),'this is new. what happening'
            read_pos = int(read_pos)
            reference_kmer = meth_ref[read_pos:read_pos+k]
 
            if len([x for x in reference_kmer if x == 'M']) >= 1:
               pos_in_kmer = [i for i,x in enumerate(list(reference_kmer)) if x == 'M'][0]
               if not mpos and pos_in_kmer >= (k-2): #disregard any methylated positions with > 1 skips at the beginning of a k-mer TODO: reject any methylated positions too close to beginning of read
                  mpos = read_pos+pos_in_kmer
                  last_pos_in_kmer = pos_in_kmer
                  last_read = read_name
                  kmer_level = [float(event_current)-float(model_current)]
                  if pos_in_kmer == k-2:
                     diff_col.append(0) #if one skip at beginning, code as 0  
                     w_skips.add(mpos)
               elif not mpos and pos_in_kmer < (k-2):
                  #print 'pos in read too low'
                  continue
               elif mpos and read_name != last_read: #skip position if read ends before context
                  #print 'read ends?'
                  mpos = None
                  continue               
               elif mpos:
                  if read_pos+pos_in_kmer <= mpos and pos_in_kmer == last_pos_in_kmer: #ref_kmer == last_kmer
                     kmer_level.append(float(event_current)-float(model_current))
                     #print 'stay',kmer_level
                  else:
                     if read_pos <= mpos+1:
                        diff_col.append(np.mean(kmer_level))   
                     last_read = read_name
                     last_pos = read_pos 
                     last_kmer = read_kmer
                     kmer_level = [float(event_current)-float(model_current)]
                     if pos_in_kmer == last_pos_in_kmer-2 and read_pos <= mpos: #add 0 for skip to list of differences
                        diff_col.append(0)
                        w_skips.add(mpos)
                     if pos_in_kmer < last_pos_in_kmer-2: #skip contexts that contain >1 skip for now
                        continue #but better to do later, otherwise multiple-mod regions don't work
                     last_pos_in_kmer = pos_in_kmer
               print mpos, reference_kmer, read_pos, pos_in_kmer, last_pos_in_kmer, diff_col, kmer_level

            #allow one skip at end of context, coded as 0
            if mpos and read_pos >= mpos+1:
               diff_col.append(np.mean(kmer_level))
               last_read = read_name
               last_pos = read_pos
               last_kmer = read_kmer
               kmer_level = [float(event_current)-float(model_current)]
               if read_pos >= mpos+2 and len(diff_col) < k:
                  #print 'last is skip'
                  diff_col.append(0)
                  last_pos_in_kmer = k-1
                  kmer_level = []
                  w_skips.add(mpos) #should save pos + read
               #print diff_col

            if mpos and read_pos >= mpos+1:
               print mpos, diff_col
               if len(diff_col) == k:
                  with open(tsv_output,'a') as outfi:
                     if rev:
                        diffs = diff_col[::-1]
                     else:
                        diffs = diff_col
                     outfi.write(last_read+'\t'+str(mpos)+'\t'+revcomp(meth_ref[mpos-k+1:mpos+k],rev)+'\t'+','.join([str(diff) for diff in diffs])+'\t'+strand(rev)+'\t'+label+'\n')
                  num_reads += 1
                  pos_set.add(last_pos)
                  read_set.add(last_read)
                  if len(read_set)%1000 == 0 and len(read_set) > last_read_num:
                      print len(read_set), 'reads examined'
                      last_read_num = len(read_set)
               else:
                  print 'next',read_pos,diff_col
                  last_pos_in_kmer = k
               if (len(reference_kmer.split('M')) < 2):
                  diff_col = []
                  last_kmer = ''
                  meth_pos_set = set()
                  mpos = None
                  last_pos_in_kmer = k
               else:
                  multi_meth_pos_set.add((last_read,mpos))
                  #print reference_kmer, mpos, diff_col
                  last_mpos = mpos
                  pos_in_kmer = len(reference_kmer.split('M')[0]) 
                  mpos = read_pos + pos_in_kmer
                  mspacing = mpos - last_mpos
                  last_pos_in_kmer = pos_in_kmer #last_pos_in_kmer - mspacing
                  diff_col = diff_col[mspacing:]
                  print reference_kmer, mpos, diff_col
                  
   print num_reads,'instances represented for'
   num_pos = len(pos_set)
   print num_pos,'positions'
   #cPickle.dump(diffs_by_read,open(tsv_input+'.pkl','wb'))
   print len(multi_meth_pos_set),'regions with multiple methylated bases'
   print len(w_skips), 'positions with skips'

