from Bio import SeqIO
import numpy as np
import sys
import gzip

def extract_read_quality(fastqfi):
   read2qual = {}
   if fastqfi.find(".gz")!=-1:
      with gzip.open(fastqfi, "rt") as handle:
         for read in SeqIO.parse(handle, "fastq"):
            rid = read.id
            rid = rid.split(':')[0].split('_')[0]
            read2qual[rid] = np.mean(read.letter_annotations["phred_quality"])
   else:
      for read in SeqIO.parse(fastqfi,"fastq"):
         rid = read.id
         rid = rid.split(':')[0].split('_')[0]
         read2qual[rid] = np.mean(read.letter_annotations["phred_quality"])
   return read2qual
