from Bio import SeqIO
import numpy as np
import sys
import gzip

def extract_read_quality(fastqfi):
   read2qual = {}
   if fastqfi.find(".gz")!=-1:
      with gzip.open(fastqfi, "rt") as handle:
         for read in SeqIO.parse(handle, "fastq"):
            read2qual[read.id] = np.mean(read.letter_annotations["phred_quality"])
   else:
      for read in SeqIO.parse(fastqfi,"fastq"):
         read2qual[read.id] = np.mean(read.letter_annotations["phred_quality"])
   return read2qual
