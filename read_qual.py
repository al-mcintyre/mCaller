from Bio import SeqIO
import numpy as np
import sys

def extract_read_quality(fastqfi):
   read2qual = {}
   for read in SeqIO.parse(fastqfi,"fastq"):
      read2qual[read.id] = np.mean(read.letter_annotations["phred_quality"])
   return read2qual
