###### &nbsp;&nbsp;&nbsp;&nbsp;H  
###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|  
###### H-C-aller  
###### &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;|  
###### &nbsp;&nbsp;&nbsp;&nbsp;H

This program is designed to call m6A from nanopore data using the differences between measured and expected currents.  

## Requirements
> - nanopolish (https://github.com/jts/nanopolish)
> - an aligner to create a bam file (has been tested with graphmap and bwa mem)
python packages
> - scikit-learn 
> - h5py
> - biopython 
> - matplotlib
> - seaborn
> - numpy
> - pysam
> - scipy
> - pandas

other
> - nanopore sequencing data (fastq format + fast5 to run nanopolish, basecalled using Albacore or another basecaller that saves event data)
> - a reference sequence file (fasta)

## Options
```
usage: mCaller_nanopolish.py [-h] (-p POSITIONS | -m MOTIF) -r REFERENCE -e
                             TSV -f FASTQ [-t THREADS] [-b BASE]
                             [-n NUM_VARIABLES] [--train] [-d MODELFILE]
                             [-s SKIP_THRESH] [-q QUAL_THRESH] [-c CLASSIFIER]
                             [-v]
```

optional arguments:
```
  -h, --help            show this help message and exit
  -p, --positions
                        file with a list of positions at which to classify
                        bases (must be formatted as space- or tab-separated
                        file with chromosome, position, strand, and label if
                        training)
  -m, --motif 
                        classify every base of type --base in the motif
                        specified instead (can be single one-mer)
  -r, --reference 
                        fasta file with reference aligned to
  -e, --tsv     tsv file with nanopolish event alignment
  -f, --fastq 
                        fastq file with nanopore reads
  -t, --threads
                        specify number of processes (default = 1) 
			multiple threads should work but output won't be sorted
  -b, --base  bases to classify as methylated or unmethylated (A or
                        C, default A)
  -n, --num_variables
                        change the length of the context used to classify
                        (default of 6 variables corresponds to 11-mer context
                        (6*2-1))
  --train               train a new model (requires labels in positions file)
  -d, --modelfile 
                        model file name
  -s, --skip_thresh 
                        number of skips to allow within an observation
                        (default 0)
  -q, --qual_thresh
                        quality threshold for reads (under development, please
                        sort your own reads for now)
  -c, --classifier 
                        use alternative classifier: options = NN (default) RF,
                        LR, or NBC
  -v, --version         print version
```

## Pipeline for methylation detection from R9 data

1. extract template strand reads from fast5 files using a method that saves the file path in the fastq header, eg.
``` 
nanopolish extract -q -t template <fast5 directory> > <filename>.fastq 
```
   or 
``` 
poretools fastq --type fwd <fast5 directory> > <filename>.fastq 
```
  for albacore version > 2.0, follow the nanopolish guidelines (https://github.com/jts/nanopolish) and use 
```
nanopolish index -d <fast5 directory> -q <filename>.fastq
```
2. align fastq reads to reference assembly (we have used both GraphMap and bwa mem, with comparable results):
``` 
bwa index <reference>.fasta 
bwa mem -x ont2d -t <num_threads> <reference>.fasta <filename>.fastq | samtools view -Sb - | samtools sort -T /tmp/<filename>.sorted -o <filename>.sorted.bam 
samtools index <filename>.sorted.bam 
```
   or 
``` 
graphmap align -r <reference>.fasta -d <filename>.fastq -o <filename>.sam 
samtools view -bS <filename>.sam | samtools sort -T /tmp/<filename>.sorted -o <filename>.sorted.bam 
samtools index <filename>.sorted.bam 
``` 
4. run nanopolish with the following command to save a tsv file and the event values scaled towards the model:
``` 
nanopolish eventalign -t <num_threads> --scale-events -n -r <filename>.fastq -b <filename>.sorted.bam -g <reference>.fasta > <filename>.eventalign.tsv
```
5. run mCaller to detect m6A:
```
mCaller_nanopolish.py <-m GATC or -p positions.txt> -r <reference>.fasta -e <filename>.eventalign.tsv -f <filename>.fastq -b A 
```
   This returns a tabbed file with chromosome, read name, genomic position, position k-mer context, features, strand, and label
6. (optionally) run summary script to generate a bed file of methylated positions:
```
make_bed.py -f <filename>.eventalign.diffs.6 -d 15 -m 0.5 
```

Results and analysis scripts for the E. coli datasets are provided in the bioRxiv folder. 

## Test data

Reference fasta, PacBio calls for m6A and a subset of A positions, and eventalign tsv + fastq are provided for a single read for testing purposes in the "testdata" folder. 

1. To run mCaller on the testdata, use:
``` 
./mCaller_nanopolish.py -p testdata/test_positions_<m6A/A>.txt -r testdata/pb_ecoli_polished_assembly.fasta -e testdata/masonread1.eventalign.tsv -d mod_NN_6_m6A.pkl -f testdata/masonread1.fastq 
```
   Can also try using -m GATC, although not all GATC positions within the read were identified as methylated on both strands using PacBio and the model is slightly weighted to accept more false negatives than false positives at the moment. 

  This will generate the output file testdata/masonread1.eventalign.diffs.6

``` 
./make_bed.py -f testdata/masonread1.eventalign.diffs.6 -d 1 -m 0.5 
```
  Will then generate the output bed file testdata/masonread.methylation.summary.bed with columns chrom, chromStart, chromEnd, context, % methylated, strand, depth of coverage

2. To train on the testdata (don't actually use this model trained on one read), try:
``` 
./mCaller_nanopolish.py -p testdata/test_positions.txt -r testdata/pb_ecoli_polished_assembly.fasta -e testdata/masonread1.eventalign.tsv -t 4 --train -f testdata/masonread1.fastq
```

  This will generate the output file model_NN_6_m6A.pkl
