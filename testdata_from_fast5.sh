filename='testdata/masonread1' 
ref='testdata/pb_ecoli_polished_assembly'
model='mod_NN_6_m6A.pkl'
#for R9.5 data test, use filename='testdata/r95_test_read' and model='r95_model_NN_6_m6A.pkl'
threads=1
nanopolish extract -q -t template $filename.fast5 -o $filename.fastq
graphmap align -r $ref.fasta -d $filename.fastq -o $filename.sam
samtools view -bS $filename.sam | samtools sort -T /tmp/$filename.sorted -o $filename.sorted.bam
rm $filename.sam
samtools index $filename.sorted.bam
nanopolish eventalign -t $threads --scale-events -n -r $filename.fastq -b $filename.sorted.bam -g $ref.fasta > $filename.eventalign.tsv
#mCaller_nanopolish.py -m GATC -b A -r $ref.fasta -e $filename.eventalign.tsv -d $model -f $filename.fastq -t $threads
mCaller_nanopolish.py -p testdata/test_positions_m6A.txt -b A -r $ref.fasta -e $filename.eventalign.tsv -d $model -f $filename.fastq -t $threads
make_bed.py -f $filename.eventalign.diffs.6 -d 1 -t 0.5 --plot
