filename='testdata/r95_test_read' #mac178083_med_cornell_edu_20160922_FN_MN17623_sequencing_run_20160922_nasa_mix_98347_ch372_read1326_strand1'
ref='testdata/pb_ecoli_polished_assembly'
#../nanopolish_dir/nanopolish extract -q -t template $filename.fast5 > $filename.fastq
#graphmap align -r $ref.fasta -d $filename.fastq -o $filename.sam
#samtools view -bS $filename.sam | samtools sort -T /tmp/$filename.sorted -o $filename.sorted.bam
#samtools index $filename.sorted.bam
#../nanopolish_dir/nanopolish eventalign -t 1 --scale-events -n -r $filename.fastq -b $filename.sorted.bam -g $ref.fasta > $filename.eventalign.tsv
./mCaller_nanopolish.py -m GACC -b A -r $ref.fasta -e $filename.eventalign.tsv -d r95_model_NN_6_m6A.pkl -f $filename.fastq -t 2
#./mCaller_nanopolish.py -p testdata/test_positions_m6A.txt -r $ref.fasta -e $filename.eventalign.tsv -d mod_NN_6_m6A.pkl -f $filename.fastq -t 3

