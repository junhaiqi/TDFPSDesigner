# Select barcodes from the entire k-mer space:
python selectBarcodeSeq.py --length 10 \
        --qsize 10000 \
        --outdir test_select_kmer \
        --threshold 10 \
        --thread-num 8 \
        --mode kmer \
        --seed 15 \
        --training-precison-cutoff 0.95 \
        --kit dna-r10-min

# Select barcodes from the given fasta file containing DNA sequences:
python selectBarcodeSeq.py --length 10 \
        --qsize 10000 \
        --outdir test_select_fasta \
        --threshold 10 \
        --thread-num 8 \
        --mode fasta \
        --fasta test_select_kmer/first_selected_barcodes.fa \
        --seed 15 \
        --training-precison-cutoff 0.95 \
        --kit dna-r10-min
