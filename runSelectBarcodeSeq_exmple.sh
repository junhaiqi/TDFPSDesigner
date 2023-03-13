# Select barcodes from the entire k-mer space:
python selectBarcodeSeq.py -l 10 \
        -q 10000 \
        -o test_kmer_mode.txt \
        -oinfo tempoutput/kmer_log.txt \
        -d 10 \
        -t 8 \
        -m kmer \
        -s 15

# Select barcodes from the given fasta file containing DNA sequences:
python selectBarcodeSeq.py -l 10 \
        -q 10000 \
        -o test_Fasta_mode.txt \
        -oinfo tempoutput/fasta_log.txt \
        -d 10 \
        -t 8 \
        -f 10mer_filter_results.fasta \
        -m fasta