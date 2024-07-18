
# python selectBarcodeSeq.py --length 20 \
#         --qsize 10000 \
#         --outdir 20_795_final_select \
#         --threshold 70 \
#         --thread-num 8 \
#         --fasta 20_795_barcodes.fa \
#         --training-recall-cutoff 0.98 \
#         --training-f1Score-cutoff 0.98 \
#         --mode fasta \
#         --kit dna-r9-min

# python selectBarcodeSeq.py --length 30 \
#         --qsize 10000 \
#         --outdir 30_2120_final_select \
#         --threshold 70 \
#         --thread-num 8 \
#         --training-num-each-barcode 100 \
#         --fasta 30_2120_barcodes.fa \
#         --mode fasta \
#         --kit dna-r9-min

# python selectBarcodeSeq.py --length 24 \
#         --qsize 10000 \
#         --outdir 24_1093_second_select_1 \
#         --threshold 50 \
#         --thread-num 8 \
#         --fasta 24_1093_barcodes.fa \
#         --mode fasta \
#         --kit dna-r9-min
