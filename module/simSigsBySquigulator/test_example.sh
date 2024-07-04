# Test ideal mode.
python simulate_nano_sigs.py \
       --fasta test/test1.fasta \
       --fast5 test/test.fast5 \
       --kit dna-r10-min \
       --ideal True \
       --txt-dir test/out_ideal

# Test ideal-amp mode.
python simulate_nano_sigs.py \
       --fasta test/test1.fasta \
       --fast5 test/test.fast5 \
       --kit dna-r10-min \
       --ideal-amp True \
       --txt-dir test/out_ideal_amp

# Test ideal-time mode.
python simulate_nano_sigs.py \
       --fasta test/test1.fasta \
       --fast5 test/test.fast5 \
       --kit dna-r10-min \
       --ideal-time True \
       --txt-dir test/out_ideal_time

# Test ideal-time and ideal-amp mode.
python simulate_nano_sigs.py \
       --fasta test/test1.fasta \
       --fast5 test/test.fast5 \
       --kit dna-r10-min \
       --ideal-time True \
       --ideal-amp True \
       --txt-dir test/out_ideal_time_amp

# Test default mode.
python simulate_nano_sigs.py \
       --fasta test/test1.fasta \
       --fast5 test/test.fast5 \
       --kit dna-r10-min \
       --txt-dir test/out_default