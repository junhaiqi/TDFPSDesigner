# This script to simulate nanopore signals by 'squigulator' (https://github.com/hasindu2008/squigulator). 

import sys
import os
from ont_fast5_api.fast5_interface import get_fast5_file
from matplotlib import pyplot as plt
import argparse

def sig2text(sigList, outFile):

    file = open(outFile, 'w')
    for sig in sigList:
        file.write('%f\n'%sig)
    file.close()

def text2sig(sigFile):

    file = open(sigFile)
    sigList = []
    for sig in file:
        sig = sig.strip('\n')
        if sig != '':
            sigList.append(float(sig))
    file.close()

    return sigList


def drawSigPlot(sigFile):
    sig = text2sig(sigFile)
    plt.plot(sig)
    plt.show()


def f2t(fast5Filepath, outSigsDir):  # This can be a single- or multi-read file, transfer fast5s to text files.

    if not os.path.exists(outSigsDir):
        os.makedirs(outSigsDir)

    with get_fast5_file(fast5Filepath, mode = "r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            readName = read.read_id.split('!')[1]
            sigFile = os.path.join(outSigsDir, '%s.txt'%readName)
            raw_data = list(raw_data)
            sig2text(raw_data, sigFile)

def simuSigs(fasta, 
             outFile, 
             mode = 'R9min', 
             ideal = False,
             ideal_amp = False,
             ideal_time = False,
             slow5_dir = 'tempoutput',
             fast5 = False):

    if '.fast' not in outFile:
        print('Out file must include .fast')
        exit(-1)

    prefixName = outFile.split('/')[-1].split('.')[0]

    if ideal_amp == False and ideal_time == False:
        if ideal:
            cmd = './bin/squigulator -x %s %s -o %s/%s.slow5 --full-contigs --ideal'%(mode, fasta, slow5_dir, prefixName)
        else:
            cmd = './bin/squigulator -x %s %s -o %s/%s.slow5 --full-contigs'%(mode, fasta, slow5_dir, prefixName)
    
    else:
        if ideal_amp and ideal_time:
            cmd = './bin/squigulator -x %s %s -o %s/%s.slow5 --full-contigs --ideal'%(mode, fasta, slow5_dir, prefixName)
    
        elif ideal_amp and ideal_time == False:
            cmd = './bin/squigulator -x %s %s -o %s/%s.slow5 --full-contigs --ideal-amp'%(mode, fasta, slow5_dir, prefixName)

        else:
            cmd = './bin/squigulator -x %s %s -o %s/%s.slow5 --full-contigs --ideal-time'%(mode, fasta, slow5_dir, prefixName)
            
    os.system(cmd)
    
    if fast5:
        cmd2 = './bin/slow5tools s2f %s/%s.slow5 -o %s' % (slow5_dir, prefixName, outFile)

        os.system(cmd2)

    # os.system('rm tempoutput/%s.slow5' % prefixName)

def initializationParameters():
    parser = argparse.ArgumentParser()

    parser.add_argument('--fasta', type = str, required = True,
                        help = 'Indicates a input fasta file that containing reference DNA sequences.')
    parser.add_argument('--fast5', type = str, required = True,
                        help = 'Indicates a output fast5 file that containing simulated nanopore signals corresponding to reference DNA sequences.')
    parser.add_argument('--kit', type = str, required = True, choices = ["dna-r9-min", "dna-r9-prom", "rna-r9-min", "rna-r9-prom", "dna-r10-min", "dna-r10-prom"],
                        help = 'Indicates the nanopore sequencing kits.')
    parser.add_argument('--ideal', type = bool, required = False, default = False, choices = [False, True],
                        help = 'Generate ideal signals with no noise. default: False')
    parser.add_argument('--ideal-amp', type = bool, required = False, default = False, choices = [False, True],
                        help = 'Generate signals with no amplitiude domain noise. default: False')
    parser.add_argument('--ideal-time', type = bool, required = False, default = False, choices = [False, True],
                        help = 'Generate signals with no time domain noise. default: False')
    parser.add_argument('--txt', type = bool, required = False, default = True, choices = [False, True],
                        help = 'Indicates the simulated signals with "txt" format will be output. default: True')
    parser.add_argument('--txt-dir', type = str, required = False, default = "outDir",
                        help = 'Indicates a floder that load the simulated signal files with "txt" format. default: outDir')

    args = parser.parse_args()
    return args

def main():
    args = initializationParameters()

    simuSigs(fasta = args.fasta, 
             outFile = args.fast5, 
             mode = args.kit,
             ideal = args.ideal,
             ideal_amp = args.ideal_amp,
             ideal_time = args.ideal_time)

    if args.txt:
        f2t(fast5Filepath = args.fast5, outSigsDir = args.txt_dir)

if __name__ == "__main__":
    main()

