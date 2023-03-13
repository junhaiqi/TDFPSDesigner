from time import time
from tqdm import tqdm
import os
import sys
from multiprocessing import Pool
from functools import partial
import warnings
warnings.filterwarnings("ignore", category=Warning)
import argparse

if 'module/initialSelection/' not in sys.path:
    sys.path.append('module/initialSelection/')

if 'module/generate_nanoTruesig/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/')

if 'module/generate_nanoTruesig/module/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/module/')

if 'module/generate_nanoTruesig/model_data/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/model_data/')

#####################################################################################################################
# Note: Given the DNA sequence, based on the 5-mer table, we can obtain its noise-free nanopore signal or simulated 
# nanopore signal containing noise. We can choose barcode from the noiseless signal space or the noisy signal space.

from generatNoiseSignal import sequence_to_true_signal # From noise signal space to select barcode signal.
# from generatNoiselessSignal import sequence_to_true_signal # From noiseless signal space to select barcode signal.
#####################################################################################################################

from make_sequence import mainFunction as initialSelection


def generateTrueNanoporeSignal(seqTupleList=[('AAATTGGTTCGCCCCCCGGCCCGGC', i) for i in range(10000)],
                               output_folder='signal', sigroot='timeSeries', threadNum=32):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    seq2signal = partial(
        sequence_to_true_signal, output_folder=output_folder, sigroot=sigroot)

    args = [seqTupleList[i]
            for i in range(len(seqTupleList))]

    pool = Pool(threadNum)
    res = list(tqdm(pool.imap(seq2signal, args)))
    pool.close()
    pool.join()


def getSeqList(file_name='24mer_filter_results.fasta'):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x != '', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    a = list(seq_list)
    return a


def getIDList(file_name='24mer_filter_results.fasta'):
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: int(x.split('|')[0][1:]), lines)
    b = list(id_list)
    return b

def fromFastaFile2Signal(fastaFilePath, threadNum):

    seqList = getSeqList(file_name=fastaFilePath)
    seqIDList = getIDList(file_name=fastaFilePath)

    output_folder = fastaFilePath.split('/')[-1].split('.')[0] + '_nanoporeSignal%d'%len(seqIDList)
    print('######%d noise nanopore signals are being generated######'% len(seqIDList))
    seqTupleList = list(zip(seqList, seqIDList))
    generate_start_time = time()
    if output_folder not in os.listdir('.'):
        generateTrueNanoporeSignal(
            seqTupleList=seqTupleList, output_folder=output_folder, sigroot='timeSeries', threadNum=threadNum)
    
        generate_end_time = time()
        print('######%d noise nanopore signals are generated! Total time: %fs######' %
            (len(seqIDList), generate_end_time - generate_start_time))
        print('\n')
    
    seqTupleDict = {}
    for item in seqTupleList:
        seqTupleDict[item[1]] = item[0]

    return seqTupleDict



def FromInitialSelectionToGetSinal(selectLength=24, selectQuantity=10000, randomSeed=0, threadNum=32):

    print('######Initial selection######')
    init_start_time = time()
    initialSelection(length=selectLength,
                     quantity=selectQuantity, seed=randomSeed)
    init_end_time = time()
    print('######End selection! Total time: %fs######' %
          (init_end_time - init_start_time))
    print('\n')

    initialSelectionFile = '%dmer_filter_results.fasta' % selectLength

    seqTupleDict = fromFastaFile2Signal(fastaFilePath = initialSelectionFile, threadNum = threadNum)

    return seqTupleDict
    
def byFPSCudaDTWFinalSelection(selectLength=24, selectQuantity=10000, randomSeed=0, \
                               slectedInfoFile='BarcodeInfo.txt', \
                               finalSelectedFile = 'finalSelectedBarcode.fasta',
                               initFastaFile = 'test.fasta',
                               mode = 'fasta', 
                               fastaFilePath='test.fasta',
                               thresFactor=0, threadNum=32):
    start_time = time()
    # print(mode)
    if mode != 'fasta':
        seqTupleDict = FromInitialSelectionToGetSinal(selectLength=selectLength, selectQuantity=selectQuantity, \
                                   randomSeed=randomSeed, threadNum=threadNum)
        fastaFilePath = '%dmer_filter_results'%selectLength
        output_folder = fastaFilePath + '_nanoporeSignal%d'%selectQuantity
    else:
        try:
            seqTupleDict = fromFastaFile2Signal(fastaFilePath = fastaFilePath, threadNum = threadNum)
            output_folder = fastaFilePath.split('/')[-1].split('.')[0] + '_nanoporeSignal%d'%selectQuantity
        except:
            print('%s is not exist! Pleace check it!'%fastaFilePath)
            return 0
    # print(seqTupleDict)
    
    print("######Final selection######")
    FPSCudaDTWCommand = './bin/FpsCudaDTWThreshold -i %s -l %d -o %s -t %f' \
                        %(output_folder, selectLength, slectedInfoFile, thresFactor)
    
    os.system(command=FPSCudaDTWCommand)
    
    with open(slectedInfoFile, 'r') as sif:
        lines = sif.readlines()
    selectedSeqIndexList = [int(item) for item in lines[-1].strip('\n').split(' ') if item != '']
    
    writeFile = open(finalSelectedFile, 'w')
    t = 0
    for _index in selectedSeqIndexList:
        writeFile.write('>%d\n'%t)
        writeFile.write(seqTupleDict[_index])
        writeFile.write('\n')
        t += 1
    writeFile.close()
    end_time = time()
    
    print('The total time for the entire barcode selection process is: %fs.'%(end_time - start_time))
    

def get_parameters():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', type=int, required=True,
                        help='Specify the length of the designed barcode.')
    
    parser.add_argument('-q', type=int, required=True,
                        help='Specify the size of the initially selected \
                        sequence space, which is recommended to be more than 100000.')
    
    parser.add_argument('-o', type=str, required=True,
                        help='Specify the output file, which contains the final barcode sequences.')

    parser.add_argument('-s', type=int, required=False, default=0,
                        help='Specify a random seed to determine the initially selected barcode signal, \
                        have a slight impact on the size of the final barcode set.')

    parser.add_argument('-oinfo', type=str, required=True,
                        help='Specify an output file to record output information.')
    
    parser.add_argument('-d', type=float, required=False, default=0, 
                        help='Specify a value to control the threshold of the TDFPS algorithm, the recommended value is 0~30.')
    
    parser.add_argument('-t', type=int, required=False, default=32,
                        help='Specify the number of threads.')

    parser.add_argument('-m', type=str, required=True, default='kmer',
                        help='Specify the selected mode. If the mode is "fasta", \
                        then -f must be followed by a file of "fasta" type.')

    parser.add_argument('-f', type=str, required=False, default=None,
                        help='Specify a file(format: fasta). The sequences \
                        contained in this file must be of the same length. \
                        When the "-m" is followed by "fasta", this parameter \
                        must be used. In other cases, it has no effect.')
    
    parser.add_argument('-introduction', type=str, required=False,
                        help='This script attempts to solve the barcode design problem in nanopore multi-sample sequencing.')
    
    args = parser.parse_args()
    
    return args

def main():
    
    args = get_parameters()
    
    byFPSCudaDTWFinalSelection(selectLength=args.l, selectQuantity=args.q, randomSeed=args.s, \
                               slectedInfoFile=args.oinfo, finalSelectedFile = args.o, \
                               thresFactor=args.d, threadNum=args.t, mode=args.m, fastaFilePath=args.f)

if __name__ == "__main__":
    main()
    
