
import argparse
from time import time
import os
import sys
from functools import partial
from multiprocessing import Pool
import warnings
import numpy as np
warnings.filterwarnings("ignore", category=Warning)
from ont_fast5_api.fast5_interface import get_fast5_file

if 'module' not in sys.path:
    sys.path.append('module')

if 'module/generate_nanoTruesig/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/')

if 'module/queryLocalSignal/module' not in sys.path:
    sys.path.append('module/queryLocalSignal/module')

if 'module/queryLocalSignal/' not in sys.path:
    sys.path.append('module/queryLocalSignal/')

if 'module/generate_nanoTruesig/module/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/module/')

if 'module/generate_nanoTruesig/model_data/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/model_data/')

if 'module/simSigsBySquigulator/' not in sys.path:
    sys.path.append('module/simSigsBySquigulator/')

from simulate_nano_sigs import simuSigs as simBySquigulator
    
from findLocalSignalPosition import fromLongRefFindShortQuery
from generatNoiselessSignal import sequence_to_true_signal

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

def generateTrueNanoporeSignal(seqTupleList=[('AAATTGGTTCGCCCCCCGGCCCGGC', i) for i in range(10000)],
                               output_folder='signal', sigroot='timeSeries', threadNum=2):

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    seq2signal = partial(
        sequence_to_true_signal, output_folder=output_folder, sigroot=sigroot)

    args = [seqTupleList[i]
            for i in range(len(seqTupleList))]

    pool = Pool(threadNum)
    res = list(pool.imap(seq2signal, args))
    pool.close()
    pool.join()

def f2t(fast5Filepath, outSigsDir, sigRoot = 'timeSeries'):  # This can be a single- or multi-read file, transfer fast5s to text files.

    def sig2text(sigList, outFile):
        file = open(outFile, 'w')
        for sig in sigList:
            file.write('%f\n'%sig)
        file.close()
        
    if not os.path.exists(outSigsDir):
        os.makedirs(outSigsDir)

    with get_fast5_file(fast5Filepath, mode = "r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            readName = read.read_id.split('!')[1]
            sigFile = os.path.join(outSigsDir, f'{sigRoot}_{readName}.txt')
            raw_data = list(raw_data)
            sig2text(raw_data, sigFile)

def squigulatorAPI(fastaFile, kit, outDir, sigRoot = 'timeSeries', fast5 = True):
    prefixName = fastaFile.split('/')[-1].split('.')[0]
    fast5Path = f'tempoutput/{prefixName}.fast5'
    simBySquigulator(fasta = fastaFile,
            outFile = fast5Path,
            mode = kit,
            ideal = False,
            ideal_amp = False,
            ideal_time = False,
            fast5 = fast5)
    f2t(fast5Filepath = fast5Path, outSigsDir = outDir, sigRoot = sigRoot)

def generateAdapterSignal(AdapterFastaFile='test.fasta', outSignalDir='AdapterSignal', kit = 'dna-r9-min', threadNum=2):
    """>0: top adapter, >1: tail adapter."""
    if not os.path.exists(outSignalDir):
        os.makedirs(outSignalDir)

    seqList = getSeqList(file_name=AdapterFastaFile)
    seqIDList = getIDList(file_name=AdapterFastaFile)
    SeqLength = len(seqList[0])

    seqTupleList = list(zip(seqList, seqIDList))
    if kit == 'dna-r9-min':
        generateTrueNanoporeSignal(
            seqTupleList=seqTupleList,
            output_folder=outSignalDir, sigroot='timeSeries', threadNum=2)
    else:
        squigulatorAPI(fastaFile = AdapterFastaFile, kit = kit, outDir = outSignalDir, fast5 = True)
    
    return SeqLength
    
def CalMNDistMatrix(timeSeriesDataDirM, timeSeriesDataDirN, distMatrixOutFile):

    mainCommand = './bin/CalDTWDistMatrixMN -i_M %s -i_N %s -o %s' \
        % (timeSeriesDataDirM, timeSeriesDataDirN, distMatrixOutFile)

    os.system(mainCommand)

    finalDistList = []
    with open(distMatrixOutFile) as f:
        lines = f.readlines()
        for line in lines:
            distValueList = [float(item) for item in line.strip(
                '\n').split(' ') if item != '']
            finalDistList.append(distValueList)

    return finalDistList

def get_signal_file(filetxt_path):

    signal_list = list()
    with open(filetxt_path, 'r') as f:
        signal_file = f.readlines()
        for line in signal_file:
            signal_Value = float(line.rstrip())
            signal_list.append(signal_Value)

    signal = np.array(signal_list)
    return signal

def findBarcodeSinalPosition(signalFile='signal.txt', outAdapterSignalDir='AdapterSignal', \
                             BarcodeLength = 24, outBarcodeSigFile='testBarcode.txt', kit = ''):

    estimateBarcodeLength = BarcodeLength * 10 + 70 # a super parameter.
    refSignal = get_signal_file(filetxt_path = signalFile)[0:1000]
    querySignalAdapterPath = outAdapterSignalDir + '/' + 'timeSeries_0.txt'
    
    querySignalAdapter = get_signal_file(filetxt_path = querySignalAdapterPath)
    position_start = fromLongRefFindShortQuery(refSignal, querySignalAdapter)[1]
    if kit == 'dna-r9-min' or kit == 'dna-r9-prom':
        barcodeSig = refSignal[position_start + 40: position_start + 40 + estimateBarcodeLength]
    elif kit == 'dna-r10-min' or kit == 'dna-r10-prom':
        barcodeSig = refSignal[position_start + 70: position_start + 70 + estimateBarcodeLength]
    with open(outBarcodeSigFile, 'w') as f:
        for item in barcodeSig:
            f.write('%s\n'%str(item))
    return barcodeSig

def splitSonSignalFromSignalFile(sigValueList, resFile):

    file = open(resFile, 'w')
    for sig in sigValueList:
        file.write(str(sig))
        file.write('\n')
    file.close()

def findOutliersBound(data, threshold = 5):
    """By 5-sigma method to find outliers."""
    mean = sum(data) / len(data)
    std_dev = (sum((x - mean) ** 2 for x in data) / len(data)) ** 0.5
    upper_bound = mean + threshold * std_dev
    return upper_bound

def demultiplexingByDistMatrix(DistMatrix=[[1, 2, 3, 4], [2, 1, 7, 9], [1, 2, 3, 4], [2, 1, 7, 9]]):

    DistMatrix = [list(item) for item in np.transpose(DistMatrix)]
    resList = [-1 for row in DistMatrix]  # '-1' indicates the barcode label is fuzzy. 
    distList = [min(row) for row in DistMatrix]
    outBound = findOutliersBound(distList)
    t = 0
    for row in DistMatrix:
        dist = distList[t]
        if dist <= outBound:
            resList[t] = row.index(dist)
        t += 1
    return resList

def main(AdapterFastaFile='testData/testAdapter.fasta', \
         barcodeFastaFile = 'testData/testAdapter.fasta', \
         sequencedNanoporeSignalDir='testData/AdapterSignal', \
         outDir = 'test', \
         FlankSeqLength = 8, \
         threadNum = 32, \
         kit = ''):
    
    # generate adapter signal.
    main_start_time = time()
    print("Start Demultiplexing...")
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    
    adapterSigDir = f'{outDir}/adapterSig'
    if not os.path.exists( adapterSigDir ):
        os.makedirs( adapterSigDir )
    generateAdapterSignal(AdapterFastaFile=AdapterFastaFile, outSignalDir=adapterSigDir, kit = kit)
    # generate barcode signal.
    outBarcodeSignalDir = f'{outDir}/barcodeTrueSig'
    if not os.path.exists( outBarcodeSignalDir ):
        os.makedirs( outBarcodeSignalDir )
    barcodelength = generateAdapterSignal(AdapterFastaFile=barcodeFastaFile, outSignalDir=outBarcodeSignalDir, threadNum=16, kit = kit)
    if sequencedNanoporeSignalDir[-1] == '/':
        sequencedNanoporeSignalDir = sequencedNanoporeSignalDir.strip('/')
    
    extractedBarcodeSigDir = f'{outDir}/extractedBarcodeSig'
    if not os.path.exists( extractedBarcodeSigDir ):
        os.makedirs( extractedBarcodeSigDir )
    
    decodeNum = len(os.listdir(sequencedNanoporeSignalDir))
    
    args1 = [sequencedNanoporeSignalDir + '/' + item
            for item in os.listdir(sequencedNanoporeSignalDir)]
    args2 = [adapterSigDir] * decodeNum
    args3 = [barcodelength - 2*FlankSeqLength] * decodeNum
    args4 = [extractedBarcodeSigDir + '/' + 'timeSeries_%d.txt'%i for i in range(decodeNum)]
    args5 = [kit for i in range(decodeNum)]
    
    args = [(args1[i], args2[i], args3[i], args4[i], args5[i]) for i in range(decodeNum)]

    start_time = time()
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    print("Extracting barcode signal from nanopore signal...")
    pool = Pool(threadNum)
    pool.starmap(findBarcodeSinalPosition, args)
    pool.close()
    pool.join()
    end_time = time()
    print('Extracting time: %fs'%(end_time - start_time))
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

    outDemultiplexingResTxt = f'{outDir}/dem_res.txt'
    outDistMatrixFile = f'{outDir}/dtw_dist_matrix.txt'
    file = open(outDemultiplexingResTxt, 'w')
    if len(os.listdir(outBarcodeSignalDir))*len(os.listdir(sequencedNanoporeSignalDir)) < 1100000000:
        signalList = os.listdir(sequencedNanoporeSignalDir)
        
        print("##########Calculating DTW distance matrix##########")
        distMatrix = CalMNDistMatrix(timeSeriesDataDirM = outBarcodeSignalDir, 
                        timeSeriesDataDirN = extractedBarcodeSigDir, 
                        distMatrixOutFile = outDistMatrixFile)
        print("###############Calculation completed###############")
        
        demultiplexingResList = demultiplexingByDistMatrix(DistMatrix=distMatrix)
        for i in range(len(signalList)):
            file.write('%s: %d\n'%(args1[i], demultiplexingResList[i]))
    else:
        print("The input data volume is too large, and does not support demultiplexing!")
        exit()
        
    file.close()
    
    main_end_time = time()
    print("End demultiplexing! Total run time: %fs"%(main_end_time - main_start_time))

def getParameters():

    parser = argparse.ArgumentParser('This script attempts to solve the demultiplexing problem in nanopore multi-sample sequencing.')
    
    parser.add_argument('--iAF', type=str, required=True,
                        help='Specify a input fasta file, which contains a adpter sequence.')

    parser.add_argument('--iBF', type=str, required=True,
                        help='Specify an input fasta file, which contains barcode sequences with flanking sequence.')
    
    parser.add_argument('--iNS', type=str, required=True,
                        help='Specify an input folder, which contains nanopore signals (.txt) to be demultiplexed.')

    parser.add_argument('--oRes', type=str, required=True,
                        help='Specifies an output folder, which contains the results of the demultiplexing.')

    parser.add_argument('--iFL', type=int, required=False, default=8,
                        help='Specify the length of the flanking sequence in the barcode sequence. \
                        It needs to be specified when the length of the top flanking sequence is \
                        the same as the length of the tail flanking sequence for better detection \
                        of barcode fragment in nanopore signal.')

    parser.add_argument('--kit', type=str, required=False, default='dna-r9-min', choices = ["dna-r9-min", "dna-r9-prom","dna-r10-min", "dna-r10-prom"], help='Specify ONT sequencing kit.')
    
    parser.add_argument('--thread-num', type=int, required=False, default=32,
                        help='Specifies the number of threads, which affects the speed of extracting barcde signals.')
    
    args = parser.parse_args()

    return args

    
def argsMain():
    args = getParameters()
    
    main(AdapterFastaFile = args.iAF, \
         barcodeFastaFile = args.iBF, \
         sequencedNanoporeSignalDir = args.iNS, \
         outDir = args.oRes, \
         FlankSeqLength = args.iFL,
         threadNum = args.thread_num,
         kit = args.kit)

if __name__ == "__main__":
    argsMain()

    

    
    

