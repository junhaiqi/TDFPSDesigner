
import argparse
from time import time
import os
import sys
from functools import partial
from multiprocessing import Pool
import warnings
import numpy as np
warnings.filterwarnings("ignore", category=Warning)

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


def generateAdapterSignal(AdapterFastaFile='test.fasta', outSignalDir='AdapterSignal', threadNum=2):
    """>0: top adapter, >1: tail adapter."""
    if not os.path.exists(outSignalDir):
        os.makedirs(outSignalDir)

    seqList = getSeqList(file_name=AdapterFastaFile)
    seqIDList = getIDList(file_name=AdapterFastaFile)
    SeqLength = len(seqList[0])

    seqTupleList = list(zip(seqList, seqIDList))
    generateTrueNanoporeSignal(
        seqTupleList=seqTupleList,
        output_folder=outSignalDir, sigroot='timeSeries', threadNum=2)
    
    return SeqLength
    
def CalMNDistMatrix(timeSeriesDataDirM, timeSeriesDataDirN, distMatrixOutFilePrefixName):

    mainCommand = './bin/CalDTWDistMatrixMN -i_M %s -i_N %s -o tempoutput/%s.txt' \
        % (timeSeriesDataDirM, timeSeriesDataDirN, distMatrixOutFilePrefixName)

    os.system(mainCommand)

    finalDistList = []
    with open('tempoutput/' + distMatrixOutFilePrefixName + '.txt') as f:
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
                             BarcodeLength = 24, outBarcodeSigFile='testBarcode.txt'):

    estimateBarcodeLength = BarcodeLength * 10 + 70 # a super parameter.
    refSignal = get_signal_file(filetxt_path = signalFile)[0:1000]
    querySignalAdapterPath = outAdapterSignalDir + '/' + 'timeSeries_0.txt'
    
    querySignalAdapter = get_signal_file(filetxt_path = querySignalAdapterPath)
    position_start = fromLongRefFindShortQuery(refSignal, querySignalAdapter)[1]
    barcodeSig = refSignal[position_start + 40: position_start + 40 + estimateBarcodeLength]

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

def demultiplexingByDistMatrix(DistMatrix=[[1, 2, 3, 4], [2, 1, 7, 9], [1, 2, 3, 4], [2, 1, 7, 9]]):

    DistMatrix = [list(item) for item in np.transpose(DistMatrix)]
    resList = [row.index(min(row)) for row in DistMatrix]

    return resList

def main(AdapterFastaFile='testData/testAdapter.fasta', \
         outAdapterSignalDir='testData/AdapterSignal', \
         barcodeFastaFile = 'testData/testAdapter.fasta', \
         outBarcodeSignalDir='testData/BarcodeSignal', \
         sequencedNanoporeSignalDir='testData/AdapterSignal', \
         outDemultiplexingResTxt = 'testData/testDemultiplexingRes.txt', \
         outDistMatrixFilePrefixName = 'test', 
         FlankSeqLength = 8, 
         threadNum = 32):
    
    # generate adapter signal.
    main_start_time = time()
    print("Start Demultiplexing...")
    generateAdapterSignal(AdapterFastaFile=AdapterFastaFile, outSignalDir=outAdapterSignalDir)
    # generate barcode signal.
    barcodelength = generateAdapterSignal(AdapterFastaFile=barcodeFastaFile, outSignalDir=outBarcodeSignalDir, threadNum=16)
    if sequencedNanoporeSignalDir[-1] == '/':
        sequencedNanoporeSignalDir = sequencedNanoporeSignalDir.strip('/')
    newDirName = '%sParsedBarcodeSig'%sequencedNanoporeSignalDir
    os.system('mkdir %s'%newDirName)
    
    decodeNum = len(os.listdir(sequencedNanoporeSignalDir))
    
    args1 = [sequencedNanoporeSignalDir + '/' + item
            for item in os.listdir(sequencedNanoporeSignalDir)]
    args2 = [outAdapterSignalDir] * decodeNum
    args3 = [barcodelength - 2*FlankSeqLength] * decodeNum
    args4 = [newDirName + '/' + 'timeSeries_%d.txt'%i for i in range(decodeNum)]
    
    args = [(args1[i], args2[i], args3[i], args4[i]) for i in range(decodeNum)]

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

    file = open(outDemultiplexingResTxt, 'w')
    if len(os.listdir(outBarcodeSignalDir))*len(os.listdir(sequencedNanoporeSignalDir)) < 1100000000:
        signalList = os.listdir(sequencedNanoporeSignalDir)
        
        print("##########Calculating DTW distance matrix##########")
        distMatrix = CalMNDistMatrix(timeSeriesDataDirM = outBarcodeSignalDir, 
                        timeSeriesDataDirN = newDirName, 
                        distMatrixOutFilePrefixName = outDistMatrixFilePrefixName)
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

    parser = argparse.ArgumentParser()
    
    parser.add_argument('-iAF', type=str, required=True,
                        help='Specify a input fasta file, which contains a adpter sequence.')

    parser.add_argument('-oASD', type=str, required=True,
                        help='Specify an output folder to load adapter signals.')

    parser.add_argument('-iBF', type=str, required=True,
                        help='Specify an input fasta file, which contains barcode sequences.')
    
    parser.add_argument('-oBSD', type=str, required=True,
                        help='Specify an output folder to load Noiseless barcode signals.')
    
    parser.add_argument('-iNS', type=str, required=True,
                        help='Specify an input folder, which contains nanopore signals to be demultiplexed.')

    parser.add_argument('-oRes', type=str, required=True,
                        help='Specifies an output file, which contains the results of the demultiplexing.')

    parser.add_argument('-iFL', type=int, required=False, default=8,
                        help='Specify the length of the flanking sequence in the barcode sequence. \
                        It needs to be specified when the length of the top flanking sequence is \
                        the same as the length of the tail flanking sequence for better detection \
                        of barcode fragment in nanopore signal.')
    
    parser.add_argument('-t', type=int, required=False, default=32,
                        help='Specifies the number of threads, which affects the speed of extracting barcde signals.')
    
    parser.add_argument('-oMN', type=str, required=False, default='DTWDistMatrix',
                        help='Specifies the prefix name of an output file, which contains a DTW distance matrix.')

    args = parser.parse_args()

    return args

    
def argsMain():
    args = getParameters()
    
    main(AdapterFastaFile = args.iAF, \
         outAdapterSignalDir = args.oASD, \
         barcodeFastaFile = args.iBF, \
         outBarcodeSignalDir = args.oBSD, \
         sequencedNanoporeSignalDir = args.iNS, \
         outDemultiplexingResTxt = args.oRes, \
         outDistMatrixFilePrefixName = args.oMN,
         FlankSeqLength = args.iFL,
         threadNum = args.t)

if __name__ == "__main__":
    argsMain()

    

    
    

