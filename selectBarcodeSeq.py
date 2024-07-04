from time import time
from tqdm import tqdm
import os
import sys
from multiprocessing import Pool
from functools import partial
import warnings
warnings.filterwarnings("ignore", category=Warning)
import argparse
from ont_fast5_api.fast5_interface import get_fast5_file
import pyslow5

if 'module/initialSelection/' not in sys.path:
    sys.path.append('module/initialSelection/')

if 'module/generate_nanoTruesig/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/')

if 'module/generate_nanoTruesig/module/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/module/')

if 'module/generate_nanoTruesig/model_data/' not in sys.path:
    sys.path.append('module/generate_nanoTruesig/model_data/')

if 'module/simSigsBySquigulator/' not in sys.path:
    sys.path.append('module/simSigsBySquigulator/')

#####################################################################################################################
# Note: Given the DNA sequence, based on the 5-mer table, we can obtain its noise-free nanopore signal or simulated 
# nanopore signal containing noise. We can choose barcode from the noiseless signal space or the noisy signal space.

from generatNoiseSignal import sequence_to_true_signal, generate_noisy_sigs # From noise signal space to select barcode signal.
from generatNoiselessSignal import generate_true_sigs
# from generatNoiselessSignal import sequence_to_true_signal # From noiseless signal space to select barcode signal.
#####################################################################################################################

from make_sequence import mainFunction as initialSelection

from simulate_nano_sigs import simuSigs as simBySquigulator


def generateTrueNanoporeSignal(seqTupleList=[   ('AAATTGGTTCGCCCCCCGGCCCGGC', i) for i in range(10000)  ],
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

def sig2text(sigList, outFile):
    file = open(outFile, 'w')
    for sig in sigList:
        file.write('%f\n'%sig)
    file.close()

def f2t(fast5Filepath, outSigsDir, mode = None):  # This can be a single- or multi-read file, transfer fast5s to text files.
        
    if not os.path.exists(outSigsDir):
        os.makedirs(outSigsDir)

    with get_fast5_file(fast5Filepath, mode = "r") as f5:
        for read in f5.get_reads():
            raw_data = read.get_raw_data()
            readName = read.read_id.split('!')[1]
            if mode == 'root':
                sigFile = os.path.join(outSigsDir, 'timeSeries_%s.txt'%readName)
            else:
                sigFile = os.path.join(outSigsDir, '%s.txt'%readName)
            raw_data = list(raw_data)
            sig2text(raw_data, sigFile)


def squigulatorAPI(fastaFile, kit, outDir, mode = None, ideal = False, ideal_amp = False, ideal_time = False, slow5_dir = 'tempoutput', fast5 = True):
    prefixName = fastaFile.split('/')[-1].split('.')[0]
    fast5Path = f'{slow5_dir}/{prefixName}.fast5'
    simBySquigulator(fasta = fastaFile,
            outFile = fast5Path,
            mode = kit,
            ideal = ideal,
            ideal_amp = ideal_amp,
            ideal_time = ideal_time,
            slow5_dir = slow5_dir,
            fast5 = fast5)
    if fast5:
        f2t(fast5Filepath = fast5Path, outSigsDir = outDir, mode = mode)
   
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
        info = list(lines)
        if '_' in info[0]:
            id_list = map(lambda x: int(x.strip('\n').split('_')[1]), lines)
            b = list(id_list)
        else:
            b = []
            for item in info:
                idx = int( item.split('>')[-1] )
                b.append( idx )
        return b

def read_slow5(slow5_file_path, out_dir_path, thread_num, sig_root = 'timeSeries'):
    s5 = pyslow5.Open(slow5_file_path, 'r')
    if not os.path.exists( out_dir_path ):
        os.makedirs( out_dir_path )
    reads = s5.seq_reads_multi(threads = thread_num, batchsize = thread_num + 1)
    count = 0
    for read in reads:
        out_sig_path = f'{out_dir_path}/{sig_root}_{count}.txt'
        sig2text( read['signal'], out_sig_path )
        count += 1

def fromFastaFile2Signal(fastaFilePath, threadNum, kit = '', output_folder = 'test'):
    seqList = getSeqList(file_name=fastaFilePath)
    # seqIDList = getIDList(file_name=fastaFilePath)
    seqIDList = [i for i in range( len(seqList)) ]
    seqTupleList = list(zip(seqList, seqIDList))
    # output_folder = fastaFilePath.split('/')[-1].split('.')[0] + '_nanoporeSignal%d'%len(seqIDList)
    selectLength = len( seqList[0] )
    # output_folder = f'{output_folder}/{selectLength}mer_init_filter_results_nanopore_sigs'
    print('######%d noise nanopore signals are being generated######'% len(seqIDList))
    generate_start_time = time()

    if kit == 'dna-r9-min':
        if not os.path.exists( output_folder ):
            os.makedirs( output_folder )
        generateTrueNanoporeSignal(
            seqTupleList=seqTupleList, output_folder=output_folder, sigroot='timeSeries', threadNum=threadNum)
    else:
        # squigulatorAPI(fastaFile = fastaFilePath, kit = kit, outDir = output_folder)
        squigulatorAPI(fastaFile = fastaFilePath, kit = kit, outDir = output_folder, slow5_dir = 'tempoutput', fast5 = False)

        prefixName = fastaFilePath.split('/')[-1].split('.')[0]
        slow_path = f'tempoutput/{prefixName}.slow5'
        read_slow5(slow_path, output_folder, threadNum, sig_root = 'timeSeries')

    generate_end_time = time()
    print('######%d noise nanopore signals are generated! Total time: %fs######' %
        (len(seqIDList), generate_end_time - generate_start_time))
    print('\n')
    
    seqTupleDict = {}
    for item in seqTupleList:
        seqTupleDict[item[1]] = item[0]

    return seqTupleDict

def FromInitialSelectionToGetSinal(selectLength=24, selectQuantity=10000, randomSeed=0, threadNum=32, filter2 = False, kit = '', outDir = ''):

    print('######Initial selection######')
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    
    initialSelectionFile = f'{outDir}/{selectLength}mer_init_filter_results.fasta'
    output_folder = f'{outDir}/{selectLength}mer_init_filter_results_nanopore_sigs'
    init_start_time = time()
    initialSelection(length=selectLength,
                     quantity=selectQuantity, seed=randomSeed, filter2 = filter2, output_path = initialSelectionFile)
    init_end_time = time()
    print('######End selection! Total time: %fs######' %
          (init_end_time - init_start_time))
    print('\n')

    seqTupleDict = fromFastaFile2Signal(fastaFilePath = initialSelectionFile, threadNum = threadNum, kit = kit, output_folder = output_folder)

    return seqTupleDict
    
def byFPSCudaDTWFinalSelection(selectLength=24, selectQuantity=10000, randomSeed=0, \
                               outDir = 'test',
                               mode = 'fasta', 
                               fastaFilePath='test.fasta',
                               thresFactor=0, threadNum=32, filter2 = False, kit = ''):
    start_time = time()
    # print(mode)
    if not os.path.exists( outDir ):
        os.makedirs( outDir )
    if mode != 'fasta':
        # fastaFilePath = f'{outDir}/{selectLength}mer_init_filter_results'
        output_folder = f'{outDir}/{selectLength}mer_init_filter_results_nanopore_sigs'
        seqTupleDict = FromInitialSelectionToGetSinal(selectLength=selectLength, selectQuantity=selectQuantity, \
                                   randomSeed=randomSeed, threadNum=threadNum, filter2 = filter2, kit = kit, outDir = outDir)
    else:
        output_folder = f'{outDir}/{selectLength}mer_fasta_results_nanopore_sigs'
        seqTupleDict = fromFastaFile2Signal(fastaFilePath = fastaFilePath, threadNum = threadNum, output_folder = output_folder, kit = kit)
    
    print("######Final selection######")
    slectedInfoFile = f'{outDir}/TDFPS.info'
    FPSCudaDTWCommand = './bin/FpsCudaDTWThreshold -i %s -l %d -o %s -t %f' \
                        %(output_folder, selectLength, slectedInfoFile, thresFactor)
    
    os.system(command=FPSCudaDTWCommand)
    
    with open(slectedInfoFile, 'r') as sif:
        lines = sif.readlines()
        selectedSeqIndexList = [int(item) for item in lines[-1].strip('\n').split(' ') if item != '']
    
    firstSelectedFile = f'{outDir}/first_selected_barcodes.fa'
    writeFile = open(firstSelectedFile, 'w')
    t = 0
    for _index in selectedSeqIndexList:
        writeFile.write('>%d\n'%t)
        writeFile.write(seqTupleDict[_index])
        writeFile.write('\n')
        t += 1
    writeFile.close()
    end_time = time()
    
    print('The total time for the entire barcode selection process is: %fs.'%(end_time - start_time))
    
from sklearn.metrics import recall_score, precision_score, f1_score
import random
import numpy as np

from sklearn.metrics import recall_score, precision_score, f1_score
import random
import numpy as np

if 'module/queryLocalSignal/module' not in sys.path:
    sys.path.append('module/queryLocalSignal/module')
if 'module/queryLocalSignal/' not in sys.path:
    sys.path.append('module/queryLocalSignal/')

from findLocalSignalPosition import fromLongRefFindShortQuery

def get_signal_file(filetxt_path):

    signal_list = list()
    with open(filetxt_path, 'r') as f:
        signal_file = f.readlines()
        for line in signal_file:
            signal_Value = float(line.rstrip())
            signal_list.append(signal_Value)

    signal = np.array(signal_list)
    return signal

def find_barcode_region(raw_signal_file, query_signal_file, 
                        bar_len = 24, out_bar_sig = '', kit = 'dna-r9-min', region_cutoff = 1000):

    # estimateBarcodeLength = BarcodeLength * 10 + 70 # a super parameter.
    est_bar_len = bar_len * 10 + 70
    raw_signal = get_signal_file(filetxt_path = raw_signal_file)[0:region_cutoff]
    query_signal = get_signal_file(filetxt_path = query_signal_file)
    
    if kit == 'dna-r9-min' or kit == 'dna-r9-prom':
        position_start = fromLongRefFindShortQuery(raw_signal, query_signal)[1]
        queryed_signal = raw_signal[position_start + 40: position_start + 40 + est_bar_len]
        # queryed_signal = raw_signal[position_start: position_start + est_bar_len]

    elif kit == 'dna-r10-min' or kit == 'dna-r10-prom':
        position_start = fromLongRefFindShortQuery(raw_signal, query_signal)[1]
        queryed_signal = raw_signal[position_start + 70: position_start + 70 + est_bar_len]

    with open(out_bar_sig, 'w') as f:
        for item in queryed_signal:
            f.write('%s\n'%str(item))
    return queryed_signal

def muti_find(raw_signa_dir, query_signal_path = '', ex_out_dir = '', bar_len = 24, 
              thread_num = 16, sig_root = 'timeSeries', kit = 'dna-r9-min', region_cutoff = 1000):
    decode_num = len(os.listdir(raw_signa_dir))
    args1 = [raw_signa_dir + '/' + '%s_%d.txt'%(sig_root, i)
            for i in range(decode_num)]
    args2 = [query_signal_path] * decode_num
    args3 = [bar_len] * decode_num
    args4 = [ex_out_dir + '/' + '%s_%d.txt'%(sig_root, i) for i in range(decode_num)]
    args5 = [kit for i in range(decode_num)]
    args6 = [region_cutoff for i in range(decode_num)]
    
    args = [(args1[i], args2[i], args3[i], args4[i], args5[i], args6[i]) for i in range(decode_num)]

    start_time = time()
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")
    print("Extracting barcode signal from nanopore signal...")
    pool = Pool(thread_num)
    pool.starmap(find_barcode_region, args)
    pool.close()
    pool.join()
    end_time = time()
    print('Extracting time: %fs'%(end_time - start_time))
    print("&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&")

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

def get_dem_res(prased_signals_path = '',
                true_barcode_sigs_path = '',
                dist_matrix = '',
                dem_file = '',
                sig_root = 'timeSeries'):

    cmd = './bin/CalDTWDistMatrixMNSigroot -i_M %s -i_N %s -sig_root %s -o %s' \
            %(true_barcode_sigs_path, prased_signals_path, sig_root, dist_matrix)

    print(cmd)
    os.system(cmd)

    finalDistList = []
    with open(dist_matrix) as f:
        lines = f.readlines()
        for line in lines:
            distValueList = [float(item) for item in line.strip(
                '\n').split(' ') if item != '']
            finalDistList.append(distValueList)

    demultiplexingResList = demultiplexingByDistMatrix(DistMatrix=finalDistList)

    idxList = [item for item in range(len(os.listdir(prased_signals_path)))]

    file = open(dem_file, 'w')
    for i in range(len(idxList)):
        file.write('%s_%d.txt: %d\n'%(sig_root, idxList[i], demultiplexingResList[i]))

    file.close()

    return demultiplexingResList

class BarcodeSet:  # inital barcode set --> final barcode set
    def __init__(self, bar_fa_file, adapter_seq, top_flank_seq, bottom_flank_seq, obj_precison_cutoff = 0.99,
    obj_recall_cutoff = 0.99, obj_F1score_cutoff = 0.99, training_num_each_barcode = 500, out_dir = 'temp', source_type = 'dna-r9-min', thread_num = 8):
        self.bar_fa_file = bar_fa_file
        self.adapter_seq = adapter_seq
        self.top_flank_seq = top_flank_seq
        self.bottom_flank_seq = bottom_flank_seq
        self.obj_precison_cutoff = obj_precison_cutoff
        self.obj_recall_cutoff = obj_recall_cutoff
        self.obj_F1score_cutoff = obj_F1score_cutoff
        self.training_num_each_barcode = training_num_each_barcode
        self.out_dir = out_dir
        self.source_type = source_type
        self.out_fa_path = ''
        self.out_sig_dir_path = ''
        self.barcode_true_sigs_path = ''
        self.out_adapter_fa = ''
        self.adapter_sig_path = ''
        self.ex_out_dir = ''
        self.bar_fa_with_flanks = ''
        self.bar_len = 24
        self.barcode_num = 0
        self.dem_res = []
        self.thread_num = thread_num
        self.flanked_barcodes = []

    def read_fasta_sequences(self, file_path):
        sequences = []
        with open(file_path, 'r') as file:
            sequence = ''
            for line in file:
                if line.startswith('>'):
                    if sequence:
                        sequences.append(sequence)
                        sequence = ''
                else:
                    sequence += line.strip()
            if sequence:
                sequences.append(sequence)
        return sequences

    def sig2text(self, sigList, outFile):
        file = open(outFile, 'w')
        for sig in sigList:
            file.write('%f\n'%sig)
        file.close()

    def read_slow5(self, slow5_file_path, out_dir_path, sig_root = 'timeSeries'):
        s5 = pyslow5.Open(slow5_file_path, 'r')
        reads = s5.seq_reads_multi(threads = self.thread_num, batchsize = self.thread_num + 1)
        count = 0
        for read in reads:
            out_sig_path = f'{out_dir_path}/{sig_root}_{count}.txt'
            self.sig2text( read['signal'], out_sig_path )
            count += 1
    
    def generate_strand_library_seqs(self, read_len = 100):
        if not os.path.exists( self.out_dir ):
            os.makedirs( self.out_dir )
        bases = ['A', 'T', 'C', 'G']
        barcode_seqs = self.read_fasta_sequences( self.bar_fa_file )  # read short barcode sequences
        self.bar_fa_with_flanks = f'{self.out_dir}/barcode_with_flanks.fa'
        _file = open(self.bar_fa_with_flanks, 'w')
        
        t = 0
        for barcode in barcode_seqs:
            _file.write(f'>{t}\n')
            flank_barcode  = self.top_flank_seq + barcode + self.bottom_flank_seq
            _file.write(f'{flank_barcode}\n')
            self.flanked_barcodes.append( flank_barcode )
            t += 1
        _file.close()

        self.bar_len = len(barcode_seqs[0])
        self.barcode_num = len( barcode_seqs )
        out_fa_path = f'{ self.out_dir }/training.fa'
        self.out_fa_path = out_fa_path
        _file = open(out_fa_path, 'w')
        count = 0
        for i in tqdm(range( len(barcode_seqs) )):
            for j in range( self.training_num_each_barcode ):
                random_read = ''.join(random.choice(bases) for _ in range( read_len ))
                multiplexed_read = self.adapter_seq + self.top_flank_seq + barcode_seqs[i] + self.bottom_flank_seq + random_read
                _file.write(f'>{count}\n')
                _file.write(f'{multiplexed_read}\n')
                count += 1
        _file.close()

        self.out_adapter_fa = f'{ self.out_dir }/adapter.fa'
        _file = open(self.out_adapter_fa, 'w')
        _file.write('>adapter\n')
        _file.write(f'{ self.adapter_seq }\n')
        _file.close()

    def generate_sim_nano_sigs(self):  # include all strandard siganls of barcodes and noisy signals of multiplexed reads
        self.out_sig_dir_path = f'{self.out_dir}/all_sim_signals_{self.source_type}'
        self.barcode_true_sigs_path = f'{self.out_dir}/barcode_true_signals_{self.source_type}'
        if not os.path.exists( self.out_sig_dir_path ):
            os.makedirs( self.out_sig_dir_path )
        if not os.path.exists( self.barcode_true_sigs_path ):
            os.makedirs( self.barcode_true_sigs_path )
        
        squigulatorAPI(fastaFile = self.out_fa_path, kit = self.source_type, outDir = self.out_sig_dir_path, slow5_dir = self.out_dir, fast5 = False)
        self.read_slow5( f'{self.out_dir}/training.slow5', self.out_sig_dir_path)

        squigulatorAPI(fastaFile = self.bar_fa_with_flanks, kit = self.source_type, outDir = self.barcode_true_sigs_path, ideal = True, ideal_amp = True, ideal_time = True, mode = 'root', slow5_dir = self.out_dir, fast5 = False)

        prefix_name = self.bar_fa_with_flanks.split('/')[-1].split('.')[0]
        self.read_slow5( f'{self.out_dir}/{prefix_name}.slow5', self.barcode_true_sigs_path)

        squigulatorAPI(fastaFile = self.out_adapter_fa, kit = self.source_type, outDir = self.out_dir, slow5_dir = self.out_dir, ideal = True, ideal_amp = True, ideal_time = True, fast5 = True)
        self.adapter_sig_path = f'{ self.out_dir }/adapter.txt'

    def ex_barcode_sigs(self):
        self.ex_out_dir = self.out_dir + '/' + 'ex_bar_sigs'
        if not os.path.exists( self.ex_out_dir ):
            os.makedirs( self.ex_out_dir )
        muti_find(raw_signa_dir = self.out_sig_dir_path, query_signal_path = self.adapter_sig_path, 
        ex_out_dir = self.ex_out_dir, bar_len = self.bar_len, thread_num = self.thread_num, sig_root = 'timeSeries', 
        kit = self.source_type, region_cutoff = 1000)

    def dem_by_ex_sigs_info(self):
        dem_res = get_dem_res(prased_signals_path = self.ex_out_dir,
                    true_barcode_sigs_path = self.barcode_true_sigs_path,
                    dist_matrix = f'{self.out_dir}/dist_matrix.txt',
                    dem_file = f'{self.out_dir}/dem_res.txt',
                    sig_root = 'timeSeries')

        self.dem_res = dem_res

    def cal_eva_indicators(self):
        y_pred = self.dem_res
        y_true = [0 for i in range( len(y_pred) )]
        idx = 0
        for i in range( self.barcode_num ):
            for j in range( self.training_num_each_barcode ):
                y_true[ idx ] = i
                idx += 1
        precision_per_class = precision_score(y_true, y_pred, average=None, zero_division=0)
        recall_per_class = recall_score(y_true, y_pred, average=None, zero_division=0)
        f1_per_class = f1_score(y_true, y_pred, average=None, zero_division=0)

        each_class_precision_list = []
        each_class_recall_list = []
        each_class_f1score_list = []
        labels = sorted(set(y_true))
        for label, precision, recall, f1 in zip(labels, precision_per_class, recall_per_class, f1_per_class):
            each_class_precision_list.append( precision )
            each_class_recall_list.append( recall )
            each_class_f1score_list.append( f1 )

        return { 'recall_list': each_class_recall_list, 'precision_list': each_class_precision_list, 'f1_list': each_class_f1score_list }

    def filter_bad_barcodes(self):
        indicators = self.cal_eva_indicators()
        final_barcode_fa = f'{self.out_dir}/final_barcodes.fa'
        final_idx_file = f'{self.out_dir}/final_barcode_idx.txt'
        _file = open( final_barcode_fa, 'w' )
        _idx_file = open( final_idx_file, 'w' )
        _len = len(indicators[ 'recall_list' ])
        t = 0
        for i in range( _len ):
            if indicators[ 'recall_list' ][i] > self.obj_recall_cutoff:
                if indicators[ 'f1_list' ][i] > self.obj_F1score_cutoff and indicators[ 'precision_list' ][i] > self.obj_precison_cutoff:
                    _file.write(f'>{t}\n')
                    _idx_file.write( f'{i}\n' )
                    _file.write(f'{ self.flanked_barcodes[i] }\n')
                    t += 1
        _file.close()
        _idx_file.close()

    def run(self):
        self.generate_strand_library_seqs()
        self.generate_sim_nano_sigs()
        self.ex_barcode_sigs()
        self.dem_by_ex_sigs_info()
        self.filter_bad_barcodes()

def get_parameters():
    
    parser = argparse.ArgumentParser('This script attempts to solve the barcode design problem in nanopore multi-sample sequencing.')
    parser.add_argument('--length', type=int, required=True,
                        help='Specify the length of the designed barcode.')
    
    parser.add_argument('--qsize', type=int, required=True,
                        help='Specify the size of the initially selected \
                        sequence space, which is recommended to be more than 100000.')
    
    parser.add_argument('--outdir', type=str, required=True,
                        help='Specify the output file, which contains the final barcode sequences.')

    parser.add_argument('--seed', type=int, required=False, default=0,
                        help='Specify a random seed to determine the initially selected barcode signal, \
                        have a slight impact on the size of the final barcode set.')
    
    parser.add_argument('--threshold', type=float, required=False, default = 0, 
                        help='Specify a value to control the threshold of the TDFPS algorithm, the recommended value is 0~30.')
    
    parser.add_argument('--thread-num', type=int, required=False, default=32,
                        help='Specify the number of threads.')

    parser.add_argument('--mode', type=str, required=True, default='kmer', choices = ["kmer", "fasta"],
                        help='Specify the selected mode. If the mode is "fasta", \
                        then -f must be followed by a file of "fasta" type.')

    parser.add_argument('--fasta', type=str, required=False, default=None,
                        help='Specify a file(format: fasta). The sequences \
                        contained in this file must be of the same length. \
                        When the "--mode" is followed by "fasta", this parameter \
                        must be used. In other cases, it has no effect.')

    parser.add_argument('--kit', type=str, required=False, default='dna-r9-min', choices = ["dna-r9-min", "dna-r9-prom","dna-r10-min", "dna-r10-prom"], help='Specify ONT sequencing kit.')

    parser.add_argument('--adapter-seq', type=str, required=False, default='GGCGTCTGCTTGGGTGTTTAACCTTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT', help='Specify ONT adapter sequence for select barcodes again.')

    parser.add_argument('--top-flank-seq', type=str, required=False, default='AAGGTTAA', help='Specify ONT top flanking sequence for select barcodes again.')

    parser.add_argument('--bottom-flank-seq', type=str, required=False, default='CAGCACCT', help='Specify ONT bottom flanking sequence for select barcodes again.')

    parser.add_argument('--training-num-each-barcode', type=int, required=False, default=100, help='Specify the training number for select barcodes again.')

    parser.add_argument('--training-precison-cutoff', type=float, required=False, default=0.999, help='Specify training precison cutoff for select barcodes again.')

    parser.add_argument('--training-recall-cutoff', type=float, required=False, default=0.99, help='Specify training recall cutoff for select barcodes again.')

    parser.add_argument('--training-f1Score-cutoff', type=float, required=False, default=0.99, help='Specify training F1-Score cutoff for select barcodes again.')

    parser.add_argument('--bio-criteria', action='store_true',
                        help='Based on biological criteria, sequences with a GC content lower than 0.4 or greater than 0.6, sequences containing reapte triples, sequences containing GGC, and self-complementary sequences were filtered out.')

    args = parser.parse_args()
    
    return args
        
def main():
    
    args = get_parameters()
    
    byFPSCudaDTWFinalSelection(selectLength=args.length, selectQuantity=args.qsize, randomSeed=args.seed, \
                               outDir = args.outdir, \
                               thresFactor=args.threshold, threadNum=args.thread_num, mode=args.mode, fastaFilePath=args.fasta, \
                               filter2 = args.bio_criteria, kit = args.kit)

    firstBarcodesPath = f'{args.outdir}/first_selected_barcodes.fa'
    barcodeSet = BarcodeSet(firstBarcodesPath, args.adapter_seq,
    args.top_flank_seq, args.bottom_flank_seq, out_dir = args.outdir, training_num_each_barcode = args.training_num_each_barcode, thread_num = args.thread_num, obj_precison_cutoff = args.training_precison_cutoff, obj_recall_cutoff = args.training_recall_cutoff, obj_F1score_cutoff = args.training_f1Score_cutoff)
    barcodeSet.run()

if __name__ == "__main__":
    main()