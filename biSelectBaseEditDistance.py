from edlib import align
from multiprocessing import Pool
import argparse

def ourAlign(seq1, seq2):
    res = align(seq1, seq2, mode = "NW")
    return res["editDistance"]

def writeSeqs2Fasta(seqList, filePath):
    file = open(filePath, 'w')
    t = 0
    for seq in seqList:
        file.write('>%d\n'%t)
        file.write('%s\n'%seq)
        t += 1
    file.close()

def editDistBetweenSeqAndSeqList(seq, seqList, threadNum):
    pool = Pool(threadNum)
    args = [(seq, item) for item in seqList]
    distList = list(pool.starmap(ourAlign, args))
    pool.close()
    pool.join()
    return min(distList)

def fpsWithEditDist(seqList, editDistThreshold, threadNum):

    selectedSeqList = []
    initSeq = seqList[0]
    seqList.pop(0)
    selectedSeqList.append(initSeq)

    for seq in seqList:
        minEditDist = editDistBetweenSeqAndSeqList(seq, selectedSeqList, threadNum)
        if minEditDist > (editDistThreshold - 1):
            selectedSeqList.append(seq)

    return selectedSeqList

def writeSeqs2Fasta(seqList, filePath):
    file = open(filePath, 'w')
    t = 0
    for seq in seqList:
        file.write('>%d\n'%t)
        file.write('%s\n'%seq)
        t += 1
    file.close()

def getSeqsFromFasta(file):
    _file = open(file)
    seqList = []
    for line in _file:
        if line[0] != '>':
            seq = line.strip('\n')
            seqList.append(seq)
    _file.close()
    return seqList

def get_parameters():
    """This script generates noiseless nanopore signals based on DNA sequences"""
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta-file', type = str, required = True,
                        help='It is the fasta file that contains short sequence fragments.')

    parser.add_argument('--edit-dist', type = int, required = True,
                        help='It is an edit distance threshold.')

    parser.add_argument('--out-file', type = str, required = True,
                        help='It is a fasta file to store the selected barcode sequence.')

    parser.add_argument('--thread-num', type=int, required=False, default = 2,
                        help='It is the number of threads used to execute the task.')

    args = parser.parse_args()

    return args


def main():
    args = get_parameters()
    seqList = getSeqsFromFasta(args.fasta_file)
    barList = fpsWithEditDist(seqList = seqList, editDistThreshold = args.edit_dist, threadNum = args.thread_num)
    writeSeqs2Fasta(barList, args.out_file)

if __name__ == "__main__":
   main()

    