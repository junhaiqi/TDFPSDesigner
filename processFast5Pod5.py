
import pod5
from ont_fast5_api.fast5_interface import get_fast5_file
import numpy as np
import os
import sys

def sig2text(sigList, outFile):
    file = open(outFile, 'w')
    for sig in sigList:
        file.write('%f\n'%sig)
    file.close()

def readFast5(fast5DirPath, outSigsDir):
    fast5Files = os.listdir(path=fast5DirPath)
    if not os.path.exists(outSigsDir):
        os.makedirs(outSigsDir)
    for fast5 in fast5Files:
        fast5Filepath = os.path.join(fast5DirPath, fast5)
        with get_fast5_file(fast5Filepath, mode = "r") as f5:
            for read in f5.get_reads():
                raw_data = read.get_raw_data()
                # readName = read.read_id.split('!')[1]
                readName = read.read_id
                sigFile = os.path.join(outSigsDir, '%s.txt'%readName)
                raw_data = list(raw_data)
                sig2text(raw_data, sigFile)

def readPOD5(POD5DirPath, outSigsDir):
    if not os.path.exists(outSigsDir):
        os.makedirs(outSigsDir)
    
    pod5Files = os.listdir(path=POD5DirPath)
    for pod5File in pod5Files:
        pod5FilePath = os.path.join(POD5DirPath, pod5File)
        with pod5.Reader(pod5FilePath) as reader:
            for read_record in reader.reads():
                readName = read_record.read_id
                sigFile = os.path.join(outSigsDir, '%s.txt'%readName)
                raw_data = list( read_record.signal )
                sig2text(raw_data, sigFile)



if __name__ == "__main__":
    if len( sys.argv ) != 4:
        print(f'Usage: python {sys.argv[0]} POD5/Fast5 InputDir OutputDir\n')
    else:
        print(f'Exacting signals from {sys.argv[2]} into {sys.argv[3]}...')
        if sys.argv[1] == 'POD5':
            readPOD5(POD5DirPath = sys.argv[2], outSigsDir = sys.argv[3])
        else:
            readFast5(fast5DirPath = sys.argv[2], outSigsDir = sys.argv[3])
        print(f'Exacting end!')