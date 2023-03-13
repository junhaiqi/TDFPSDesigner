import os
import h5py
import numpy as np
from tqdm import tqdm
import sys
import time
import datetime
import argparse
############################################################################################################################

def get_seq_list(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    seq_list = filter(lambda x: x != '', lines)
    seq_list = filter(lambda x: '>' not in x, seq_list)
    a = list(seq_list)
    return a

def get_id_list(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    id_list = map(lambda x: x.split('|')[0][1:], lines)
    b = list(id_list)
    return b


def get_id_list2(file_name):
    '''
    Note: There must be no newline characters in each sequence.
    '''
    with open(file_name, 'r') as f:
        text = f.read()
        lines = text.splitlines()
    lines = filter(lambda x: '>' in x, lines)
    #id_list = map(lambda x: x.split('|')[0][1:], lines)
    #b = list(id_list)
    b = list(lines)
    return b

def groupmake(seq_file):
    seq_list = get_seq_list(seq_file)
    id_list = get_id_list(seq_file)
    #group_num = len(barcode_list)
    a = zip(seq_list, id_list)
    in_list = list(a)
    return in_list


