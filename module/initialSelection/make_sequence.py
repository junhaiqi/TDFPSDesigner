import numpy as np
import pandas as pd
import random
import os
import time
import datetime
import sys
import argparse
from seq_processing import *
#################################################################################################################################
#---------------------------------------------------Part1: Make index------------------------------------------------------#
def get_interval_size(length, quantity):
    total = 4**length
    interval_size = total//quantity
    return interval_size

def get_Index_list(length, quantity):
    Interval_size = get_interval_size(length, quantity)
    Index_list = list()

    for i in range(0, quantity):
        n = i * Interval_size
        num = random.randint(n, n + Interval_size)
        Index_list.append(num)

    return Index_list

#---------------------------------------------------Part2: Make sequence------------------------------------------------------#
def turn_base(key):
    selfComple_dict = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    base = selfComple_dict[key]
    return base

def calculate_key(number, power):
    a = 4**power
    Remainder = number%a
    key = (number - Remainder)/a
    key = int(key)%4
    return key

def make_sequence(index, length):
    sequence = ''
    for i in range(0, length):
        key_base = calculate_key(index, i)
        real_base = turn_base(key_base)
        sequence = real_base + sequence
    return sequence

def check_one(file_path, index, sequence):
    seq_list = get_seq_list(file_path)
    if seq_list[index] != sequence:
        print('Error!!!!!!!!\n')
        print('Error index: ' + str(index))

def check_all(file_path, mer):
    seq_list = get_seq_list(file_path)

    for i in range(0, 4**mer):
        seq = make_sequence(i, mer)
        if seq_list[i] != seq:
            print('Error!!!!!!!!')
            print('Error index: ' + str(i) + '\n')

def get_filter_seqlist(length, quantity, output_path):

    l = length
    quan = quantity

    index_list = get_Index_list(l, quan)

    seqs_list = list()
    ids_list = list()

    for index in index_list:
        ids_list.append(index)
        seq = make_sequence(index, l)
        seqs_list.append(seq)

    all_list = list(zip(ids_list, seqs_list))
    #print(len(all_list))

    outfile = output_path + str(l) + 'mer_filter_results.fasta'

    with open(outfile, 'w') as o:
        for i in range(len(all_list)):
            o.write('>' + str(i) + '|  original_index: ' + str(all_list[i][0]) + '\n')
            o.write(all_list[i][-1] + '\n')

# -----------------------------------------------Part3: main function -----------------------------------------------#
def initialization_parameters():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', action='store', dest='length', type=int, required=True, help='Sequence length')
    parser.add_argument('-q', action='store', dest='quantity', type=int, required=True, help='Filter size')
    parser.add_argument('-s', action='store', dest='seed', type=int, default=0, help='Random seed size')
    parser.add_argument('-o', action='store', dest='output', type=str, required=True, help='The output file')
    args = parser.parse_args()
    return args
###########################################################################################################################
# ----------------------------------------------------- main -----------------------------------------------------#

def mainFunction(length = 24, quantity = 100000, output = '', seed = 0):
    random.seed(seed)
    get_filter_seqlist(length, quantity, output)

def main():
    time_start = time.time()

    args = initialization_parameters()

    length = args.length
    quantity = args.quantity
    s = args.seed
    output = args.output

    random.seed(s)
    get_filter_seqlist(length, quantity, output)

    time_end = time.time() 
    time_sum = time_end - time_start 
    print('Program running time: {}'.format(time_sum))

if __name__ == "__main__":
    main()
