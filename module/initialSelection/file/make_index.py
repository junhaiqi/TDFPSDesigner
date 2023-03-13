import numpy as np
import pandas as pd
import random
import os
import time
import datetime
import sys
import argparse
from seq_processing import *
##########################################################################################################################

random.seed(0)
'''
a = random.randint(0,104)
b = random.randint(0,104)
print(a)
print(b)
'''
'''
c = 1048576

index_list = list()

for i in range(0,10000):
    n = i*100
    num = random.randint(n,n+100)
    index_list.append(num)

print(index_list)interval
print(len(index_list))'''

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


'''
d = get_interval_size(24,1000000)
print(d)'''

g = get_Index_list(24,1000000)
print(g)
print(len(g))