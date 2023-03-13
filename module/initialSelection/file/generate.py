from itertools import product
import time
import os

###############################################################################################################################


p = "/data1/UMI_TEST/dex_test/filter/test_data/"


time_start = time.time()  # 记录开始时间

Base = ['A', 'C', 'G', 'T']
# new_barcode = []
c = 12
filename = p + 'seq_' + str(c) + 'mer.txt'
barcodelist = product(Base, repeat = c)
file_object = open(filename, 'w')
i = 0
for items in barcodelist:
    item = list(items)
    i = i + 1
    barcode = ''
    for line in item:
        barcode = barcode + line
    file_object.write( '>' + str(i-1) + '\n')
    file_object.write(barcode  + '\n')
    #print(barcode)
    #print(item)

file_object.close()

time_end = time.time()  # 记录结束时间
time_sum = time_end - time_start  # 计算的时间差为程序的执行时间，单位为秒/s
print(time_sum)