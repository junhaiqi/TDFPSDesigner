#!/bin/sh

nvcc -std=c++11 -o ./bin/CalDTWDistMatrixMNSigroot CalDTWDistMatrixMNSigroot.cu -Xcompiler -fopenmp

nvcc -std=c++11 -o ./bin/CalDTWDistMatrixMN CalDTWDistMatrixMN.cu -Xcompiler -fopenmp

nvcc -std=c++11 -o ./bin/FpsCudaDTWThreshold FpsCudaDTWThreshold.cu -Xcompiler -fopenmp

echo compile  finished!
