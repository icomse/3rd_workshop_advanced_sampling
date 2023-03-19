#! /usr/bin/env python

import sys
import argparse
import numpy as np


help_text = '\
'

parser = argparse.ArgumentParser(description=help_text)
parser.add_argument('--colvar-file', dest='filename_colvar', required=True)
parser.add_argument('--output-file', dest='filename_output', required=True)
parser.add_argument('--time-min', dest='time_min', type=float, required=False)
parser.add_argument('--time-max', dest='time_max', type=float, required=False)
args = parser.parse_args()

filename_colvar = args.filename_colvar
print(' Input colvar file:  ',filename_colvar)
filename_output = args.filename_output
print(' Output colvar file: ',filename_output)
if(args.time_min):
    time_min = args.time_min
else:
    time_min = 0
if(args.time_max):
    time_max = args.time_max
else:
    time_max = 1.0e30

colvar_file = open(filename_colvar,'r')
output_file = open(filename_output,'w')
output_file.write(colvar_file.readline())
iline=0
for line in colvar_file:
    if (line[0:1] != '#' and line[0:1] != '@' and line.split() != []):
        time = float(line.split()[0])
        in_range = time_min <= time <= time_max
        if (in_range):
            output_file.write(line)
            iline = iline + 1
colvar_file.close()
output_file.close()
print(' Number of lines:    ',iline)
