#!/usr/bin/python3

import numpy as np
from read_the_lines_of_a_file import *
from split_a_row_in_columns import *

## ------------------------------------------------------------------------------------------------------
def readLineFileSringsAndDoubles(fname, sepChar, linei, vecColsDouble, vecColsChars):

   lines = read_the_lines_of_a_file(fname)
   
   nd = len(vecColsDouble)
   nc = len(vecColsChars)

   doubleData = np.zeros(nd)
   textData   = np.empty(nc, dtype=object)
   textData[:] = ''

   columns = split_a_row_in_columns(lines[linei-1], sepChar)
   j = 0
   for col in vecColsDouble:
      doubleData[j] = columns[col-1]
      j = j + 1

   j = 0
   for col in vecColsChars:
      textData[j] = columns[col-1]
      j = j + 1

   return [doubleData, textData]