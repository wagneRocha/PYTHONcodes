#!/usr/bin/python3

import numpy as np
from read_the_lines_of_a_file import *
from split_a_row_in_columns import *

## ------------------------------------------------------------------------------------------------------
def readfileSringAndDoubles(fname, sepChar, vecColsDouble, vecColsChars):

   lines = read_the_lines_of_a_file(fname)

   m  = len(lines)
   nd = len(vecColsDouble)
   nc = len(vecColsChars)

   doubleData = np.zeros((m, nd))
   textData   = np.empty((m, nc), dtype=object)
   textData[:] = ''

   i = 0;
   for line in lines:
      columns = split_a_row_in_columns(line, sepChar)
      j = 0
      for col in vecColsDouble:
         doubleData[i,j] = columns[col-1]
         j = j + 1

      j = 0
      for col in vecColsChars:
         textData[i,j] = columns[col-1]
         j = j + 1

      i = i + 1

   return [doubleData, textData]