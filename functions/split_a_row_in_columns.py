#!/usr/bin/python3

import numpy as np
from findAll import *

## ------------------------------------------------------------------------------------------------------
def split_a_row_in_columns(line, sepChar):
   vecSepChar = findAll(line, sepChar)
   barraN = line.find("\n")

   n = len(findAll(line, sepChar)) + 1
   col = line[0 : vecSepChar[0]].replace(" ", "")
   for j in range(1, n-1):
      col = np.append(col, line[vecSepChar[j-1]+1 : vecSepChar[j]].replace(" ", ""))
   col = np.append(col, line[vecSepChar[n-2]+1 : barraN].replace(" ", ""))

   return col