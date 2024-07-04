#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_CA1_3(dirname, PDBchain, case, AA):
   CA = np.zeros(3)
   # #############################
   # -----------------------------
   vecColsDouble = np.array([3])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', 1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}trans/{AA}/All/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   d2  = fixedSC[0]
   Ed2 = EX[1]
   # ***************************************************************
   if(case[2] == '0'):
      CA[1] = -d2
   # ***************************************************************
   elif(case[2] == '1'):
      CA[1] = -Ed2
   # ***************************************************************
   else:
      print('ERROR in CA1\n')
   # ***************************************************************
   # #############################
   
   return CA