#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_N1_3(dirname, PDBchain, case, AA):
   N = np.zeros(3)
   # #############################
   # -----------------------------
   vecColsDouble = np.array([2, 3, 6])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', 1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}trans/{AA}/All/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   d1   = fixedSC[0]
   d2   = fixedSC[1]
   th1  = fixedSC[2]*np.pi/180
   Ed1  = EX[0]
   Ed2  = EX[1]
   Eth1 = EX[4]*np.pi/180
   # ***************************************************************
   if((case[1] == '0') & (case[2] == '0')):
      N[0] = -d1*np.sin(th1)
      N[1] = d1*np.cos(th1) - d2
   # ***************************************************************
   elif((case[1] == '0') & (case[2] == '1')):
      N[0] = -Ed1*np.sin(th1)
      N[1] = Ed1*np.cos(th1) - Ed2
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '0')):
      N[0] = -d1*np.sin(Eth1)
      N[1] = d1*np.cos(Eth1) - d2
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '1')):
      N[0] = -Ed1*np.sin(Eth1)
      N[1] = Ed1*np.cos(Eth1) - Ed2
   # ***************************************************************
   else:
      print('ERROR in N1\n')
   # ***************************************************************
   # #############################
   
   return N