#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_N1(dirname, PDBchain, case, AA):
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
   fname = f'{dirname}trans/{AA}/All/trigonometric_constants.dat'
   EcosORsinX = np.loadtxt(fname)
   # -----------------------------
   fname = f'{dirname}trans/{AA}/All/constants.dat'
   const = np.loadtxt(fname)
   # -----------------------------
   d1        = fixedSC[0]
   d2        = fixedSC[1]
   th1       = fixedSC[2]*np.pi/180
   Ed1       = EX[0]
   Ed2       = EX[1]
   Ecosth1   = EcosORsinX[0]
   Esinth1   = EcosORsinX[1]
   Ed1cosTh1 = const[0]
   Ed1sinTh1 = const[1]
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
      N[0] = -d1*Esinth1
      N[1] = d1*Ecosth1 - d2
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '1')):
      N[0] = -Ed1sinTh1
      N[1] = Ed1cosTh1 - Ed2
   # ***************************************************************
   else:
      print('ERROR in N1\n')
   # ***************************************************************
   # #############################

   return N