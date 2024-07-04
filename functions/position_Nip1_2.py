#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_Nip1_2(dirname, PDBchain, case, N, CA, C, psi, i, pbt, AA, RR):
   v1 = C - CA
   v2 = N - CA

   d1 = np.linalg.norm(v2)
   d2 = np.linalg.norm(v1)

   dNC = np.linalg.norm(N - C)

   d1d1 = d1*d1
   d1cosTh1 = 0.5*(d1d1 + d2*d2 - dNC*dNC)/d2
   d1sinTh1 = np.sqrt(d1d1 - d1cosTh1*d1cosTh1)

   # #############################
   # -----------------------------
   vecColsDouble = np.array([4, 7])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', i+1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/trigonometric_constants.dat'
   EcosORsinX = np.loadtxt(fname)
   # -----------------------------
   d3      = fixedSC[0]
   th2     = fixedSC[1]*np.pi/180
   Ed3     = EX[2]
   EcosTh2 = EcosORsinX[2]
   EsinTh2 = EcosORsinX[3]
   # ***************************************************************
   if((case[1] == '0') & (case[2] == '0')):
      abar = d3*np.cos(th2)
      dbar = d3*np.sin(th2)
   # ***************************************************************
   elif((case[1] == '0') & (case[2] == '1')):
      abar = Ed3*np.cos(th2)
      dbar = Ed3*np.sin(th2)
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '0')):
      abar = d3*EcosTh2
      dbar = d3*EsinTh2
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '1')):
      abar = Ed3*EcosTh2
      dbar = Ed3*EsinTh2
   # ***************************************************************
   else:
      print('ERROR in N{:3.0f}\n'.format(i+1))
   # ***************************************************************
   # #####################

   a = abar/d2
   b = d2
   c = d1cosTh1
   d = dbar/(d2*d1sinTh1)

   Nip1 = (C - a*v1) + (d*b*v2 - d*c*v1)*np.cos(psi) + d*np.cross(v1, v2)*np.sin(psi)

   return Nip1