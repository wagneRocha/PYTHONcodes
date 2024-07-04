#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_C_3(dirname, PDBchain, case, Cim1, N, CA, phi, i, pbt, AA, RR):
   v1 = CA   - N
   v2 = Cim1 - N
   
   d1 = np.linalg.norm(v1)
   d3 = np.linalg.norm(v2)
   
   dCim1CA = np.linalg.norm(CA - Cim1)

   d3d3 = d3*d3
   d3cosTh3 = 0.5*(d3d3 + d1*d1 - dCim1CA*dCim1CA)/d1
   d3sinTh3 = np.sqrt(d3d3 - d3cosTh3*d3cosTh3)

   # #############################
   # -----------------------------
   vecColsDouble = np.array([3, 6])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', i+1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   d2   = fixedSC[0]
   th1  = fixedSC[1]*np.pi/180
   Ed2  = EX[1]
   Eth1 = EX[4]*np.pi/180
   # ***************************************************************
   if((case[1] == '0') & (case[2] == '0')):
      abar = d2*np.cos(th1)
      dbar = d2*np.sin(th1)
   # ***************************************************************
   elif((case[1] == '0') & (case[2] == '1')):
      abar = Ed2*np.cos(th1)
      dbar = Ed2*np.sin(th1)
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '0')):
      abar = d2*np.cos(Eth1)
      dbar = d2*np.sin(Eth1)
   # ***************************************************************
   elif((case[1] == '1') & (case[2] == '1')):
      abar = Ed2*np.cos(Eth1)
      dbar = Ed2*np.sin(Eth1)
   # ***************************************************************
   else:
      print('ERROR in C{:3.0f}\n'.format(i+1))
   # ***************************************************************
   # #####################

   a = abar/d1
   b = d1
   c = d3cosTh3
   d = dbar/(d1*d3sinTh3)

   C = (CA - a*v1) + (d*b*v2 - d*c*v1)*np.cos(phi) + d*np.cross(v1,v2)*np.sin(phi)

   return C