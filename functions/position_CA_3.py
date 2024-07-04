#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_CA_3(dirname, PDBchain, case, CAim1, Cim1, N, omega, i, pbt, AA, RR):
   v1 = N     - Cim1 
   v2 = CAim1 - Cim1

   d2 = np.linalg.norm(v2)
   d3 = np.linalg.norm(v1)
   
   dCAim1N = np.linalg.norm(N - CAim1)

   d2d2 = d2*d2
   d2cosTh2 = 0.5*(d2d2 + d3*d3 - dCAim1N*dCAim1N)/d3
   d2sinTh2 = np.sqrt(d2d2 - d2cosTh2*d2cosTh2)

   # #############################
   # -----------------------------
   vecColsDouble = np.array([2, 8])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', i+1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   d1     = fixedSC[0]
   th3    = fixedSC[1]*np.pi/180
   Ed1    = EX[0]
   Eth3   = EX[6]*np.pi/180
   Eomega = EX[8]*np.pi/180
   # ***************************************************************
   if((case[0] == '0') & (case[1] == '0') & (case[2] == '0')):
   #if((case[1] == '0') & (case[2] == '0')):
      abar  = d1*np.cos(th3)
      dbar0 = d1*np.sin(th3)*np.cos(omega)
      dbar1 = d1*np.sin(th3)*np.sin(omega)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '0') & (case[2] == '1')):
   #elif((case[1] == '0') & (case[2] == '1')):
      abar  = Ed1*np.cos(th3)
      dbar0 = Ed1*np.sin(th3)*np.cos(omega)
      dbar1 = Ed1*np.sin(th3)*np.sin(omega)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '1') & (case[2] == '0')):
   #elif((case[1] == '1') & (case[2] == '0')):
      abar  = d1*np.cos(Eth3)
      dbar0 = d1*np.sin(Eth3)*np.cos(omega)
      dbar1 = d1*np.sin(Eth3)*np.sin(omega)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '1') & (case[2] == '1')):
   #elif((case[1] == '1') & (case[2] == '1')):
      abar  = Ed1*np.cos(Eth3)
      dbar0 = Ed1*np.sin(Eth3)*np.cos(omega)
      dbar1 = Ed1*np.sin(Eth3)*np.sin(omega)
   # ***************************************************************
   #else:
   #   print('ERROR in CA{:3.0f}\n'.format(i+1))
   # ***************************************************************
   # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '0')):
      abar  = d1*np.cos(th3)
      dbar0 = d1*np.sin(th3)*np.cos(Eomega)
      dbar1 = d1*np.sin(th3)*np.sin(Eomega)
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '1')):
      abar  = Ed1*np.cos(th3)
      dbar0 = Ed1*np.sin(th3)*np.cos(Eomega)
      dbar1 = Ed1*np.sin(th3)*np.sin(Eomega)
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '0')):
      abar  = d1*np.cos(Eth3)
      dbar0 = d1*np.sin(Eth3)*np.cos(Eomega)
      dbar1 = d1*np.sin(Eth3)*np.sin(Eomega)
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '1')):
      abar  = Ed1*np.cos(Eth3)
      dbar0 = Ed1*np.sin(Eth3)*np.cos(Eomega)
      dbar1 = Ed1*np.sin(Eth3)*np.sin(Eomega)
   # ***************************************************************
   else:
      print('ERROR in CA{:3.0f}\n'.format(i+1))
   # *************************************************************** '''
   # #####################

   a  = abar/d3
   b  = d3
   c  = d2cosTh2
   dc = dbar0/(d3*d2sinTh2)
   ds = dbar1/(d3*d2sinTh2)

   CA = (N - a*v1) + (dc*b*v2 - dc*c*v1) + ds*np.cross(v1,v2)

   return CA