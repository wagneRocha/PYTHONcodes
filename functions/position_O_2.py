#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_O_2(dirname, PDBchain, case, CA, C, Nip1, sigma, i, pbt, AA, RR):
   v1 = Nip1 - C 
   v2 = CA   - C

   d2 = np.linalg.norm(v2)
   d3 = np.linalg.norm(v1)
   
   dCANip1 = np.linalg.norm(Nip1 - CA)

   d2d2 = d2*d2
   d2cosTh2 = 0.5*(d2d2 + d3*d3 - dCANip1*dCANip1)/d3
   d2sinTh2 = np.sqrt(d2d2 - d2cosTh2*d2cosTh2)

   # #############################
   # -----------------------------
   vecColsDouble = np.array([5, 9])
   fname = f'{dirname}stereochemistry/{PDBchain}'
   [fixedSC, emptyVec] = readLineFileSringsAndDoubles(fname, '\t', i+1, vecColsDouble, np.array([]))
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/random_variables.dat'
   EX = np.loadtxt(fname)
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/trigonometric_constants.dat'
   EcosORsinX = np.loadtxt(fname)
   # -----------------------------
   d4        = fixedSC[0]
   th4       = fixedSC[1]*np.pi/180
   Ed4       = EX[3]
   EcosTh4   = EcosORsinX[6]
   EsinTh4   = EcosORsinX[7]
   EcosSigma = EcosORsinX[10]
   EsinSigma = EcosORsinX[11]
   # ***************************************************************
   if((case[0] == '0') & (case[1] == '0') & (case[2] == '0')):
      abar  = d4*np.cos(th4)
      dbarC = d4*np.sin(th4)*np.cos(sigma)
      dbarS = d4*np.sin(th4)*np.sin(sigma)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '0') & (case[2] == '1')):
      abar  = Ed4*np.cos(th4)
      dbarC = Ed4*np.sin(th4)*np.cos(sigma)
      dbarS = Ed4*np.sin(th4)*np.sin(sigma)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '1') & (case[2] == '0')):
      abar  = d4*EcosTh4
      dbarC = d4*EsinTh4*np.cos(sigma)
      dbarS = d4*EsinTh4*np.sin(sigma)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '1') & (case[2] == '1')):
      abar  = Ed4*EcosTh4
      dbarC = Ed4*EsinTh4*np.cos(sigma)
      dbarS = Ed4*EsinTh4*np.sin(sigma)
   # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '0')):
      abar  = d4*np.cos(th4)
      dbarC = d4*np.sin(th4)*EcosSigma
      dbarS = d4*np.sin(th4)*EsinSigma
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '1')):
      abar  = Ed4*np.cos(th4)
      dbarC = Ed4*np.sin(th4)*EcosSigma
      dbarS = Ed4*np.sin(th4)*EsinSigma
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '0')):
      abar  = d4*EcosTh4
      dbarC = d4*EsinTh4*EcosSigma
      dbarS = d4*EsinTh4*EsinSigma
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '1')):
      abar  = Ed4*EcosTh4
      dbarC = Ed4*EsinTh4*EcosSigma
      dbarS = Ed4*EsinTh4*EsinSigma
   # ***************************************************************
   else:
      print('ERROR in O{:3.0f}\n'.format(i+1))
   # ***************************************************************
   # #####################   

   a  = abar/d3
   b  = d3
   c  = d2cosTh2
   dc = dbarC/(d3*d2sinTh2)
   ds = dbarS/(d3*d2sinTh2)
   
   O = (C + a*v1) + (dc*b*v2 - dc*c*v1) + ds*np.cross(v1,v2)
   
   return O