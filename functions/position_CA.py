#!/usr/bin/python3

import numpy as np
from readLineFileSringsAndDoubles import *

## ------------------------------------------------------------------------------------------------------
def position_CA(dirname, PDBchain, case, CAim1, Cim1, N, omega, i, pbt, AA, RR):
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
   fname = f'{dirname}{pbt}/{AA}/{RR}/trigonometric_constants.dat'
   EcosORsinX = np.loadtxt(fname)
   # -----------------------------
   fname = f'{dirname}{pbt}/{AA}/{RR}/constants.dat'
   const = np.loadtxt(fname)
   # -----------------------------
   d1                = fixedSC[0]
   th3               = fixedSC[1]*np.pi/180
   Ed1               = EX[0]
   EcosTh3           = EcosORsinX[4]
   EsinTh3           = EcosORsinX[5]
   EcosOmega         = EcosORsinX[8]
   EsinOmega         = EcosORsinX[9]
   Ed1cosTh3         = const[4]
   Ed1sinTh3         = const[5]
   Ed1cosOmega       = const[10]
   Ed1sinOmega       = const[11]
   EsinTh3cosOmega   = const[14]
   EsinTh3sinOmega   = const[15]
   Ed1sinTh3cosOmega = const[18]
   Ed1sinTh3sinOmega = const[19]
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
      abar  = d1*EcosTh3
      dbar0 = d1*EsinTh3*np.cos(omega)
      dbar1 = d1*EsinTh3*np.sin(omega)
   # ***************************************************************
   elif((case[0] == '0') & (case[1] == '1') & (case[2] == '1')):
   #elif((case[1] == '1') & (case[2] == '1')):
      abar  = Ed1cosTh3
      dbar0 = Ed1sinTh3*np.cos(omega)
      dbar1 = Ed1sinTh3*np.sin(omega)
   # ***************************************************************
   #else:
   #   print('ERROR in CA{:3.0f}\n'.format(i+1))
   # ***************************************************************
   # xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '0')):
      abar  = d1*np.cos(th3)
      dbar0 = d1*np.sin(th3)*EcosOmega
      dbar1 = d1*np.sin(th3)*EsinOmega
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '0') & (case[2] == '1')):
      abar  = Ed1*np.cos(th3)
      dbar0 = Ed1cosOmega*np.sin(th3)
      dbar1 = Ed1sinOmega*np.sin(th3)
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '0')):
      abar  = d1*EcosTh3
      dbar0 = d1*EsinTh3cosOmega
      dbar1 = d1*EsinTh3sinOmega
   # ***************************************************************
   elif((case[0] == '1') & (case[1] == '1') & (case[2] == '1')):
      abar  = Ed1cosTh3
      dbar0 = Ed1sinTh3cosOmega
      dbar1 = Ed1sinTh3sinOmega
   # ***************************************************************
   else:
      print('ERROR in CA{:3.0f}\n'.format(i+1))
   # ***************************************************************'''
   # #####################

   a  = abar/d3
   b  = d3
   c  = d2cosTh2
   dc = dbar0/(d3*d2sinTh2)
   ds = dbar1/(d3*d2sinTh2)

   CA = (N - a*v1) + (dc*b*v2 - dc*c*v1) + ds*np.cross(v1,v2)

   return CA