#!/usr/bin/python3

import os, sys, math
import numpy as np

from all_functions import *

## ------------------------------------------------------------------------------------------------------
# main

# Open file with all pdb ids
dataset = np.genfromtxt('pdb_id.dat', delimiter='\t', dtype=str)
dirname = f'fixed_values/'

nStructures = dataset.shape[0]
RMSD = np.zeros((nStructures,8))

cases0 = np.array([0])
cases = np.array([0,1])

outputDir0 = f'output/'
if not os.path.exists(outputDir0):
   os.mkdir(outputDir0)

for k in range(nStructures):
   print(k+1)
   ## -------------------------------------------------------------------------------
   pdbID   = dataset[k, 0]
   chainID = dataset[k, 1]
   PDBchain = f'{pdbID}_{chainID}.dat'

   #print(PDBchain)
   ## -------------------------------------------------------------------------------
   fname0 = f'dat_structures/{PDBchain}'
   [Xdouble, Xchar] = readfileSringAndDoubles(fname0, "\t", np.array([2, 4, 5, 6]), np.array([1, 3]))

   X = Xdouble[:,1:4]
   n = X.shape[0]

   ind = np.arange(0, n, 4)
   nRes = n//4

   RES = Xchar[:,0]
   res = RES[ind]
   ## -------------------------------------------------------------------------------
   fname1 = f'{dirname}stereochemistry/{PDBchain}'
   cs = 0
   for RRmethod in cases0:
      [torsionAngles, RM_HC] = readfileSringAndDoubles(fname1, "\t", np.array([10, 11, 12, 13]), np.array([1]))
      if(RRmethod == 0):
         RM_HC = '*'*nRes
      omega = torsionAngles[:,0]*np.pi/180;
      phi   = torsionAngles[:,1]*np.pi/180;
      psi   = torsionAngles[:,2]*np.pi/180;
      sigma = torsionAngles[:,3]*np.pi/180;
      ## -------------------------------------------------------------------------------
      for case1 in cases:
         for case2 in cases:
            for case3 in cases:
               case = str(case1)+str(case2)+str(case3)
               ## -----------------------
               outputDir0 = f'output/pdb/'
               if not os.path.exists(outputDir0):
                  os.mkdir(outputDir0)
               outputDir0 = f'output/pdb/{case}/'
               if not os.path.exists(outputDir0):
                  os.mkdir(outputDir0)
               outputDir1 = f'output/dat/'
               if not os.path.exists(outputDir1):
                  os.mkdir(outputDir1)
               outputDir1 = f'output/dat/{case}/'
               if not os.path.exists(outputDir1):
                  os.mkdir(outputDir1)
               ## -----------------------
               N  = np.zeros((nRes, 3))
               CA = np.zeros((nRes, 3))
               C  = np.zeros((nRes, 3))
               O  = np.zeros((nRes, 3))
               ## -------------------------------------------------------------------------------
               #print(1)
               Npos = position_N1_2(dirname, PDBchain, case, res[0])
               #Npos = X[0,:]
               ## -----------------------
               CApos = position_CA1_2(dirname, PDBchain, case, res[0])
               #CApos = X[1,:]
               ## -----------------------
               Cpos = C[0,:]
               #Cpos = X[2,:]
               ## -----------------------
               ## -----------------------
               ## -----------------------
               RR = RM_HC[0][0]
               if(RR == "*"):
                  RR = "All"
               Np1pos = position_Nip1_2(dirname, PDBchain, case, Npos, CApos, Cpos, psi[0], 0, 'trans', res[0], RR)
               ## -----------------------
               Opos   = position_O_2(dirname, PDBchain, case, CApos, Cpos, Np1pos, sigma[0], 0, 'trans', res[0], RR)
               ## -----------------------
               N[0,:]  = Npos
               CA[0,:] = CApos
               C[0,:]  = Cpos
               N[1,:]  = Np1pos
               O[0,:]  = Opos
               ## -------------------------------------------------------------------------------
               ## -------------------------------------------------------------------------------
               CAim1pos = CApos
               Cim1pos  = Cpos
               Npos     = Np1pos
               ## -------------------------------------------------------------------------------
               for i in range(1, nRes-1):
                  #print(i+1)
                  ## -------------------------------------------------------------------------------
                  if(abs(omega[i]) >= np.pi/2):
                     pbtype = "trans"
                  else:
                     pbtype = "cis"

                  RR = RM_HC[i][0]
                  if(RR == "*"):
                     RR = "All"
                  ## -------------------------------------------------------------------------------
                  CApos  = position_CA_2(dirname, PDBchain, case, CAim1pos, Cim1pos, Npos, omega[i], i, pbtype, res[i], RR)
                  ## -----------------------
                  Cpos   = position_C_2(dirname, PDBchain, case, Cim1pos, Npos, CApos, phi[i], i, pbtype, res[i], RR)
                  ## -----------------------
                  Np1pos = position_Nip1_2(dirname, PDBchain, case, Npos, CApos, Cpos, psi[i], i, pbtype, res[i], RR)
                  ## -----------------------
                  Opos   = position_O_2(dirname, PDBchain, case, CApos, Cpos, Np1pos, sigma[i], i, pbtype, res[i], RR)
                  ## -----------------------
                  CA[i,:]  = CApos
                  C[i,:]   = Cpos
                  N[i+1,:] = Np1pos
                  O[i,:]   = Opos
                  ## -------------------------------------------------------------------------------
                  ## -------------------------------------------------------------------------------
                  CAim1pos = CApos
                  Cim1pos  = Cpos
                  Npos     = Np1pos
               ## -------------------------------------------------------------------------------
               #print(nRes)
               ## -------------------------------------------------------------------------------
               CApos = position_CA_2(dirname, PDBchain, case, CAim1pos, Cim1pos, Npos, omega[nRes-1], nRes-1, 'trans', res[nRes-1], RR)
               ## -----------------------
               Cpos  = position_C_2(dirname, PDBchain, case, Cim1pos, Npos, CApos, phi[nRes-1], nRes-1, 'trans', res[nRes-1], RR)
               ## -----------------------
               CA[nRes-1,:] = CApos
               C[nRes-1,:]  = Cpos
               ## -------------------------------------------------------------------------------
               ## -------------------------------------------------------------------------------
               Xc = np.zeros_like(X)
                              
               Xc[ind+0,:] = N
               Xc[ind+1,:] = CA
               Xc[ind+2,:] = C
               Xc[ind+3,:] = O
               ## -------------------------------------------------------------------------------
               #exc = np.array(range(5, n+1)) - 1
               #X  = np.delete(X, exc, axis=0)
               #Xc = np.delete(Xc, exc, axis=0)
               X0 = np.delete(X, n-1, axis=0)
               Xc = np.delete(Xc, n-1, axis=0)
               ## -------------------------------------------------------------------------------
               Xbar, Xcbar, RMSD[k,cs] = RMSDcalc(X0, Xc)
               
               #print(RMSD[k,cs])

               fname_k = f'{outputDir0}/{pdbID}_{chainID}_X_{case}_2.pdb'
               save_coordinates_pdb_format(fname_k, pdbID, chainID, Xchar[:,0], Xdouble[:,0], Xchar[:,1], Xc, case)
               fname_k = f'{outputDir1}{pdbID}_{chainID}_X_{case}_2.dat'
               np.savetxt(fname_k, Xc)

               cs+= 1

   #input("Press Enter to continue...")

np.savetxt('rmsd_2.dat', RMSD)