#!/usr/bin/python3

import numpy as np

## ------------------------------------------------------------------------------------------------------
def RMSDcalc(X0, Xr):
   
   n = X0.shape[0]
   X0bar = np.zeros((n,3))
   Xrbar = np.zeros((n,3))
   
   cmX0 = np.sum(X0, axis=0)/n
   cmXr = np.sum(Xr, axis=0)/n

   for i in range(n):
      X0bar[i,:] = X0[i,:] - cmX0
      Xrbar[i,:] = Xr[i,:] - cmXr

   U,S,Vt = np.linalg.svd(np.dot(np.transpose(Xrbar), X0bar))

   Q = np.dot(np.transpose(Vt), np.transpose(U))

   Xrbar = np.dot(Xrbar, np.transpose(Q))

   RMSD = np.linalg.norm(X0bar - Xrbar, ord='fro')/np.sqrt(n)

   return [X0bar, Xrbar, RMSD]