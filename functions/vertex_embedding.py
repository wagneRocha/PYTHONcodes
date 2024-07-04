#!/usr/bin/python3

import numpy as np

## ------------------------------------------------------------------------------------------------------
def vertex_embedding(x1, x2, x3, d34, th24, tau14):

   v1 = x1 - x2
   v2 = x3 - x2

   xhat = v2/np.linalg.norm(v2);
   
   z = np.cross(v2, v1)
   zhat = z/np.linalg.norm(z)

   yhat = np.cross(zhat, xhat)

   d34sinTh2 = d34*np.sin(th24)

   c1 = -d34*np.cos(th24)
   c2 =  d34sinTh2*np.cos(tau14)
   c3 =  d34sinTh2*np.sin(tau14)
   
   a = c1*xhat
   b = c2*yhat
   c = c3*zhat

   x4 = x3 + a + b + c

   return x4