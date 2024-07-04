#!/usr/bin/python3

## ------------------------------------------------------------------------------------------------------
def save_protein_chain_datfile(filename, chain_name, res_name, res_seq, atom_name, X):

   n = len(X)

   with open(filename, 'w') as new_file:
      for k in range(n):
         name       = atom_name[k]
         resName    = res_name[k]
         chainID    = chain_name
         resSeq     = int(res_seq[k])
         x          = X[k,0]
         y          = X[k,1]
         z          = X[k,2]

         new_file.write('{:^4s} {:3s} {:1s} {:4d} {:9.4f} {:9.4f} {:9.4f}\n'.format(name, resName, chainID, resSeq, x, y, z))