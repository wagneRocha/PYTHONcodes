#!/usr/bin/python3

from datetime import date

## ------------------------------------------------------------------------------------------------------
def save_coordinates_pdb_format(filename, pdb_id, chain_name, res_name, res_seq, atom_name, X, method):

   n = len(X)

   with open(filename, 'w') as new_file:
      new_file.write('HEADER                                                {}   {}         \n'.format(date.today(), pdb_id))

      if method != 'PDB':
         new_file.write('TITLE     {} solution                                                          \n'.format(method))
      else:
         new_file.write('TITLE     PDB structure                                                         \n')

      new_file.write('REMARK   2                                                                      \n')
      new_file.write('REMARK   2 RESOLUTION.    2.00 ANGSTROMS.                                       \n')
      new_file.write('MODEL        1                                                                  \n')

      for k in range(n):
         atomType   = "ATOM  "
         serial     = k+1
         name       = atom_name[k]
         altLoc     = " "
         resName    = res_name[k]
         chainID    = chain_name
         resSeq     = int(res_seq[k])
         iCode      = " "
         x          = X[k,0]
         y          = X[k,1]
         z          = X[k,2]
         occupancy  = 1
         tempFactor = 1
         element    = atom_name[k][0]
         charge     = "  "

         new_file.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format(atomType, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge))
         #new_file.write('ATOM   {:4.0f} {} {} {}{:4.0f}    {:8.3f}{:8.3f}{:8.3f} {:5.2f} {:5.2f}           {}  \n'.format('ATOM', k+1, a_name, res_name[k], chain_name, res_seq[k], X[k,0], X[k,1], X[k,2], 1, 1, atom_name[k][0]))

      atomType = "TER   "
      serial   = k+1
      name     = "    "
      altLoc   = " "
      resName  = res_name[k]
      chainID  = chain_name
      resSeq   = int(res_seq[k])
      iCode    = " "

      new_file.write('{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}\n'.format(atomType, serial, name, altLoc, resName, chainID, resSeq, iCode))
      #new_file.write('TER    {:4.0f} {} {} {}{:4.0f}\n'.format(k+1, '    ', res_name[k], chain_name, res_seq[k]))