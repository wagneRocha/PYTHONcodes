#!/usr/bin/python3

import numpy as np

## ------------------------------------------------------------------------------------------------------
def read_the_lines_of_a_file(name):
   file = open(name,"r")
   lines = np.array(file.readlines())
   file.close()
   return lines