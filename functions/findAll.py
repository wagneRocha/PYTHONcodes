#!/usr/bin/python3

## ------------------------------------------------------------------------------------------------------
def findAll(vec, chr):
   vecChr = [pos for pos, char in enumerate(vec) if char == chr]
   return vecChr