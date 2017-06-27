#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':
  ndet = 21 # !!!! hard-code
  data = np.zeros(ndet)
  with open('results.dat','r') as f:
    f.readline() # skip header
    for line in f:
      tokens = line.split()
      idet = int( tokens[0].split('/')[0].replace('det','') )
      energy = float( tokens[3] )
      data[idet] = energy
    # end for line
  # end with
  np.savetxt('diff.dat',data-data[0],fmt='%6.4f')

