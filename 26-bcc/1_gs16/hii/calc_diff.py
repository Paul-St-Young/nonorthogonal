#!/usr/bin/env python
import os
import numpy as np

if __name__ == '__main__':
  ref0 = -7.76679008886357
  os.system('qmca -q ev --sac det*/*scalar.dat > results.dat')
  ndet = 5 # !!!! hard-code
  data = np.zeros([ndet,2])
  with open('results.dat','r') as f:
    f.readline() # skip header
    for line in f:
      tokens = line.split()
      idet = int( tokens[0].split('/')[0].replace('det','') )
      energy = float( tokens[3] )
      error  = float( tokens[5] )
      data[idet,0] = energy
      data[idet,1] = error
    # end for line
  # end with
  diff = np.zeros(data.shape)
  #diff[:,0] = data[:,0] - data[0,0]
  #diff[:,1] = np.sqrt( data[:,1]**2. + data[0,1]**2. )
  diff[:,0] = data[:,0] - ref0
  diff[:,1] = data[:,1]
  np.savetxt('diff.dat',diff,fmt='%6.4f')
# __main__
