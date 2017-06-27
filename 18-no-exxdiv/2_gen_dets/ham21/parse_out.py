#!/usr/bin/env python

if __name__ == '__main__':
  import numpy as np
  import sys
  sys.path.insert(0,'..')
  from one_det import read_next_det

  from mmap import mmap
  ndet = 21 # !!!! hard-code number of determinants

  with open('ham21.out','r+') as f:
    mm = mmap(f.fileno(),0)
  # end with

  idx0,hmmt0 = read_next_det(mm,prefix="Hamiltonian")
  idx,hmmt = read_next_det(mm,prefix="Hamiltonian")
  idx,ovmt = read_next_det(mm,prefix="Overlap")

  nx = int(np.ceil(np.sqrt(len(hmmt))))
  assert nx==ndet

  hmat = hmmt.reshape(nx,nx,order='F')
  omat = ovmt.reshape(nx,nx,order='F')

  nuc = -12.66920340013197 # !!!! hard-code ion-ion repusion
  diag = (hmat.diagonal()/omat.diagonal()).real + nuc

  np.savetxt('diag.dat',diag-diag[0],fmt='%6.4f')
# end __main__
