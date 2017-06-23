#!/usr/bin/env python
import numpy as np

if __name__ == '__main__':
  import sys 
  sys.path.insert(0,'../../15-pyscf2qmcpack/1_workflow/')
  from diamond_carbon import run_carbon

  # run pyscf and extract Kohn-Sham eigensystem
  # ================================================
  mf = run_carbon(verbose=3)
  print(mf.e_tot)

  from datetime import datetime
  from pyscf import tools, ao2mo
  from functools import reduce
  cell = mf.cell
  c = mf.mo_coeff
  nmo = c.shape[1]
  print(datetime.now())
  print('evaluating 1-electron integrals')
  h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
  print(datetime.now())
  print('evaluating 2-electron integrals')
  eri = mf.with_df.ao2mo(c) # 2-electron integrals
  print(datetime.now())
  print('restore 8-fold symmetry of the integral table')
  eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
  print(datetime.now())
  print('dumping the integral table')
  tools.fcidump.from_integrals('fcidump.dat',
    h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)
  print(datetime.now())

# end __main__
