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

  c = mf.mo_coeff
  nmo = c.shape[1]
  h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
  eri = mf.with_df.ao2mo(c) # 2-electron integrals
  eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
  tools.fcidump.from_integrals('fcidump.dat',
    h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)

# end __main__
