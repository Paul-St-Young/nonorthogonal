def step2_dump_integral_table(mf):
  import os
  import numpy as np
  from datetime import datetime

  from pyscf import tools, ao2mo
  from functools import reduce
  cell = mf.cell
  c = mf.mo_coeff
  nmo = c.shape[1]

  print('evaluating 1-electron integrals')
  print(datetime.now())
  h1e = reduce(np.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
  print(datetime.now())

  print('evaluating 2-electron integrals')
  print(datetime.now())
  eri = mf.with_df.ao2mo(c) # 2-electron integrals
  print(datetime.now())

  print('restore 8-fold symmetry of the integral table')
  eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
  print('dumping the integral table')
  tools.fcidump.from_integrals('fcidump.dat',
    h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)
