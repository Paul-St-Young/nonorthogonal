#!/usr/bin/env python

import h5py
import numpy
from functools import reduce
import lib
from pbc import gto, scf, dft
import tools,ao2mo

alat0 = 3.6

cell = gto.Cell()
cell.a = (numpy.ones((3,3))-numpy.eye(3))*alat0/2.0
cell.atom = (('C',0,0,0),('C',numpy.array([0.25,0.25,0.25])*alat0))
cell.basis = 'gth-szv'
cell.pseudo = 'gth-pade'
cell.gs = [5]*3  # use gs directly, instead of ke
cell.verbose = 4
cell.build()

mf = scf.RHF(cell,exxdiv=None)#,exxdiv=0)
#ehf = mf.kernel()
ehf = mf.scf()

c = mf.mo_coeff
nmo = c.shape[1]
h1e = reduce(numpy.dot, (c.T, mf.get_hcore(), c)) # 1-electron integrals
eri = mf.with_df.ao2mo(c) # 2-electron integrals
eri = ao2mo.restore('s8',eri,nmo) # use 8-fold symmetry of integrals of real orbitals [ij|kl] = [kj|il\ etc.
tools.fcidump.from_integrals('fcidump.dat', 
  h1e, eri, nmo, cell.nelectron,nuc=cell.energy_nuc(), ms=0)

