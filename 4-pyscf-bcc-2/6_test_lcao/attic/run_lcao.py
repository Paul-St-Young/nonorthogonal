#!/usr/bin/env python
import os
import numpy as np
from pyscf.pbc import gto as pbcgto
from pyscf.pbc import dft as pbcdft

import sys
sys.path.insert(0,'../../../5-hdf5/2_system_from_cell/')
from system_from_cell import bcc2

if __name__ == '__main__':

    alat  = 3.77945227
    cell  = bcc2(alat)

    abs_kpts = cell.make_kpts([2,2,2])        # kpoints in 2pi/alat units

    # build simulation
    chk_file = 'pvdz.chk'
    kmf = pbcdft.KRKS(cell,abs_kpts)
    kmf.xc = 'lda,lda'
    kmf.verbose = 3
    if os.path.isfile(chk_file):
        kmf.init_guess = chk_file
    else:
        kmf.chkfile = chk_file
    # end if
    kmf.scf()

# end __main__
