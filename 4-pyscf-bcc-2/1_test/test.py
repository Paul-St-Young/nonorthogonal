#!/usr/bin/env python
import numpy as np
from pyscf.pbc import gto as pbcgto
import pyscf.pbc.tools.pyscf_ase as pyscf_ase

def build_cell(ase_atom,unit='B',ke=20.0,gsmax=None,basis='cc-pVDZ',
        pseudo=None,dimension=3):
    """ construct pyscf simulation cell from ase_atom """
    cell      = pbcgto.Cell()
    cell.unit = unit
    cell.atom = pyscf_ase.ase_atoms_to_pyscf(ase_atom)
    cell.a    = ase_atom.cell

    cell.basis  = basis
    cell.pseudo = pseudo
    cell.dimension = dimension

    if gsmax is not None: # different way to specify fft grid cutoff
        cell.gs = np.array([gsmax,gsmax,gsmax])
    else:
        cell.ke_cutoff = ke
    # end if

    cell.build()
    return cell
# end def build_cell

def run_dft(cell):
    from pyscf.pbc import dft as pbcdft
    sim = pbcdft.RKS(cell)
    #sim.xc = 'lda' # what is vwn?
    sim.verbose = 3
    sim.scf()
    return sim
# end def run_dft

if __name__ == '__main__':

    from ase.build import bulk
    ase_atom = bulk('H','bcc',a=3.77945227,cubic=True)
    cell = build_cell(ase_atom,basis='cc-pVDZ')

    test = run_dft(cell)

# end __main__
