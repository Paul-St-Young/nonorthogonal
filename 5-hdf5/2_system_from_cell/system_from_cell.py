#!/usr/bin/env python
import numpy as np
from pyscf.pbc import gto as pbcgto
from pyscf.pbc import dft as pbcdft

def atom_text(elem,pos):
    assert len(elem) == len(pos)
    lines = []
    for iatom in range(len(elem)):
        mypos = pos[iatom]
        line = '%5s  %10.6f  %10.6f  %10.6f' % (elem[iatom],mypos[0],mypos[1],mypos[2])
        lines.append(line)
    atext = ';\n'.join(lines)
    return atext
# end def

def read_atoms(atext):
    # inverse of atom_text
    lines = atext.split(';')
    elem  = []
    pos   = []
    for line in lines:
        tokens = line.strip('\n').strip(' ').split()
        elem.append(tokens[0])
        pos.append( map(float,tokens[-3:]) )
    # end for line
    return elem,np.array(pos)
# end def read_atoms

def bcc2(alat,basis='cc-pVDZ',ke=20.):

    # define system
    axes  = alat*np.eye(3)
    elem  = ['H','H']
    upos  = np.array([[0,0,0],[0.5,0.5,0.5]])
    pos   = np.dot(upos,axes)
    atext = atom_text(elem,pos)

    # build simulation cell
    cell = pbcgto.Cell()
    cell.build(
        a    = axes,
        atom = atext,
        unit = 'B',
        basis = basis,
        ke_cutoff = ke
    )
    return cell
# end def bcc2

if __name__ == '__main__':

    alat = 3.77945227
    cell = bcc2(alat)

    abs_kpts = cell.make_kpts([2,2,2])        # kpoints in 2pi/alat units
    rkpts    = cell.get_scaled_kpts(abs_kpts) # kpoints in crystal units

    """
    # build simulation
    kmf = pbcdft.KRKS(cell,abs_kpts)
    kmf.xc = 'lda,lda'
    kmf.verbose = 3
    kmf.scf()
    """

    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    ref.system_from_cell(new,cell)

# end __main__
