#!/usr/bin/env python
import os
import numpy as np
import pandas as pd
from pyscf import lib
from pyscf.pbc import gto, scf

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

def run_carbon():
    alat0 = 3.6 # angstrom?
    axes  = (np.ones((3,3))-np.eye(3))*alat0/2.0
    elem  = ['C','C']
    pos   = np.array([[0,0,0],[0.25,0.25,0.25]])*alat0
    atoms = atom_text(elem,pos)
    fname = 'szv.h5'

    cell = gto.Cell()
    cell.build(a=axes,atom=atoms,basis='gth-szv',pseudo='gth-pade'
      ,gs=np.array([5]*3), verbose=3)

    mf = scf.RHF(cell,exxdiv=None)
    if os.path.isfile(fname):
      mf.__dict__.update(lib.chkfile.load(fname,'scf'))
    else:
      mf.chkfile = fname
      mf.scf()
    # end if
    return mf
# end def

def main():
    mf = run_carbon()

    from pyscf.pbc.dft import gen_grid, numint
    coords = gen_grid.gen_uniform_grids(mf.cell)
    aoR = numint.eval_ao(mf.cell,coords)
    nao = aoR.shape[1]
    rgrid_shape = 2*mf.cell.gs+1
    assert np.prod(rgrid_shape)==aoR.shape[0]

    ci_coeff = np.loadtxt('../1_ref/ci_coeff.dat').view(complex)
    detlist  = np.loadtxt('../1_ref/detlist.dat').view(complex)
    ndet = len(detlist)
    assert len(ci_coeff) == ndet
    assert detlist.shape[1] == nao*nao

    eig_fname = 'eigsys.json'
    gfile     = 'gvectors.dat'

    # generate gvectors (no ke_cutoff)
    nx,ny,nz = mf.cell.gs
    from itertools import product
    int_gvecs = np.array([gvec for gvec in product(
      range(-nx,nx+1),range(-ny,ny+1),range(-nz,nz+1))])
    npw = len(int_gvecs)

    # turn detlist into a dataframe containing the eigensystem
    nfill = 4 # 4 filled orbitals
    ikpt=ispin=0 # only do this for RHF wavefunction at Gamma
    
    data = []
    for idet in range(ndet):
      det = detlist[idet].reshape(nao,nao)
      moR = np.dot(aoR,det)
      for iorb in range(nfill):
        rgrid = moR[:,iorb].reshape(rgrid_shape)
        moG   = np.fft.fftn(rgrid)/np.prod(rgrid_shape)*mf.cell.vol
        psig  = np.zeros([npw,2]) # store real & complex
        for igvec in range(npw):
          comp_val = moG[tuple(int_gvecs[igvec])]
          psig[igvec,:] = comp_val.real,comp_val.imag
        # end for igvec
        istate = idet*nfill+iorb
        entry  = {'ikpt':ikpt,'ispin':ispin,'istate':istate,
          'reduced_k':[0,0,0],'evalue':iorb,'evector':psig}
        data.append(entry)
      # end for iorb
    # end for idet

    df = pd.DataFrame(data)
    df.to_json(eig_fname)
    np.savetxt(gfile,int_gvecs)

    eig_df    = pd.read_json(eig_fname).set_index(
      ['ikpt','ispin','istate'],drop=True).sort_index()
    int_gvecs = np.loadtxt(gfile)
    # write pwscf.h5
    import h5py
    from pwscf_h5 import PwscfH5
    new = h5py.File('pwscf.pwscf.h5','w')
    ref = PwscfH5()
    nelecs = ref.system_from_cell(new,mf.cell,pseudized_charge={'C':2})
    ref.create_electrons_group(new,int_gvecs,eig_df,nelecs)

    # transfer version info.
    new.create_dataset('application/code',data=['pyscf'])
    new.create_dataset('application/version',data=['69d4b826c01950437f8e16663120942d5709c5e3'])
    new.create_dataset('format',data=['ES-HDF'])
    new.create_dataset('version',data=[2,1,0])
# end def main

if __name__ == '__main__':

    main()

# end __main__
