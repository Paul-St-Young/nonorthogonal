#!/usr/bin/env python
import numpy as np
from pwscf_h5 import PwscfH5
import sys
sys.path.insert(0,'../../3-psir2hdf5/batch_convert/')
from convert_psir_to_psig import isosurf

def psig_to_psir(gvecs,psig,rgrid_shape,vol):
    """ contruct orbital given in planewave basis
     Inputs: 
      gvecs: gvectors in reciprocal lattice units i.e. integers
      psig: planewave coefficients, should have the same length as gvecs
      vol: simulation cell volume, used to normalized fft
     Output:
      rgrid: orbital on a real-space grid """
    assert len(gvecs) == len(psig)

    kgrid = np.zeros(rgrid_shape,dtype=complex)
    for igvec in range(len(gvecs)):
        kgrid[tuple(gvecs[igvec])] = psig[igvec]
    # end for
    rgrid = np.fft.ifftn(kgrid) * np.prod(rgrid_shape)/vol
    return rgrid
# end def psig_to_psir

def simulate_bcc2(alat=3.77945227,basis='cc-pVDZ',ke=20.,kgrid=[2,2,2],chkfile='pvdz.h5',from_scratch=False):
  import os
  import pyscf.pbc.gto as pbcgto
  import pyscf.pbc.dft as pbcdft

  # define system
  axes = alat*np.eye(3) 
  atom_text = 'H 0 0 0; H %f %f %f' % tuple(alat*np.array([0.5,0.5,0.5]))
  cell = pbcgto.Cell()
  cell.build(
    a         = axes,
    atom      = atom_text,
    unit      = 'B',
    basis     = basis,
    ke_cutoff = ke
  )

  # define simulation
  abs_kpts = cell.make_kpts(kgrid)
  kmf = pbcdft.KRKS(cell,abs_kpts)
  if os.path.isfile(chkfile) and (not from_scratch):
    from pyscf.pbc import lib
    kmf.__dict__.update(lib.chkfile.load(chkfile,'scf'))
  else:
    kmf.xc = 'lda,lda'
    kmf.verbose = 3
    kmf.chkfile = chkfile
    kmf.scf()
  # end if

  return cell,kmf
# end def simulate_bcc2

def get_psir(pwscf_h5_fname,ikpt=0,ispin=0,istate=0):
    wf = PwscfH5()
    wf.read(pwscf_h5_fname)

    # get lattice vectors
    axes = wf.get('axes')
    vol  = np.dot(np.cross(axes[0],axes[1]),axes[2])

    # get gvectors
    gvecs = wf.get('gvectors')

    # get eigenvector
    psig_arr = wf.psig(ikpt=ikpt,ispin=ispin,istate=istate)
    psig = psig_arr[:,0] + 1j*psig_arr[:,1]
    # determine real-space grid size (QMCPACK 3.0.0 convention)
    #  ref: QMCWaveFunctions/Experimental/EinsplineSetBuilder.cpp::ReadGvectors_ESHDF()
    mesh_factor = 1.0
    rgrid_shape = map(int, np.ceil(gvecs.max(axis=0)*4*mesh_factor) )
    rgrid_shape = [9,9,9] # !!!! override grid size

    psir = psig_to_psir(gvecs,psig,rgrid_shape,vol)
    return psir
# end def get_psir

def get_pyscf_psir(kmf,cell,ikpt=0,ispin=0,istate=0):
    # !!!! ispin is not used for KRKS
    # get molecular orbital
    from pyscf.pbc.dft import gen_grid,numint
    coords = gen_grid.gen_uniform_grids(cell)
    aoR = numint.eval_ao(cell,coords)
    rgrid_shape = cell.gs*2+1 # shape of real-space grid
    assert np.prod(rgrid_shape) == aoR.shape[0]

    moR = np.dot(aoR,kmf.mo_coeff[ikpt,:,istate])
    return moR.reshape(rgrid_shape)
# end def get_pyscf_psir

if __name__ == '__main__':

    # run DFT
    kgrid = [2,2,2]
    cell,kmf = simulate_bcc2(
      alat=3.77945227, basis='cc-pVDZ',chkfile='pvdz.h5',ke=20.,kgrid=kgrid)
    abs_kpts = cell.make_kpts(kgrid)
    rkpts0   = cell.get_scaled_kpts(abs_kpts)

    # define orbital to compare
    loc = {'ikpt':2,'ispin':0,'istate':5}

    psir0 = get_pyscf_psir(kmf,cell,**loc)
    val0  = np.absolute(psir0)
    print val0.min(),val0.max()

    # compare to pwscf.pwscf.h5
    fname ='../../5-hdf5/5_pyscf_hdf5/pwscf.pwscf.h5'
    wf = PwscfH5()
    wf.read(fname)
    rkpt_path = 'electrons/kpoint_%d/reduced_k'%loc['ikpt']
    rkpt  = wf.fp[rkpt_path].value
    assert np.allclose(rkpts0[loc['ikpt']],rkpt)
    psir  = get_psir(fname,**loc)
    val  = np.absolute(psir)
    print val.min(),val.max()

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax0 = fig.add_subplot(121,projection='3d',aspect=1)
    ax  = fig.add_subplot(122,projection='3d',aspect=1)
    isosurf(ax0,val0)
    isosurf(ax,val)
    plt.show()

# end __main__
