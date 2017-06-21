#!/usr/bin/env python
from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../6-verify_hdf5/3_check_no_ecut')
from no_ecut_psig_to_moR import isosurf
sys.path.insert(0,'../1_eigsys')
from carbon import ao_on_grid

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

def get_psir(pwscf_h5_fname,ikpt=0,ispin=0,istate=0):
    wf = PwscfH5()
    wf.read(pwscf_h5_fname)

    # get lattice vectors
    axes = wf.get('axes')
    vol  = np.dot(np.cross(axes[0],axes[1]),axes[2])

    # get gvectors
    gvecs = wf.get('gvectors').astype(int)

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

def get_pyscf_psir(mo_coeff,cell):
    # get molecular orbital
    aoR = ao_on_grid(cell)
    rgrid_shape = np.array(cell.gs)*2+1 # shape of real-space grid
    assert np.prod(rgrid_shape) == aoR.shape[0]

    moR = np.dot(aoR,mo_coeff)
    return moR.reshape(rgrid_shape)
# end def get_pyscf_psir

def check_orb_mf(mf,pwscf_h5_fname,ikpt=0,ispin=0,istate=0,thres=0.5):
    psir0 = get_pyscf_psir(mf.mo_coeff[:,istate],mf.cell)
    val0  = np.absolute(psir0)

    # compare to pwscf.pwscf.h5
    wf = PwscfH5()
    wf.read(pwscf_h5_fname)
    rkpt_path = 'electrons/kpoint_%d/reduced_k'%loc['ikpt']
    rkpt  = wf.fp[rkpt_path].value
    #assert np.allclose(rkpts0[loc['ikpt']],rkpt)
    psir  = get_psir(pwscf_h5_fname,**loc)
    val  = np.absolute(psir)

    # see if orbitals are the same
    percent_difference = np.linalg.norm(val0-val)/np.linalg.norm(val0)*100.
    if percent_difference > thres:
      #if istate > 3:
      #  assert 1==0
      #print val0.min(),val0.max()
      #print val.min(),val.max()
      print 'ikpt=%d istate=%d'%(ikpt,istate), percent_difference,'% different'
      import matplotlib.pyplot as plt
      fig = plt.figure()
      ax0 = fig.add_subplot(121,projection='3d',aspect=1)
      ax  = fig.add_subplot(122,projection='3d',aspect=1)
      isosurf(ax0,val0)
      isosurf(ax,val)
      plt.show()
    # end if
# end def check_orb

if __name__ == '__main__':
  sys.path.insert(0,'../1_eigsys')
  from carbon import run_carbon
  mf = run_carbon()

  ikpt = 0
  istate = 0
  for istate in range(len(mf.mo_energy)):
    print istate
    loc = {'ikpt':ikpt,'ispin':0,'istate':istate}
    check_orb_mf(mf,'pwscf.pwscf.h5',thres=1e-6,**loc)
  # end for istate
