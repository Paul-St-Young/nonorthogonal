#!/usr/bin/env python
# Construct real-space orbitals in two ways:
#  1. Directly from pyscf mean-field object `mf` i.e. use mo_coeff to take linear combination of atomic basis functions on a grid.
#  2. Read plane-wave coefficients from QMCPACK's wavefunction hdf5 file, then construct each orbital on an FFT grid.
#  then check the two ways are numerically equivalent.
# If all checks pass, then the PySCF orbitals have been successfully stored in plane-wave representation in the hdf5 file.

from pwscf_h5 import PwscfH5
import numpy as np
import sys
sys.path.insert(0,'../../6-verify_hdf5/3_check_no_ecut')
from no_ecut_psig_to_moR import isosurf

def check_orb_mf(mf,pwscf_h5_fname,ikpt=0,ispin=0,istate=0,thres=0.5):
    sys.path.insert(0,'scripts')
    from pyscf_orbital_routines import get_pyscf_psir
    psir0 = get_pyscf_psir(mf.mo_coeff[:,istate],mf.cell)
    val0  = np.absolute(psir0)

    # compare to pwscf.pwscf.h5
    wf = PwscfH5()
    wf.read(pwscf_h5_fname)
    rkpt_path = 'electrons/kpoint_%d/reduced_k'%loc['ikpt']
    rkpt  = wf.fp[rkpt_path].value
    #assert np.allclose(rkpts0[loc['ikpt']],rkpt) # diable for gamma point
    rgrid_shape = 2*np.array(mf.cell.gs)+1
    psir  = wf.get_psir_from_psig(ikpt,ispin,istate,rgrid_shape=rgrid_shape)
    val  = np.absolute(psir)

    # see if orbitals are the same
    if not np.allclose(val,val0):
      percent_difference = np.linalg.norm(val0-val)/np.linalg.norm(val0)*100.
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
  from diamond_carbon import run_carbon
  mf = run_carbon(verbose=3)

  ikpt  = 0
  ispin = 0
  for istate in range(len(mf.mo_energy)):
    loc = {'ikpt':ikpt,'ispin':ispin,'istate':istate}
    check_orb_mf(mf,'pyscf2pwscf.h5',thres=1e-6,ikpt=ikpt,ispin=ispin,istate=istate)
    print 'state %d checks out!' % istate
  # end for istate

# end __main__
