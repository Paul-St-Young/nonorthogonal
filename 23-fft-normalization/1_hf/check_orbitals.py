#!/usr/bin/env python
import os
import numpy as np

def fft_grid_shape(gvecs):
  nmin  = gvecs.min(axis=0)
  nmax  = gvecs.max(axis=0)
  grid_shape = nmax-nmin+1
  return grid_shape
# end def fft_grid_shape

def numpy_fft_state(wf_fname,ikpt,ispin,istate):
  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read(wf_fname)

  # get size information
  gvecs = wf.get('gvectors')
  npw   = len(gvecs)
  grid_shape = fft_grid_shape(gvecs)

  # get plane-wave coefficients
  path = wf.state_path(ikpt,ispin,istate)
  psig_arr = wf.val( os.path.join(path,'psi_g') )
  psig = psig_arr[:,0] + 1j*psig_arr[:,1]
  assert len(psig) == npw

  # fft to real-space
  fftbox = np.zeros(grid_shape,dtype=complex)
  for ig in range(npw):
    fftbox[tuple(gvecs[ig])] = psig[ig]
  # end for

  before_fftbox = fftbox.copy()
  rgrid = np.fft.ifftn(fftbox) * np.prod(grid_shape)
  after_fftbox = rgrid

  # check previously written function
  rgrid0 = wf.psig_to_psir(gvecs,psig,grid_shape,1.0)
  assert np.allclose(rgrid,rgrid0)

  return gvecs,psig,before_fftbox,after_fftbox
# end def numpy_fft_state

if __name__ == '__main__':

  import sys
  sys.path.insert(0,'../../3-psir2hddf5/batch_convert')
  from convert_psir_to_psig import isosurf
  sys.path.insert(0,'../../15-pyscf2qmcpack/3_qmcpack')
  from show_fftbox import read_comp

  fname = 'pyscf2pwscf.h5'
  ikpt   = 0
  ispin  = 0
  istate = 0
  pw_fftbox   = False # show orbitals in reciprocal space
  real_fftbox = False # show orbitals in real space
  # !!!! hacked QMCPACK to dump plane-wave coefficients cG in psig_k{kpt}_s{state}.dat and fftbox in 'fftbox{state}_before.dat' and after FFT

  # real-space grid from python
  pgvecs,ppsig,pbefore_fftbox,pafter_fftbox = numpy_fft_state(fname
    ,ikpt,ispin,istate)
  npw   = len(pgvecs)
  grid_shape = fft_grid_shape(pgvecs)

  # check plane-wave coefficients
  qpsig_arr = read_comp('psig_k%d_s%d.dat'%(ikpt,istate))
  qpsig = qpsig_arr[:,0] + 1j*qpsig_arr[:,1]
  assert np.allclose(qpsig,ppsig)

  # check reciprocal-space fft box
  qbefore_fftbox_arr = read_comp('fftbox%d_before.dat'%istate)
  qbefore_fftbox = (qbefore_fftbox_arr[:,0]+1j*qbefore_fftbox_arr[:,1]).reshape(grid_shape)

  if pw_fftbox:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    pax  = fig.add_subplot(121,projection='3d',aspect=1)
    isosurf(pax, np.absolute(pbefore_fftbox) )
    qax  = fig.add_subplot(122,projection='3d',aspect=1)
    isosurf(qax, np.absolute(qbefore_fftbox) ,level=0.01)
    plt.show()
  # end if
  assert np.allclose(pbefore_fftbox,qbefore_fftbox)

  # check real-space fft box
  qrgrid_arr = read_comp('fftbox%d_after.dat' % istate)
  qafter_fftbox = (qrgrid_arr[:,0]+1j*qrgrid_arr[:,1]).reshape(grid_shape)

  if real_fftbox:
    import matplotlib.pyplot as plt
    fig = plt.figure()
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    pax  = fig.add_subplot(121,projection='3d',aspect=1)
    isosurf(pax, np.absolute(pafter_fftbox) )
    qax  = fig.add_subplot(122,projection='3d',aspect=1)
    isosurf(qax, np.absolute(qafter_fftbox) )
    plt.show()
  # end if
  assert np.allclose(pafter_fftbox.real,qafter_fftbox.real)
  assert np.allclose(pafter_fftbox.imag,qafter_fftbox.imag)

  import matplotlib.pyplot as plt
  fig,ax = plt.subplots(2,1)
  ax[0].plot( pafter_fftbox.real.flatten() )
  ax[0].plot( qafter_fftbox.real.flatten() )
  ax[1].plot( pafter_fftbox.imag.flatten() )
  ax[1].plot( qafter_fftbox.imag.flatten() )
  plt.show()

# end __main__
