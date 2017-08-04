#!/usr/bin/env python
import os 
import numpy as np

def fft_grid_shape(gvecs):
  nmin  = gvecs.min(axis=0)
  nmax  = gvecs.max(axis=0)
  grid_shape = nmax-nmin+1
  return grid_shape.astype(int)
# end def fft_grid_shape

def numpy_fft_state(wf_fname,ikpt,ispin,istate):
  from pwscf_h5 import PwscfH5
  wf = PwscfH5()
  wf.read(wf_fname)

  # get size information
  gvecs = wf.get('gvectors').astype(int)
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

def overlap(rgrid1,rgrid2):
  num = (rgrid1.conj()*rgrid2).sum()
  den = np.sqrt( (rgrid1.conj()*rgrid1).sum()*(rgrid2.conj()*rgrid2).sum() )
  return num/den
# end def overlap

if __name__ == '__main__':

  pyscf_chkfile = '2_dz_gs40/a-uhf/gs40.h5'
  qmcpack_wf_h5 = '2_dz_gs40/a-uhf/pyscf2qmcpack.h5'
  ikpt   = 0
  nstate = 21
  nspin  = 2

  # construct mf from chkfile
  from pyscf.pbc import scf
  cell,scf_rec = scf.chkfile.load_scf(pyscf_chkfile)
  mf = scf.RHF(cell)
  mf.__dict__.update(scf_rec)

  # load pyscf atomic basis functions
  from pyscf.pbc.dft import gen_grid,numint
  coords = gen_grid.gen_uniform_grids(cell)
  aoR    = numint.eval_ao(cell,coords)

  # specify orbital to check
  for ispin in range(nspin):
    # build pyscf orbital
    moR = np.dot(aoR,mf.mo_coeff[ispin])

    for istate in range(nstate):

      print('checking (ikpt,ispin,istate) = (%d,%d,%d)' % (ikpt,ispin,istate))

      # load orbital from wavefunction file
      h5_gvecs,h5_psig,h5_before_box,h5_after_box = numpy_fft_state(qmcpack_wf_h5,ikpt,ispin,istate)
      npw = len(h5_gvecs)
      grid_shape = fft_grid_shape(h5_gvecs)
      assert np.allclose(grid_shape,h5_after_box.shape)
      gs = np.array(cell.gs)
      assert np.allclose(grid_shape,2*gs+1)

      py_orb =  moR[:,istate].reshape(grid_shape,order='C')
      h5_orb =  h5_after_box
      ov = overlap(py_orb,h5_orb).real
      print(ov)

      if ov < 0.99:
        py_orb = py_orb.real
        h5_orb = h5_orb.real

        import sys
        sys.path.insert(0,'../utils')
        from no_plotting import isosurf

        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()

        # show orbital from wavefunction hdf5
        h5_ax   = fig.add_subplot(1,2,1,projection='3d',aspect=1)
        h5_ax.set_title('qmcpack')
        h5_mesh = isosurf(h5_ax,h5_orb)

        # show orbital from pyscf
        py_ax   = fig.add_subplot(1,2,2,projection='3d',aspect=1)
        py_ax.set_title('pyscf')
        isosurf(py_ax,py_orb)

        plt.show()
      # end if
    # end for istate
  # end for ispin

# end __main__
