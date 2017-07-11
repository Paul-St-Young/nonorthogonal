#!/usr/bin/env python
import os
import numpy as np

if __name__ == '__main__':

  gs = (16,16,16)

  import sys
  sys.path.insert(0,'../1_gs16/scripts')
  from step1_run_pyscf import step1_run_pyscf
  mf = step1_run_pyscf(gs)
  nao,nmo = mf.mo_coeff.shape

  detlist = np.loadtxt('../1_gs16/det_list.dat').view(complex)

  ikpt = 0 
  ispin = 0

  nfill = 4

  #fname = '../1_gs16/dets.h5'

  from check_orbitals import numpy_fft_state
  sys.path.insert(0,'../../15-pyscf2qmcpack/3_qmcpack')
  from show_fftbox import read_comp

  for idet in range(5):
    for istate in range(nfill):
      print( 'checking orbital %d in determinant %d' % (istate,idet) )
      from pyscf_orbital_routines import get_pyscf_psir
      # rotate MO orbitals
      det = detlist[idet].reshape(nmo,nmo)
      new_mocoeff = np.dot(mf.mo_coeff,det)
      new_movec   = new_mocoeff[:,istate]

      # get orbital from PySCF
      prgrid = get_pyscf_psir(new_movec,mf.cell)
      grid_shape = prgrid.shape

      iorb = idet*nfill + istate
      # get orbital from QMCPACK
      #print('getting orbitals from multideterminant QMCPACK')
      #qrgrid_arr = read_comp('fftbox%d_after.dat' % iorb)
      sdet_loc = './single/det%d/fftbox%d_after.dat' % (idet,istate)
      qrgrid_arr = read_comp(sdet_loc)
      qafter_fftbox = (qrgrid_arr[:,0]+1j*qrgrid_arr[:,1]).reshape(grid_shape)

      ## get orbital from hdf5
      #print('getting orbitals from hdf5 using python')
      #pgvecs,ppsig,pbefore_fftbox,pafter_fftbox = numpy_fft_state(fname
      #  ,ikpt,ispin,iorb)
      #qafter_fftbox = pafter_fftbox

      # compare orbitals
      import matplotlib.pyplot as plt

      # determinant
      #plt.matshow(det.imag)
      #plt.show()

      ## 2D
      #match_real = np.allclose( prgrid.real,qafter_fftbox.real )
      #match_imag = np.allclose( prgrid.imag,qafter_fftbox.imag )

      #if not match_real or not match_imag:
      #  fig,ax = plt.subplots(2,1)
      #  ax[0].plot( prgrid.real.flatten() )
      #  ax[0].plot( qafter_fftbox.real.flatten() )
      #  ax[1].plot( prgrid.imag.flatten() )
      #  ax[1].plot( qafter_fftbox.imag.flatten() )
      #  plt.show()
      ## end if

      # 3D
      if istate == -1: # visualize a particular state for each determinant, use -1 to disable
        pval = prgrid.real #np.absolute(prgrid)
        qval = qafter_fftbox.real #np.absolute(qafter_fftbox)
        print pval.min(),pval.max()
        print qval.min(),qval.max()
        sys.path.insert(0,'../../3-psir2hdf5/batch_convert')
        from convert_psir_to_psig import isosurf
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure()
        pax = fig.add_subplot(121,projection='3d',aspect=1)
        isosurf(pax, pval)
        qax = fig.add_subplot(122,projection='3d',aspect=1)
        isosurf(qax, qval)
        plt.show()
      # end if
